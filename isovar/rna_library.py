"""Immutable original-RNA indexes shared by reconstruction workflows.

The BAM remains authoritative. The index stores virtual offsets and conservative
placement bounds, not a regional replacement for the original library. Exact
CIGAR and sequence evidence is assessed by the existing reconstruction layers.
"""

from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from contextlib import ExitStack, closing, contextmanager
from dataclasses import dataclass
import errno
from hashlib import sha256
import json
from os import fsdecode
from pathlib import Path
import sqlite3
from tempfile import TemporaryDirectory
from typing import Any

import pysam
from varcode import ContigResolver

from .audit_bundle import software_identity
from .read_collector import ReadCollector
from .read_identity import segment_identity
from .rna_evidence import canonical_json, content_identifier


Fragment = tuple[str, str]
Segment = tuple[str, str, int]
Region = tuple[str, int, int]
_SCHEMA = "isovar.rna_library.v1"


def _digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _file_state(path: Path) -> tuple[int, int, int, int, int]:
    stat = path.stat()
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns


def _read_at(reader: pysam.AlignmentFile, offset: int) -> pysam.AlignedSegment:
    """Seek only when the next record is not already in the decompressed block."""
    current = reader.tell()
    if current != offset:
        if current < offset and current >> 16 == offset >> 16:
            while reader.tell() < offset:
                next(reader)
        else:
            reader.seek(offset)
    return next(reader)


def _eligibility_policy(collector: ReadCollector, read_filter_identity: str | None) -> dict[str, Any]:
    """Pin only original-record eligibility; event/path settings remain separate."""
    custom_policy = getattr(collector.alignment_filter_reason, "__func__", None) is not ReadCollector.alignment_filter_reason
    custom_filter = collector.read_filter is not None
    if (custom_filter or custom_policy) and not read_filter_identity:
        raise ValueError("A reusable RNA index with custom eligibility requires an explicit read_filter_identity")
    if not (custom_filter or custom_policy) and read_filter_identity is not None:
        raise ValueError("read_filter_identity requires a custom read_filter or alignment_filter_reason")
    return dict(use_secondary_alignments=collector.use_secondary_alignments,
                use_duplicate_reads=collector.use_duplicate_reads,
                min_mapping_quality=collector.min_mapping_quality,
                use_reads_without_base_qualities=collector.use_reads_without_base_qualities,
                custom_alignment_policy=custom_policy, custom_read_filter=custom_filter,
                read_filter_identity=read_filter_identity)


def _uses_custom_eligibility(collector: ReadCollector) -> bool:
    return (collector.read_filter is not None or
            getattr(collector.alignment_filter_reason, "__func__", None) is not ReadCollector.alignment_filter_reason)


def _connect(path: Path, *, readonly: bool = False) -> sqlite3.Connection:
    db = sqlite3.connect(path.as_uri() + "?mode=ro&immutable=1", uri=True) if readonly else sqlite3.connect(path)
    db.execute("PRAGMA temp_store=FILE")
    db.execute("PRAGMA cache_size=-32768")
    return db


def _index_records(source: Path, db: sqlite3.Connection, collector: ReadCollector) -> dict[str, Any]:
    """Index original offsets once; retain excluded/duplicate rows for accounting."""
    db.execute("CREATE TABLE records (offset INTEGER PRIMARY KEY, rg TEXT NOT NULL, name TEXT NOT NULL, "
               "segment INTEGER NOT NULL, excluded INTEGER NOT NULL, duplicate_of INTEGER, rid BLOB)")
    # Floating R-tree bounds round outward. This is a conservative acquisition
    # index: final CIGAR checks remain exact, including skipped introns.
    db.execute("CREATE VIRTUAL TABLE placements USING rtree(offset, ref_min, ref_max, start_min, end_max)")
    excluded: Counter[str] = Counter()
    scanned = 0
    with pysam.AlignmentFile(source, "rb") as reader:
        if not reader.is_bam:
            raise ValueError("A reusable RNA library index requires a seekable BAM")
        header = reader.header.to_dict()
        while True:
            offset = reader.tell()
            read = next(reader, None)
            if read is None:
                break
            scanned += 1
            reason = "unnamed_record" if read.query_name is None else collector.alignment_filter_reason(read)
            if reason:
                excluded[reason] += 1
            group, name, segment = segment_identity(read)
            if not isinstance(group, str):
                raise ValueError("SAM RG tags must be strings")
            db.execute("INSERT INTO records VALUES (?, ?, ?, ?, ?, NULL, NULL)",
                       (offset, group, name or "", segment, int(reason is not None)))
            if not read.is_unmapped and read.reference_start >= 0 and read.reference_end is not None:
                db.execute("INSERT INTO placements VALUES (?, ?, ?, ?, ?)",
                           (offset, read.reference_id, read.reference_id, read.reference_start, read.reference_end))
            if scanned % 100000 == 0:
                db.commit()
    db.execute("CREATE INDEX record_templates ON records (rg, name)")
    db.execute("CREATE TEMP TABLE repeated (rg TEXT, name TEXT, PRIMARY KEY (rg, name)) WITHOUT ROWID")
    db.execute("INSERT INTO repeated SELECT rg, name FROM records GROUP BY rg, name HAVING count(*) > 1")
    db.execute("CREATE TEMP TABLE payloads (rid BLOB PRIMARY KEY, offset INTEGER UNIQUE)")
    duplicates = 0
    with pysam.AlignmentFile(source, "rb") as reader, pysam.AlignmentFile(source, "rb") as previous:
        for offset, in db.execute("SELECT r.offset FROM records r JOIN repeated USING (rg, name) ORDER BY r.offset"):
            sam = _read_at(reader, offset).to_string()
            rid = sha256(sam.encode()).digest()
            inserted = db.execute("INSERT OR IGNORE INTO payloads VALUES (?, ?)", (rid, offset)).rowcount
            original = None
            if not inserted:
                original, = db.execute("SELECT offset FROM payloads WHERE rid=?", (rid,)).fetchone()
                if _read_at(previous, original).to_string() != sam:
                    raise ValueError("Conflicting original SAM payloads share a record digest")
                duplicates += 1
            db.execute("UPDATE records SET rid=?, duplicate_of=? WHERE offset=?", (rid, original, offset))
    db.execute("CREATE TABLE templates (rg TEXT, name TEXT, first_offset INTEGER, eligible_segments INTEGER, "
               "PRIMARY KEY (rg, name)) WITHOUT ROWID")
    db.execute("INSERT INTO templates SELECT rg, name, min(offset), "
               "sum(DISTINCT CASE WHEN excluded=0 AND duplicate_of IS NULL THEN (1 << (segment >> 6)) ELSE 0 END) "
               "FROM records GROUP BY rg, name")
    db.execute("CREATE INDEX template_traversal ON templates (first_offset)")
    fragments, = db.execute("SELECT count(*) FROM templates").fetchone()
    segments, = db.execute("SELECT count(*) FROM (SELECT rg, name, segment FROM records "
                           "WHERE excluded=0 AND duplicate_of IS NULL GROUP BY rg, name, segment)").fetchone()
    db.commit()
    if db.execute("PRAGMA quick_check").fetchall() != [("ok",)]:
        raise ValueError("RNA library index integrity check failed")
    if db.execute("SELECT rtreecheck('placements')").fetchone() != ("ok",):
        raise ValueError("RNA library placement index integrity check failed")
    return dict(header=header, records_scanned=scanned, unique_records=scanned - duplicates,
                fragments_scanned=fragments, eligible_segments=segments, excluded_records=dict(sorted(excluded.items())))


@dataclass(frozen=True)
class PreparedRnaLibrary:
    """Paths and identity of a validated, immutable index for one original BAM."""

    source: Path
    directory: Path
    _identity_json: str
    _statistics_json: str
    source_state: tuple[int, int, int, int, int]
    database_state: tuple[int, int, int, int, int]

    @property
    def identity(self) -> dict[str, Any]:
        """Return a detached identity snapshot; caller edits cannot alter the index."""
        return json.loads(self._identity_json)

    @property
    def statistics(self) -> dict[str, Any]:
        """Return detached original-record counts and the source SAM header."""
        return json.loads(self._statistics_json)

    def open(self) -> "RnaLibraryIndex":
        """Open an independent reader; use it as a context manager."""
        return RnaLibraryIndex(self)


def prepare_rna_library(source: str | Path, cache_dir: str | Path, *, read_collector: ReadCollector | None = None,
                       read_filter_identity: str | None = None, reference_name: str | None = None) -> PreparedRnaLibrary:
    """Validate or atomically prepare a content-addressed original-RNA index.

    Parameters
    ----------
    source : str or pathlib.Path
        Seekable original BAM. Records need not be name-collated or indexed.
    cache_dir : str or pathlib.Path
        Directory containing immutable indexes, separated by BAM bytes,
        eligibility policy and implementation identity.
    read_collector : ReadCollector, optional
        Original-record eligibility policy. Path/assembly budgets do not change
        the library index and remain specific to each event.
    read_filter_identity : str, optional
        Explicit version/provenance for a custom record predicate, including
        closed-over data. Required for a custom ``read_filter`` or overridden
        ``alignment_filter_reason``. Callers must update it when behavior changes.
    reference_name : str, optional
        Assembly context for Varcode's contig resolver. Actual source names
        always take priority; annotation aliases are not a coordinate liftover.

    Returns
    -------
    PreparedRnaLibrary
        An index descriptor. ``open()`` gives an independent context-managed
        reader, reusable across event views. Source bytes are hashed once here;
        in-place changes are checked during use. A changed source/policy/code
        creates a new index. Incomplete preparation is never published.
    """
    source, cache_dir = Path(source).resolve(), Path(cache_dir).resolve()
    collector = ReadCollector() if read_collector is None else read_collector
    before = _file_state(source)
    identity = dict(schema=_SCHEMA, source_sha256=_digest(source), reference_name=reference_name,
                    eligibility=_eligibility_policy(collector, read_filter_identity), software=software_identity())
    if _file_state(source) != before:
        raise ValueError("Original RNA changed while computing its index identity")
    identifier = content_identifier("rna_library", identity)
    cache_dir.mkdir(parents=True, exist_ok=True)
    directory = cache_dir / identifier
    if not directory.exists():
        with TemporaryDirectory(prefix=".sv-library-building-", dir=cache_dir) as temporary:
            building = Path(temporary)
            database = building / "records.sqlite"
            with closing(_connect(database)) as db:
                statistics = _index_records(source, db, collector)
                db.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
                db.executemany("INSERT INTO metadata VALUES (?, ?)",
                    [("identity", canonical_json(identity)), ("statistics", canonical_json(statistics))])
                db.commit()
            if _file_state(source) != before:
                raise ValueError("Original RNA changed during RNA library preparation")
            receipt = dict(identity=identity, statistics=statistics, database_sha256=_digest(database))
            (building / "receipt.json").write_text(canonical_json(receipt))
            try:
                building.rename(directory)
            except OSError as error:
                if error.errno not in (errno.EEXIST, errno.ENOTEMPTY):
                    raise
                # Another preparer may have published the same immutable key.
                # Its complete receipt and database are verified below.
    receipt = json.loads((directory / "receipt.json").read_text())
    if receipt["identity"] != identity or receipt["database_sha256"] != _digest(directory / "records.sqlite"):
        raise ValueError("RNA library cache identity or database digest differs")
    with closing(_connect(directory / "records.sqlite", readonly=True)) as db:
        metadata = {key: json.loads(value) for key, value in db.execute("SELECT key, value FROM metadata")}
    if metadata != dict(identity=identity, statistics=receipt["statistics"]):
        raise ValueError("RNA library receipt differs from its pinned database metadata")
    if _file_state(source) != before:
        raise ValueError("Original RNA changed while opening its RNA library index")
    return PreparedRnaLibrary(source, directory, canonical_json(identity), canonical_json(receipt["statistics"]), before,
                             _file_state(directory / "records.sqlite"))


class RnaLibraryIndex(Mapping[Segment, list[pysam.AlignedSegment]]):
    """An independent read-only index/BAM context for complete-template queries.

    The context is reusable across views, but not concurrently across threads.
    Open separate contexts for concurrent readers. Region overlap only selects
    candidate templates; every eligible original mate/placement is returned.
    """

    def __init__(self, prepared: PreparedRnaLibrary) -> None:
        self.prepared = prepared
        self.database = prepared.directory / "records.sqlite"
        self._database_state = prepared.database_state
        self._closed = False
        self.validate()
        self.connection = _connect(self.database, readonly=True)
        try:
            self.reader = pysam.AlignmentFile(prepared.source, "rb")
        except Exception:
            self.connection.close()
            raise
        self.header = self.reader.header
        self.resolver = ContigResolver(self.reader.references, reference_name=prepared.identity["reference_name"])
        self._statistics = prepared.statistics

    @property
    def statistics(self) -> dict[str, Any]:
        """Return an independent snapshot of the pinned original-library counts."""
        return self.prepared.statistics

    def validate_eligibility(self, collector: ReadCollector, *, read_filter_identity: str | None = None) -> None:
        """Reject changed collection policy before reusing original eligibility."""
        self.validate()
        if _eligibility_policy(collector, read_filter_identity) != self.prepared.identity["eligibility"]:
            raise ValueError("RNA library index eligibility differs from the current read collector")

    def validate(self) -> None:
        """Reject closed contexts or source/cache mutation during this analysis."""
        if self._closed:
            raise ValueError("RNA library index is closed")
        if (_file_state(self.prepared.source) != self.prepared.source_state or
                _file_state(self.database) != self._database_state):
            raise ValueError("Original RNA or its RNA library index changed during use")

    def __enter__(self) -> "RnaLibraryIndex":
        self.validate()
        return self

    def __exit__(self, exc_type: object, exc_value: object, traceback: object) -> None:
        try:
            if exc_type is None:
                self.validate()
        finally:
            self.close()

    def close(self) -> None:
        """Release this reader without modifying the shared index."""
        if not self._closed:
            self.reader.close()
            self.connection.close()
            self._closed = True

    def __len__(self) -> int:
        return self._statistics["eligible_segments"]

    def __contains__(self, identity: object) -> bool:
        self.validate()
        if not isinstance(identity, tuple) or len(identity) != 3:
            return False
        group, name, segment = identity
        if segment not in (0, 64, 128, 192):
            return False
        return self.connection.execute("SELECT 1 FROM templates WHERE rg=? AND name=? AND eligible_segments & ? != 0",
                                       (group, name, 1 << (segment >> 6))).fetchone() is not None

    def __iter__(self) -> Iterator[Segment]:
        self.validate()
        for group, name, bits in self.connection.execute("SELECT rg, name, eligible_segments FROM templates ORDER BY rg, name"):
            for segment in (0, 64, 128, 192):
                if bits & (1 << (segment >> 6)):
                    yield group, name, segment

    def __getitem__(self, identity: Segment) -> list[pysam.AlignedSegment]:
        records = self.template(identity[:2]).get(identity)
        if not records:
            raise KeyError(identity)
        return records

    def template(self, fragment: Fragment, *, eligible_only: bool = True) -> dict[Segment, list[pysam.AlignedSegment]]:
        """Retrieve all distinct original mates and placements of one template."""
        self.validate()
        restriction = "AND excluded=0" if eligible_only else ""
        records: dict[Segment, list[pysam.AlignedSegment]] = defaultdict(list)
        for offset, in self.connection.execute("SELECT offset FROM records WHERE rg=? AND name=? "
                "AND duplicate_of IS NULL " + restriction + " ORDER BY segment, rid, offset", fragment):
            read = _read_at(self.reader, offset)
            records[segment_identity(read)].append(read)
        self.validate()
        return dict(records)

    @contextmanager
    def _template_selection(self, regions: Sequence[Region] | None, eligible_only: bool,
                            fragments: Iterable[Fragment] = ()) -> Iterator[sqlite3.Connection]:
        """Keep potentially large regional selections in a private temporary table."""
        with closing(_connect(self.database, readonly=True)) as db:
            if regions is None:
                db.execute("CREATE TEMP VIEW selected AS SELECT rg, name FROM main.templates "
                           + ("WHERE eligible_segments != 0" if eligible_only else ""))
                yield db
                return
            db.execute("CREATE TEMP TABLE selected (rg TEXT, name TEXT, "
                       "PRIMARY KEY (rg, name)) WITHOUT ROWID")
            db.executemany("INSERT OR IGNORE INTO selected SELECT rg, name FROM templates WHERE rg=? AND name=?",
                           fragments)
            restriction = "AND r.excluded=0" if eligible_only else ""
            for contig, start, end in regions or ():
                if type(start) is not int or type(end) is not int or not 0 <= start < end:
                    raise ValueError("RNA index queries require nonnegative half-open integer intervals")
                resolution = self.resolver.resolve(contig)
                if resolution.resolved_name is None:
                    raise ValueError(resolution.diagnostic)
                reference = self.reader.get_tid(resolution.resolved_name)
                db.execute("INSERT OR IGNORE INTO selected SELECT r.rg, r.name FROM placements p "
                    "JOIN records r USING (offset) "
                    "WHERE p.ref_min<=? AND p.ref_max>=? AND p.start_min<? AND p.end_max>? "
                    "AND r.duplicate_of IS NULL " + restriction,
                    (reference, reference, end, start))
            yield db

    def templates(self, regions: Sequence[Region] | None = None, *, eligible_only: bool = True
                  ) -> Iterator[list[pysam.AlignedSegment]]:
        """Yield complete selected templates in stable (read-group, query-name) order."""
        self.validate()
        with self._template_selection(regions, eligible_only) as selected:
            for fragment in selected.execute("SELECT rg, name FROM selected ORDER BY rg, name"):
                yield [read for records in self.template(fragment, eligible_only=eligible_only).values() for read in records]

    def records_for_templates(self, fragments: Iterable[Fragment] = (), regions: Sequence[Region] = ()
                              ) -> Iterator[pysam.AlignedSegment]:
        """Yield all physical records of selected templates in original BAM order.

        Select the union of explicit template identities and region-overlapping
        templates. An empty union is empty, not the whole library. Preserve
        excluded records and duplicate payloads for the assessment layer, which
        retains its own exact sequence/placement deduplication semantics.
        """
        if regions is None:
            raise ValueError('Selected original records require explicit regions; use original_records for the whole library')
        self.validate()
        with self._template_selection(regions, False, fragments) as selected:
            offsets = selected.execute("SELECT r.offset FROM records r JOIN selected USING (rg, name) ORDER BY r.offset")
            for offset, in offsets:
                yield _read_at(self.reader, offset)
        self.validate()

    def original_records(self) -> Iterator[tuple[int, pysam.AlignedSegment]]:
        """Scan the unchanged original library, including excluded and unmapped RNA."""
        self.validate()
        with pysam.AlignmentFile(self.prepared.source, "rb") as reader:
            while True:
                offset = reader.tell()
                read = next(reader, None)
                if read is None:
                    break
                yield offset, read
        self.validate()


@contextmanager
def _original_library_context(bam: pysam.AlignmentFile, collector: ReadCollector, directory: str | Path, *,
                              library_index: RnaLibraryIndex | None, reference_name: str | None,
                              read_filter_identity: str | None) -> Iterator[RnaLibraryIndex]:
    """Use the same original-template index for persistent and private contexts."""
    if library_index is not None:
        if not isinstance(library_index, RnaLibraryIndex):
            raise TypeError('library_index must be an open RnaLibraryIndex context')
        if not bam.is_bam or not bam.filename or Path(fsdecode(bam.filename)).resolve() != library_index.prepared.source:
            raise ValueError('library_index must refer to the same original BAM')
        if library_index.prepared.identity['reference_name'] not in (None, reference_name):
            raise ValueError('library_index reference context differs from the event')
        library_index.validate_eligibility(collector, read_filter_identity=read_filter_identity)
        yield library_index
        library_index.validate()
        return
    with ExitStack() as stack:
        if bam.is_bam and bam.filename and Path(fsdecode(bam.filename)).is_file():
            source = Path(fsdecode(bam.filename))
        else:
            source = Path(directory) / 'original.bam'
            reader = bam
            if bam.is_sam and bam.filename and Path(fsdecode(bam.filename)).is_file():
                reader = stack.enter_context(pysam.AlignmentFile(bam.filename))
            elif not bam.is_sam:
                bam.reset()
            # pysam requires byte-valued format options; its current annotation
            # incorrectly declares str (the SAM/BAM conversion test covers this).
            with pysam.AlignmentFile(source, 'wb', header=bam.header, format_options=[b'level=1']) as output:  # type: ignore[list-item]
                for read in reader.fetch(until_eof=True):
                    output.write(read)
        if read_filter_identity is None and _uses_custom_eligibility(collector):
            # This index has one private lifetime and is never reused by another
            # analysis. Public persistent indexes require a declared identity.
            read_filter_identity = 'private-context-custom-eligibility'
        prepared = prepare_rna_library(source, Path(directory) / 'library', read_collector=collector,
            reference_name=reference_name, read_filter_identity=read_filter_identity)
        yield stack.enter_context(prepared.open())
