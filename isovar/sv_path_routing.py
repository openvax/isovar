"""Route one original library to all retained SV paths before view assessment.

Seed matches select candidate templates only. Existing exact sequence, exon,
junction and alternative-placement checks still determine scientific evidence.
"""

from collections import defaultdict
from collections.abc import Iterable, Iterator, Sequence
from contextlib import closing, contextmanager
from dataclasses import dataclass
from pathlib import Path
import sqlite3

import numpy as np

from .fusion import FusionReference
from .read_identity import read_group
from .rna_library import RnaLibraryIndex
from .sv_rna_assessment import PathEvidence, PathEvidenceScope, RawRead, evidence_from_records
from .sv_rna_path_support import reverse_complement


@dataclass(frozen=True)
class PathEvidenceQuery:
    """One view's retained sequences, exact reference models and thresholds."""

    identifier: str
    sequences: tuple[str, ...]
    references: tuple[FusionReference, ...]
    min_path_match: int = 20
    min_exon_overlap: int = 15

    def __post_init__(self) -> None:
        if not isinstance(self.identifier, str) or not self.identifier:
            raise ValueError('Path evidence queries require a nonempty identity')
        if any(type(value) is not int or value < 1 for value in (self.min_path_match, self.min_exon_overlap)):
            raise ValueError('Path evidence thresholds must be positive integers')
        if not isinstance(self.sequences, tuple) or not all(isinstance(s, str) for s in self.sequences):
            raise TypeError('Path evidence sequences must be a tuple of strings')
        if not isinstance(self.references, tuple) or not all(isinstance(r, FusionReference) for r in self.references):
            raise TypeError('Path evidence references must be a tuple of FusionReference models')


class _SeedRouter:
    """Use exact fixed-width byte searches; never interpret ambiguous DNA bases."""

    def __init__(self, queries: Sequence[PathEvidenceQuery]) -> None:
        self.size = min((min(15, query.min_path_match) for query in queries), default=15)
        owners: dict[str, set[str]] = defaultdict(set)
        for query in queries:
            for sequence in query.sequences:
                if len(sequence) < query.min_path_match:
                    continue
                # Any exact min_path_match interval contains this many possible
                # seed starts. One seed per such interval therefore cannot drop
                # a qualifying exact match, and bounds inventory-wide seed size.
                stride = query.min_path_match - self.size + 1
                for start in range(0, len(sequence) - self.size + 1, stride):
                    word = sequence[start:start + self.size]
                    owners[word].add(query.identifier)
                    owners[reverse_complement(word)].add(query.identifier)
        self.owners = dict(owners)
        self.words = sorted(owners)
        self.ascii = all(word.isascii() and '\0' not in word for word in self.words)
        self.encoded = np.array(self.words, dtype=f'S{self.size}') if self.ascii else None

    def matches(self, sequence: str) -> set[str]:
        """Return conservative view IDs, preserving literal legacy matching."""
        if len(sequence) < self.size or not self.words:
            return set()
        sequence = sequence.upper()
        if not self.ascii or not sequence.isascii() or '\0' in sequence:
            return {identifier for start in range(len(sequence) - self.size + 1)
                    for identifier in self.owners.get(sequence[start:start + self.size], ())}
        assert self.encoded is not None  # ASCII route words were encoded at construction.
        encoded = np.frombuffer(sequence.encode('ascii'), dtype='S1')
        windows = np.lib.stride_tricks.sliding_window_view(encoded, self.size).view(f'S{self.size}').reshape(-1)
        positions = np.searchsorted(self.encoded, windows)
        within = positions < len(self.words)
        positions, windows = positions[within], windows[within]
        exact = np.unique(positions[self.encoded[positions] == windows])
        return {identifier for position in exact for identifier in self.owners[self.words[position]]}


@dataclass(frozen=True)
class PathEvidenceIndex:
    """Read a completed source-specific route database within its owning context."""

    library: RnaLibraryIndex
    connection: sqlite3.Connection
    queries: tuple[PathEvidenceQuery, ...]

    def evidence(self, identifier: str) -> PathEvidence:
        """Collect original-order records and all view-specific exon alignments."""
        query = next((query for query in self.queries if query.identifier == identifier), None)
        if query is None:
            raise KeyError(identifier)
        self.library.validate()
        fragments = self.connection.execute('SELECT rg, name FROM matches WHERE query=?', (identifier,))
        # The exhaustive collector compares model contigs with exact SAM names.
        # Models on absent contigs therefore have no exon alignments; do not
        # introduce alias-derived placements as part of this performance change.
        names = set(self.library.header.references)
        regions = sorted({(reference.contig, start, end) for reference in query.references
                          if reference.contig in names for start, end in reference.exons})
        raw = self.connection.execute('SELECT rg, name, mate, sequence FROM raw_matches WHERE query=?', (identifier,))
        counts = self.library.statistics
        evidence = evidence_from_records(self.library.records_for_templates(fragments, regions), query.references,
            raw, query.min_exon_overlap, source_records=counts['records_scanned'],
            source_fragments=counts['fragments_scanned'])
        evidence['scope'] = PathEvidenceScope(query.references, query.min_exon_overlap,
                                               query.sequences, query.min_path_match)
        return evidence


@contextmanager
def route_path_evidence(library: RnaLibraryIndex, queries: Sequence[PathEvidenceQuery],
                         database: str | Path, *, raw_reads: Iterable[RawRead] = ()) -> Iterator[PathEvidenceIndex]:
    """Scan original BAM/FASTQ sequences once and retain matched identities on disk.

    Parameters
    ----------
    library : RnaLibraryIndex
        Open context for one unchanged original library. Its eligibility filter
        does not restrict this assessment scan or retrieved template records.
    queries : sequence of PathEvidenceQuery
        Unique view IDs, retained paths, reference models and match thresholds.
        A short minimum across views can overselect; exact assessment is unchanged.
    database : str or pathlib.Path
        New private SQLite path. Closed on exit; its parent owns file cleanup.
    raw_reads : iterable of tuple, optional
        Original ``(read_group, query_name, mate, sequence)`` ends, consumed once.
        These add sequence witnesses, never genomic placements or library counts.

    Yields
    ------
    PathEvidenceIndex
        Context-bound per-view evidence. No-path views still retain all reference
        exon counts and original-library denominators without a sequence scan.
    """
    queries = tuple(queries)
    if len({query.identifier for query in queries}) != len(queries):
        raise ValueError('Duplicate path evidence query identity')
    reference_name = library.prepared.identity['reference_name']
    if reference_name is not None and any(reference.reference_name != reference_name
            for query in queries for reference in query.references):
        raise ValueError('Path evidence reference assembly differs from its original-library context')
    database = Path(database)
    if database.exists():
        raise ValueError('Path evidence routing requires a new private database')
    library.validate()
    router = _SeedRouter(queries)
    with closing(sqlite3.connect(database)) as connection:
        connection.execute('PRAGMA temp_store=FILE')
        connection.execute('PRAGMA cache_size=-32768')
        connection.execute('CREATE TABLE matches (query TEXT, rg TEXT, name TEXT, '
                           'PRIMARY KEY (query, rg, name)) WITHOUT ROWID')
        connection.execute('CREATE TABLE raw_matches (query TEXT, rg TEXT, name TEXT, mate INTEGER, sequence TEXT, '
                           'PRIMARY KEY (query, rg, name, mate, sequence)) WITHOUT ROWID')
        if router.words:
            for count, (_, read) in enumerate(library.original_records(), 1):
                if read.query_sequence:
                    identifiers = router.matches(read.query_sequence)
                    if identifiers:
                        group, name = read_group(read), read.query_name or ''
                        connection.executemany('INSERT OR IGNORE INTO matches VALUES (?, ?, ?)',
                                               ((identifier, group, name) for identifier in identifiers))
                if count % 100000 == 0:
                    connection.commit()
        for group, name, mate, sequence in raw_reads:
            if mate not in (0, 1, 2) or not sequence:
                raise ValueError('Raw RNA needs a sequence and a mate number in 0, 1, 2')
            identifiers = router.matches(sequence)
            connection.executemany('INSERT OR IGNORE INTO raw_matches VALUES (?, ?, ?, ?, ?)',
                ((identifier, group, name, mate, sequence.upper()) for identifier in identifiers))
            connection.executemany('INSERT OR IGNORE INTO matches VALUES (?, ?, ?)',
                ((identifier, group, name) for identifier in identifiers))
        connection.commit()
        library.validate()
        yield PathEvidenceIndex(library, connection, queries)
        library.validate()
