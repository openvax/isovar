"""Disk-backed, bounded discovery followed by complete-input SV support.

Only candidate witnesses are retained in memory during support reporting. The
SQLite index keeps unrelated background and arbitrarily distant pieces on disk.
"""

from collections import Counter, defaultdict
from bisect import bisect_left
from collections.abc import Mapping
from contextlib import ExitStack, closing
from hashlib import sha256
from os import fsdecode
from pathlib import Path
import sqlite3
from tempfile import TemporaryDirectory

import pysam

from .cell_umi import CellUmiEvidence
from .read_identity import segment_identity
from .read_lineage import ReadLineage
from .read_metadata import record_evidence
from .rna_evidence import rna_support
from .sv_rna import (
    _Event, _GeometryCache, _Model, _Observations, _build_segment, _flip, _ineligible, _join_class, _merge, _record_gaps, _record_id,
    _path_result, _possible_breakpoint_clip, reconstruct_sv_rna,
)
from . import sv_rna_stream
from .sv_rna_relations import LINKAGE_STATUSES


def _read_at(reader, offset):
    """Read one original record, reusing an already decompressed BGZF block."""
    current = reader.tell()
    if current != offset:
        if current < offset and current >> 16 == offset >> 16:
            while reader.tell() < offset:
                next(reader)
        else:
            reader.seek(offset)
    return next(reader)


class _DiskGroups(Mapping):
    """Original eligible BAM records, indexed without copying their payloads."""

    def __init__(self, connection, reader, header, segment_count):
        self.connection, self.reader, self.header = connection, reader, header
        self.segment_count = segment_count

    def __getitem__(self, identity):
        records = self.template(identity[:2]).get(identity)
        if not records:
            raise KeyError(identity)
        return records

    def __contains__(self, identity):
        rg, name, segment = identity
        return self.connection.execute(
            'SELECT 1 FROM templates WHERE rg=? AND name=? AND segments & ? != 0',
            (rg, name, 1 << (segment >> 6))).fetchone() is not None

    def __iter__(self):
        for rg, name, segments in self.connection.execute(
                'SELECT rg, name, segments FROM templates ORDER BY rg, name'):
            for segment in (0, 64, 128, 192):
                if segments & (1 << (segment >> 6)):
                    yield rg, name, segment

    def __len__(self):
        return self.segment_count

    def templates(self):
        # Traverse first appearances in physical input order. Singletons need
        # neither a second SQL lookup nor a random seek for every record.
        with pysam.AlignmentFile(self.reader.filename, 'rb') as reader:
            for rg, name, offset, count in self.connection.execute(
                    'SELECT rg, name, offset, records FROM templates ORDER BY offset'):
                if count == 1:
                    yield [_read_at(reader, offset)]
                else:
                    yield [r for records in self.template((rg, name)).values() for r in records]

    def template(self, fragment):
        """All eligible original mates/placements of one consulted template."""
        groups = defaultdict(list)
        for offset, in self.connection.execute(
                'SELECT o.offset FROM offsets o LEFT JOIN record_ids r USING (offset) '
                'WHERE o.rg=? AND o.name=? AND o.excluded=0 ORDER BY o.segment, r.rid', fragment):
            read = _read_at(self.reader, offset)
            groups[segment_identity(read)].append(read)
        return groups


def _spool(bam, db, collector, directory):
    """Index original BAM offsets; serialize only potential duplicate payloads.

    QNAME/RG identify the only records that could share an original SAM payload.
    Singletons therefore need no SAM serialization, hash, or copied BAM. All
    records in repeated templates are hashed in physical file order; identical
    originals are deduplicated while exclusions still count every input record.
    """
    directory = Path(directory)
    if bam.is_bam and bam.filename and Path(fsdecode(bam.filename)).is_file():
        source = Path(fsdecode(bam.filename))
    else:
        source = directory / 'original.bam'
        with ExitStack() as stack:
            reader = bam
            if bam.is_sam and bam.filename and Path(fsdecode(bam.filename)).is_file():
                reader = stack.enter_context(pysam.AlignmentFile(bam.filename))
            elif not bam.is_sam:
                bam.reset()
            output = stack.enter_context(pysam.AlignmentFile(
                source, 'wb', header=bam.header, format_options=[b'level=1']))
            for read in reader.fetch(until_eof=True):
                output.write(read)
    db.execute('CREATE TABLE offsets (offset INTEGER PRIMARY KEY, rg TEXT, name TEXT, segment INTEGER, excluded INTEGER)')
    scanned, excluded = 0, Counter()
    with pysam.AlignmentFile(source, 'rb') as reader:
        while True:
            offset = reader.tell()
            read = next(reader, None)
            if read is None:
                break
            scanned += 1
            reason = 'unnamed_record' if read.query_name is None else _ineligible(read, collector)
            if reason:
                excluded[reason] += 1
            rg, name, segment = segment_identity(read)
            db.execute('INSERT INTO offsets VALUES (?, ?, ?, ?, ?)',
                       (offset, rg, name or '', segment, int(reason is not None)))
    db.execute('CREATE INDEX fragments ON offsets (rg, name)')
    db.execute('CREATE TABLE repeated (rg TEXT, name TEXT, PRIMARY KEY (rg, name)) WITHOUT ROWID')
    db.execute('INSERT INTO repeated SELECT rg, name FROM offsets GROUP BY rg, name HAVING count(*) > 1')
    db.execute('CREATE TABLE record_ids (rid BLOB PRIMARY KEY, offset INTEGER UNIQUE)')
    db.execute('CREATE TABLE duplicates (offset INTEGER PRIMARY KEY)')
    duplicates = 0
    with pysam.AlignmentFile(source, 'rb') as reader, pysam.AlignmentFile(source, 'rb') as previous:
        for offset, in db.execute('SELECT o.offset FROM offsets o JOIN repeated USING (rg, name) ORDER BY o.offset'):
            sam = _read_at(reader, offset).to_string()
            rid = sha256(sam.encode()).digest()
            row = db.execute('INSERT OR IGNORE INTO record_ids VALUES (?, ?)', (rid, offset))
            if not row.rowcount:
                old_offset, = db.execute('SELECT offset FROM record_ids WHERE rid=?', (rid,)).fetchone()
                if _read_at(previous, old_offset).to_string() != sam:
                    raise ValueError('Conflicting original SAM records share SV record identity %s' % rid.hex())
                db.execute('INSERT INTO duplicates VALUES (?)', (offset,))
                duplicates += 1
    db.execute('DELETE FROM offsets WHERE offset IN (SELECT offset FROM duplicates)')
    db.execute('DROP TABLE duplicates')
    db.execute('DROP TABLE repeated')
    db.execute('CREATE TABLE templates (rg TEXT, name TEXT, offset INTEGER, records INTEGER, '
               'segments INTEGER, PRIMARY KEY (rg, name)) WITHOUT ROWID')
    db.execute('INSERT INTO templates SELECT rg, name, min(offset), count(*), '
               'sum(DISTINCT (1 << (segment >> 6))) FROM offsets WHERE excluded=0 GROUP BY rg, name')
    db.execute('CREATE INDEX traversal ON templates (offset)')
    segment_count, = db.execute('SELECT count(*) FROM (SELECT rg, name, segment FROM offsets '
                                'WHERE excluded=0 GROUP BY rg, name, segment)').fetchone()
    db.commit()
    return scanned, dict(sorted(excluded.items())), scanned - duplicates, segment_count, source


def _search_regions(bam, options):
    lengths = dict(zip(bam.references, bam.lengths))
    regions = list(options['regions']) + [(r.contig, a, b) for r in options['references'] for a, b in r.exons]
    for breakpoint in (options['donor'], options['acceptor']):
        if breakpoint.contig not in lengths:
            raise ValueError('Breakpoint contig %r is absent from the alignments' % breakpoint.contig)
        regions.append((breakpoint.contig, max(0, breakpoint.position - options['breakpoint_window']),
                        min(lengths[breakpoint.contig], breakpoint.position + options['breakpoint_window'])))
    spans = defaultdict(list)
    for contig, start, end in regions:
        if (contig not in lengths or type(start) is not int or type(end) is not int
                or not 0 <= start < end <= lengths[contig]):
            raise ValueError('Invalid or unavailable SV search interval: %r' % ((contig, start, end),))
        spans[contig].append((start, end))
    return spans


def _priority(records, collector, event, annotated):
    """Placed event-related joins, breakpoint clips, other joins, then context."""
    clip = False
    possible_clip = collector.use_soft_clipped_bases and any(
        _possible_breakpoint_clip(read, [(b.contig, b.position) for b in (event.donor, event.acceptor)])
        for read in records)
    if len(records) > 1 or possible_clip:
        # Actual supplementary/alternative records stay together. An SA hint
        # alone neither raises priority nor fabricates its missing record.
        store = _Observations(records, collector, 1, {}, '', '')
        for identity in store.groups:
            for o in store.build(identity):
                clip |= bool(event.clips(o.positions, len(o.sequence)))
                if any(_join_class(event, annotated, o.positions[i], o.positions[j], j - i - 1, kind)[0]
                       in ('event', 'lazy') for i, j, kind in o.breaks):
                    return 0
    priority = 3
    for read in records:
        for left, right, kind, inserted in _record_gaps(read):
            classes = [_join_class(event, annotated, left, right, len(inserted), kind),
                       _join_class(event, annotated, (right[0], right[1], '-'),
                                   (left[0], left[1], '-'), len(inserted), kind)]
            if any(use in ('event', 'lazy') for use, _ in classes):
                return 0
            if not any(use in ('annotated', 'wobble') for use, _ in classes) and kind != 'D':
                priority = 2
    return 1 if clip else priority


def _group_signatures(identity, records):
    """Full placement bins subdivided by every observed base, without QUAL/labels.

    Actual alternative/supplementary pieces, SA linkage hints and trimming tags
    remain in the key. No sketch, local window or consensus replaces an allele.
    Original qualities, names and labels remain in the spool, not in this key.
    """
    placed = sorted((r.reference_name, r.reference_start, r.flag, r.cigarstring,
                     r.get_tag('SA') if r.has_tag('SA') else '',
                     repr(r.get_tag('pt')) if r.has_tag('pt') else '', r.query_sequence)
                    for r in records)
    structure = (identity[0], identity[2], [r[:-1] for r in placed])
    return (sha256(repr(structure).encode()).digest(),
            sha256(repr((identity[0], identity[2], placed)).encode()).digest())


def _select(groups, db, options, collector, spans, event, annotated):
    priorities = {}
    db.execute('CREATE TABLE alleles (group_id INTEGER PRIMARY KEY, priority INTEGER, signature BLOB, '
               'structure BLOB, UNIQUE (priority, signature))')
    db.execute('CREATE TABLE discovery (rg TEXT, name TEXT, segment INTEGER, group_id INTEGER, size INTEGER)')
    for template in groups.templates():
        members = defaultdict(list)
        for read in template:
            members[segment_identity(read)].append(read)
        for identity, records in members.items():
            if not any(a < end and start < b for read in records for a, b in read.get_blocks()
                       for start, end in spans.get(read.reference_name, ())):
                continue
            structure, signature = _group_signatures(identity, records)
            # QUAL and labels do not affect junction/clip priority. Six pieces
            # have at most 64 subsets, below the segment-grouping truncation
            # boundary; larger bundles must be validated separately because
            # ordering under that boundary can depend on original record IDs.
            cacheable = len(records) <= 6
            cached = priorities.get(signature) if cacheable else None
            if cached is None:
                priority = _priority(records, collector, event, annotated)
                cursor = db.execute('INSERT OR IGNORE INTO alleles (priority, signature, structure) VALUES (?, ?, ?)',
                                    (priority, signature, structure))
                group_id = (cursor.lastrowid if cursor.rowcount else db.execute(
                    'SELECT group_id FROM alleles WHERE priority=? AND signature=?', (priority, signature)).fetchone()[0])
                if cacheable:
                    if len(priorities) >= 4096:
                        priorities.pop(next(iter(priorities)))
                    priorities[signature] = group_id
            else:
                group_id = cached
            db.execute('INSERT INTO discovery VALUES (?, ?, ?, ?, ?)', (*identity, group_id, len(records)))
    db.commit()
    total = db.execute('SELECT coalesce(sum(size), 0) FROM discovery').fetchone()[0]
    selected, details = _select_group_representatives(groups, db, options)
    return selected, total, details


def _select_group_representatives(groups, db, options):
    """Visit distinct alleles first, retaining abundance inside the record cap.

    Representative-only sampling would flatten major/minor allele ratios. Keep
    the original voters; binning shares calculations, not biological counts.
    """
    query = '''SELECT group_id, rg, name, segment, size FROM (
        SELECT *, row_number() OVER (PARTITION BY group_id ORDER BY rg, name, segment) AS copy
        FROM discovery JOIN alleles USING (group_id)) ORDER BY priority, (copy - 1) / ?, signature, copy'''
    rows = db.execute(query, (options['min_alternative_fragments'],))
    selected, selected_groups, count, record_limited = [], set(), 0, False
    for group_id, rg, name, segment, size in rows:
        if count + size > options['max_records']:
            record_limited = True
            continue
        selected.extend(groups[rg, name, segment])
        count += size
        selected_groups.add(group_id)
    db.execute('CREATE TABLE selected_groups (group_id INTEGER PRIMARY KEY)')
    db.executemany('INSERT INTO selected_groups VALUES (?)', ((group_id,) for group_id in selected_groups))
    total_groups = db.execute('SELECT count(*) FROM alleles').fetchone()[0]
    represented = db.execute('SELECT coalesce(sum(size), 0) FROM discovery '
                             'JOIN selected_groups USING (group_id)').fetchone()[0]
    return selected, dict(
        eligible_structural_bins=db.execute('SELECT count(DISTINCT structure) FROM alleles').fetchone()[0],
        eligible_regional_groups=total_groups, selected_groups=len(selected_groups),
        omitted_groups=total_groups - len(selected_groups), represented_records=represented,
        groups_complete=total_groups == len(selected_groups), record_limit_reached=record_limited,
        initial_representatives_per_group=options['min_alternative_fragments'])


def _positions(path):
    positions = [None] * len(path['sequence'])
    for block in path['blocks']:
        bases = range(block['reference_start'], block['reference_end'])
        if block['strand'] == '-':
            bases = reversed(bases)
        for q, base in zip(range(block['query_start'], block['query_end']), bases):
            positions[q] = (block['contig'], base, block['strand'])
    return tuple(positions)


def _may_support(template, joins, collector, event):
    """Conservative CIGAR screen: false positives are built and checked normally.

    Mates never form a split read. More than one actual record of a segment
    needs full validation, including alternative placements and missing pieces.
    """
    identities = set()
    for read in template:
        identity = segment_identity(read)
        if identity in identities:
            return True
        identities.add(identity)
        if collector.use_soft_clipped_bases and _possible_breakpoint_clip(
                read, [(b.contig, b.position) for b in (event.donor, event.acceptor)]):
            return True
        for left, right, _, inserted in _record_gaps(read):
            for a, b in ((left, right), (_flip(right), _flip(left))):
                if (a, b) in joins or event.relation(a, b, len(inserted))[0] == 'breakpoint_junction':
                    return True
    return False


def _overlaps_placements(read, intervals):
    """Whether a CIGAR block overlaps the disjoint union of candidate blocks."""
    return any((i := bisect_left(intervals, (b,))) and intervals[i - 1][1] > a
               for a, b in read.get_blocks()) if intervals else False


def _support(groups, result, collector, event, annotated, options):
    """Stream complete templates, keeping only candidate junction witnesses."""
    header = groups.header.to_dict()
    joins = {(tuple(j['left']), tuple(j['right'])) for path in result['paths'] for j in path['junctions']
             if j['left'] is not None and j['right'] is not None}
    placements = defaultdict(list)
    for path in result['paths']:
        for block in path['blocks']:
            placements[block['contig']].append((block['reference_start'], block['reference_end']))
    placements = {contig: _merge(intervals) for contig, intervals in placements.items()}
    stream = sv_rna_stream._SupportStream(result['paths'], [_positions(p) for p in result['paths']], event, annotated)
    reasons, geometry_cache, batch, bases, intern = Counter(), _GeometryCache(), [], 0, {}
    for template in groups.templates() if result['paths'] else ():
        if not any(_overlaps_placements(read, placements.get(read.reference_name, ())) for read in template):
            continue
        if not _may_support(template, joins, collector, event):
            continue
        members = defaultdict(list)
        # Share overlapping coordinates without accumulating an input-wide
        # interning table. An indivisible template can exceed this bound.
        if len(intern) >= sv_rna_stream._SUPPORT_BATCH_BASES:
            intern.clear()
        for read in template:
            members[segment_identity(read)].append(read)
        observations = [o for identity, records in members.items()
                        for o in _build_segment(identity, records, collector, reasons, intern, geometry_cache)]
        hit = any((o.positions[i], o.positions[j]) in joins
                  or event.relation(o.positions[i], o.positions[j], j - i - 1)[0] == 'breakpoint_junction'
                  for o in observations for i, j, _ in o.breaks)
        hit |= any(event.clips(o.positions, len(o.sequence)) for o in observations)
        if hit:
            size = sum(len(o.sequence) for o in observations)
            if batch and bases + size > sv_rna_stream._SUPPORT_BATCH_BASES:
                stream.consume(batch)
                batch, bases = [], 0
            batch.extend(observations)
            bases += size
            if bases >= sv_rna_stream._SUPPORT_BATCH_BASES:
                stream.consume(batch)
                batch, bases = [], 0
    stream.consume(batch)
    batch.clear()
    observations, intern = [], {}
    del geometry_cache
    # Only output metadata persists. Parents and mates resolve from the original
    # disk-backed input, including records outside every candidate's intervals.
    store = _Observations([], collector, options['max_extension_segments'], header,
                          options['sample_id'], options['source'])
    store.lineage = ReadLineage(groups, header)
    store.cell_umi = CellUmiEvidence(groups, header, options['sample_id'], options['source'],
                                   template_loader=groups.template)
    junction_cache, linked = stream.summaries(store.cell_umi, store.lineage)
    used = {}

    def support(witnesses):
        identities = {tuple(stream.rows[w['observation']]['identity']) for w in witnesses}
        return dict(**store.cell_umi.support(identities), read_lineage=store.lineage.support(identities),
                    missing_quality_reads=len({tuple(stream.rows[w['observation']]['identity']) for w in witnesses
                                               if stream.rows[w['observation']]['missing_qualities']}),
                    scope='complete_input_direct_junction_observations', witnesses=witnesses)

    for index, path in enumerate(result['paths']):
        positions = stream.positions[index]
        reported, evidence = _path_result(path['sequence'], positions, set(), store, event, annotated, [], {},
                                          lambda *args: {}, lambda *args: None, junction_cache=junction_cache,
                                          linked_intervals=linked[index])
        path['junctions'] = reported['junctions']
        used.update(evidence)
        path['sequence_evidence']['voting_support_scope'] = 'bounded_discovery'
        path['full_interval_support'] = support(stream.interval_witnesses(index, 0, len(positions), path['junctions']))
        for candidate in path['exploratory_orfs']['candidates']:
            candidate['start_evidence']['observation_scope'] = 'bounded_discovery'
            start, end = candidate['query_interval']
            crossed = [path['junctions'][i] for i in candidate['crossed_junctions']]
            candidate['full_interval_support'] = support(stream.interval_witnesses(index, start, end, crossed))
        stream.witnesses[index].clear()
    # Discovery-only context also needs original visible mates/parents. Its
    # voting counts remain explicitly scoped to the bounded discovery set.
    consulted = {tuple(row['identity']) for row in result['cell_umi_evidence']['reads']}
    original = CellUmiEvidence(groups, header, options['sample_id'], options['source'], template_loader=groups.template)
    for identity in sorted(consulted):
        original.segment(identity)
    labels = {tuple(row['identity']): row for row in original.evidence()['reads']}
    labels.update((tuple(row['identity']), row) for row in store.cell_umi.evidence()['reads'])
    result['cell_umi_evidence']['reads'] = [labels[key] for key in sorted(labels)]
    for path in result['paths']:
        voting = {tuple(result['observations'][key]['identity'])
                  for key in path['sequence_evidence']['voting_observations']}
        path['sequence_evidence']['voting_support'] = rna_support(voting, labels.__getitem__)
        path['sequence_evidence']['label_resolution_scope'] = 'complete_eligible_input'
    for row in result['read_lineage']['reads']:
        store.lineage.segment(tuple(row['identity']))
    lineage = {tuple(row['identity']): row for row in store.lineage.evidence()}
    result['read_lineage'] = dict(scope='complete_eligible_input', reads=[lineage[key] for key in sorted(lineage)])
    result['observations'].update(sorted(used.items()))
    cited = {rid for row in used.values() for rid in row['records']}
    for identity in {tuple(row['identity']) for row in used.values()}:
        for read in groups[identity]:
            rid = _record_id(read)
            if rid in cited and rid not in result['original_records']:
                result['original_records'][rid] = read.to_string()
                result['record_evidence'][rid] = record_evidence(read)
    return dict(reasons)


def reconstruct_dense_sv_rna(bam, options):
    """Discover bounded hypotheses and recount their witnesses across the input."""
    options = dict(options)
    scratch = options.pop('scratch_dir')
    options['dense_support'] = False
    collector = options['read_collector']
    spans = _search_regions(bam, options)
    models = [_Model(r) for r in options['references']]
    annotated = set().union(*(m.junctions() for m in models))
    event = _Event(options['donor'], options['acceptor'], options['max_breakpoint_shift'], annotated,
                   options['annotated_junction_tolerance'])

    with TemporaryDirectory(prefix='isovar-sv-', dir=scratch) as directory:
        with closing(sqlite3.connect(str(Path(directory) / 'records.sqlite'))) as db, ExitStack() as stack:
            # This private spool is rebuilt from the unchanged input and is
            # discarded on every exit. It never needs crash recovery or a
            # rollback journal containing a second copy of the original SAM.
            db.execute('PRAGMA journal_mode=OFF')
            db.execute('PRAGMA synchronous=OFF')
            db.execute('PRAGMA temp_store=FILE')
            db.execute('PRAGMA cache_size=-65536')
            scanned, excluded, unique, segment_count, eligible = _spool(bam, db, collector, directory)
            reader = stack.enter_context(pysam.AlignmentFile(eligible, 'rb'))
            groups = _DiskGroups(db, reader, bam.header, segment_count)
            selected, available, grouping = _select(groups, db, options, collector, spans, event, annotated)
            unsorted, discovery = str(Path(directory) / 'discovery.unsorted.bam'), str(Path(directory) / 'discovery.bam')
            with pysam.AlignmentFile(unsorted, 'wb', header=bam.header) as output:
                for read in selected:
                    output.write(read)
            selected_count = len(selected)
            del selected
            pysam.sort('--no-PG', '-o', discovery, unsorted)
            pysam.index(discovery)
            with pysam.AlignmentFile(discovery) as discovery_bam:
                result = reconstruct_sv_rna(discovery_bam, **options)
            support_notes = _support(groups, result, collector, event, annotated, options)
            result['parameters']['dense_support'] = True
            result['acquisition']['input_scope'] = 'bounded_discovery_representatives'
            result['support_acquisition'] = dict(scope='all_original_input_records', records_complete=True,
                                                  complete=not support_notes.get('segment_grouping_limit'),
                                                  records_scanned=scanned, unique_records=unique,
                                                  eligible_segments=len(groups),
                                                  excluded_records=excluded, read_path_notes=support_notes)
            result['discovery'] = dict(policy='event_first_full_allele_groups.v1', **grouping,
                                       eligible_regional_records=available, selected_records=selected_count,
                                       omitted_records=available - selected_count,
                                       complete=available == selected_count and not result['acquisition']['limitations'])
            if grouping['record_limit_reached']:
                result['limitations'] = sorted(set(result['limitations']) | {'discovery_record_limit'})
            if support_notes.get('segment_grouping_limit'):
                result['limitations'] = sorted(set(result['limitations']) | {'support_segment_grouping_limit'})
            result['inclusion_competing_splice_scope'] = 'bounded_discovery'
            result['seeds_scope'] = result['competing_annotated_junctions_scope'] = 'bounded_discovery'
            result['observation_counts_scope'] = 'bounded_discovery'
            rank = {status: i for i, status in enumerate(LINKAGE_STATUSES)}
            result['paths'].sort(key=lambda path: (
                rank.get(path['event_linkage']['status'], len(rank)),
                -max((j['direct_support']['fragments'] for j in path['junctions']), default=0),
                path['sequence'], path['path_id']))
            return result
