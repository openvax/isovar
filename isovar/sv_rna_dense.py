"""Disk-backed, bounded discovery followed by complete-input SV support.

Only candidate witnesses are retained in memory during support reporting. The
SQLite spool keeps unrelated background and arbitrarily distant pieces on disk.
"""

from collections import Counter, defaultdict
from collections.abc import Mapping
from contextlib import ExitStack, closing
from hashlib import sha256
from itertools import groupby
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
    _Event, _GeometryCache, _Model, _Observations, _build_segment, _flip, _ineligible, _join_class, _record_gaps, _record_id,
    _path_result, _possible_breakpoint_clip, reconstruct_sv_rna,
)
from .sv_rna_orfs import _witnesses
from .sv_rna_relations import LINKAGE_STATUSES


class _DiskGroups(Mapping):
    """Original eligible BAM records with one seek per consulted template."""

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
        with pysam.AlignmentFile(self.reader.filename, 'rb') as reader:
            for _, records in groupby(reader.fetch(until_eof=True), key=lambda r: segment_identity(r)[:2]):
                yield list(records)

    def template(self, fragment):
        """All eligible original mates/placements of one consulted template."""
        row = self.connection.execute('SELECT offset, records FROM templates WHERE rg=? AND name=?',
                                      fragment).fetchone()
        groups = defaultdict(list)
        if row is not None:
            offset, count = row
            self.reader.seek(offset)
            for _ in range(count):
                read = next(self.reader)
                groups[segment_identity(read)].append(read)
        return groups


def _spool(bam, db, collector, directory):
    """Collate binary originals; index template offsets, not expanded SAM text.

    Equal SAM payloads are deduplicated only inside their original query-name
    group. Sorting cannot separate identical records, so no genome-sized set of
    record hashes or text copy is needed. Exclusion counts still count every
    input record; unique counts include excluded originals.
    """
    directory = Path(directory)
    name_sorted, eligible = directory / 'name-sorted.bam', directory / 'eligible.bam'
    if (bam.is_bam or bam.is_sam) and bam.filename and Path(fsdecode(bam.filename)).is_file():
        source = bam.filename
    else:
        source = str(directory / 'original.bam')
        if not bam.is_sam:
            bam.reset()
        with pysam.AlignmentFile(source, 'wb', header=bam.header, format_options=[b'level=1']) as original:
            for read in bam.fetch(until_eof=True):
                original.write(read)
    # Normal collation retains supplementary, secondary and unpaired records;
    # its primary-only fast mode would discard evidence and is never used.
    pysam.collate('--no-PG', '-o', str(name_sorted), source, str(directory / 'collate'))
    db.execute('CREATE TABLE templates (rg TEXT, name TEXT, offset INTEGER, records INTEGER, '
               'segments INTEGER, PRIMARY KEY (rg, name)) WITHOUT ROWID')
    scanned, unique, segment_count, excluded = 0, 0, 0, Counter()
    with pysam.AlignmentFile(name_sorted, 'rb') as reader, pysam.AlignmentFile(
            eligible, 'wb', header=bam.header, format_options=[b'level=1']) as output:
        for _, name_records in groupby(reader.fetch(until_eof=True), key=lambda r: r.query_name):
            name_records = list(name_records)
            originals, fragments = {}, defaultdict(list)
            for read in name_records:
                scanned += 1
                reason = 'unnamed_record' if read.query_name is None else _ineligible(read, collector)
                if reason:
                    excluded[reason] += 1
                rid = None
                if len(name_records) > 1:
                    sam = read.to_string()
                    rid = sha256(sam.encode()).hexdigest()
                    if rid in originals:
                        if originals[rid] != sam:
                            raise ValueError('Conflicting original SAM records share SV record identity %s' % rid)
                        continue
                    originals[rid] = sam
                unique += 1
                if reason is None:
                    fragments[segment_identity(read)[:2]].append((rid, read))
            for fragment, records in sorted(fragments.items()):
                offset, segments = output.tell(), set()
                for _, read in sorted(records, key=lambda row: row[0]):
                    output.write(read)
                    segments.add(segment_identity(read)[2])
                segment_count += len(segments)
                db.execute('INSERT INTO templates VALUES (?, ?, ?, ?, ?)',
                           (*fragment, offset, len(records), sum(1 << (segment >> 6) for segment in segments)))
    db.commit()
    name_sorted.unlink()
    return scanned, dict(sorted(excluded.items())), unique, segment_count, eligible


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


def _support(groups, result, collector, event, annotated, options):
    """Stream complete templates, keeping only candidate junction witnesses."""
    header = groups.header.to_dict()
    joins = {(tuple(j['left']), tuple(j['right'])) for path in result['paths'] for j in path['junctions']
             if j['left'] is not None and j['right'] is not None}
    placements = {(b['contig'], b['reference_start'], b['reference_end'])
                  for path in result['paths'] for b in path['blocks']}
    records, reasons, geometry_cache, prebuilt, intern = [], Counter(), _GeometryCache(), {}, {}
    for template in groups.templates() if result['paths'] else ():
        # No per-base path building for unrelated genomic background.
        if not any(read.reference_name == contig and a < high and low < b
                   for read in template for a, b in read.get_blocks() for contig, low, high in placements):
            continue
        if not _may_support(template, joins, collector, event):
            continue
        members = defaultdict(list)
        for read in template:
            members[segment_identity(read)].append(read)
        built = {identity: _build_segment(identity, group, collector, reasons, intern, geometry_cache)
                 for identity, group in members.items()}
        observations = [o for group in built.values() for o in group]
        hit = any((o.positions[i], o.positions[j]) in joins
                  or event.relation(o.positions[i], o.positions[j], j - i - 1)[0] == 'breakpoint_junction'
                  for o in observations for i, j, _ in o.breaks)
        hit |= any(event.clips(o.positions, len(o.sequence)) for o in observations)
        if hit:
            records.extend(template)  # Include visible mates when resolving labels.
            prebuilt.update(built)
    store = _Observations(records, collector, options['max_extension_segments'], header,
                          options['sample_id'], options['source'], prebuilt=prebuilt)
    # Parents can be outside the candidates' intervals; keep original producer
    # and ancestry resolution over the complete eligible input, without loading it.
    store.lineage = ReadLineage(groups, header)
    adjacency = [entry for (left, right), entries in store.joins.items() for entry in entries
                 if event.relation(left, right, entry[2] - entry[1] - 1)[0] == 'breakpoint_junction']
    clips = defaultdict(list)
    for o in store.observations.values():
        for i, j, assignment in event.clips(o.positions, len(o.sequence)):
            clips['donor' if assignment[0] == 0 else 'acceptor'].append((o, i, j, 'clip'))
    used, junction_cache, witness_groups = {}, {}, {}

    def support(witnesses):
        identities = {store.observations[w['observation']].identity for w in witnesses}
        return dict(**store.cell_umi.support(identities), read_lineage=store.lineage.support(identities),
                    missing_quality_reads=len({store.observations[w['observation']].identity for w in witnesses
                                               if store.observations[w['observation']].missing_qualities}),
                    scope='complete_input_direct_junction_observations', witnesses=witnesses)

    for path in result['paths']:
        positions = _positions(path)
        reported, evidence = _path_result(path['sequence'], positions, set(), store, event, annotated, adjacency, clips,
                                          lambda *args: {}, lambda *args: None, junction_cache=junction_cache)
        path['junctions'] = reported['junctions']
        used.update(evidence)
        path['sequence_evidence']['voting_support_scope'] = 'bounded_discovery'
        path['full_interval_support'] = support(_witnesses(path['sequence'], positions, 0, len(positions),
                                                         path['junctions'], evidence, group_cache=witness_groups))
        for candidate in path['exploratory_orfs']['candidates']:
            candidate['start_evidence']['observation_scope'] = 'bounded_discovery'
            start, end = candidate['query_interval']
            crossed = [path['junctions'][i] for i in candidate['crossed_junctions']]
            candidate['full_interval_support'] = support(_witnesses(path['sequence'], positions, start, end,
                                                                  crossed, evidence, group_cache=witness_groups))
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
    for key, o in sorted(used.items()):
        result['observations'][key] = dict(identity=list(o.identity), records=list(o.records),
                                          reverse_complement=o.reverse_complement,
                                          original_query_interval=list(o.query_interval),
                                          missing_qualities=o.missing_qualities, secondary=o.secondary)
    cited = {rid for o in used.values() for rid in o.records}
    for identity in {o.identity for o in used.values()}:
        for read in store.groups[identity]:
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
                                       eligible_regional_records=available, selected_records=len(selected),
                                       omitted_records=available - len(selected),
                                       complete=available == len(selected) and not result['acquisition']['limitations'])
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
