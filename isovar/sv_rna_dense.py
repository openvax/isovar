"""Disk-backed, bounded discovery followed by complete-input SV support.

Only candidate witnesses are retained in memory during support reporting. The
SQLite spool keeps unrelated background and arbitrarily distant pieces on disk.
"""

from collections import Counter, defaultdict
from collections.abc import Mapping
from contextlib import closing
from hashlib import sha256
from itertools import groupby
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
    _Event, _Model, _Observations, _ineligible, _join_class, _record_gaps, _record_id,
    _path_result, reconstruct_sv_rna,
)
from .sv_rna_orfs import _witnesses
from .sv_rna_relations import LINKAGE_STATUSES


class _DiskGroups(Mapping):
    """Eligible original segment records, retrieved on demand for ancestry."""

    def __init__(self, connection, header):
        self.connection, self.header = connection, header

    def __getitem__(self, identity):
        rows = self.connection.execute(
            'SELECT sam FROM records WHERE rg=? AND name=? AND segment=? AND excluded IS NULL ORDER BY rid',
            identity)
        records = [pysam.AlignedSegment.fromstring(sam, self.header) for sam, in rows]
        if not records:
            raise KeyError(identity)
        return records

    def __iter__(self):
        yield from self.connection.execute(
            'SELECT DISTINCT rg, name, segment FROM records WHERE excluded IS NULL ORDER BY rg, name, segment')

    def __len__(self):
        return self.connection.execute(
            'SELECT count(*) FROM (SELECT DISTINCT rg, name, segment FROM records WHERE excluded IS NULL)').fetchone()[0]

    def templates(self):
        rows = self.connection.execute(
            'SELECT rg, name, segment, sam FROM records WHERE excluded IS NULL ORDER BY rg, name, segment, rid')
        for _, template in groupby(rows, key=lambda row: row[:2]):
            yield [pysam.AlignedSegment.fromstring(row[3], self.header) for row in template]

    def template(self, fragment):
        """All eligible original mates/placements of one consulted template."""
        groups = defaultdict(list)
        for segment, sam in self.connection.execute(
                'SELECT segment, sam FROM records WHERE rg=? AND name=? AND excluded IS NULL ORDER BY segment, rid',
                fragment):
            groups[(*fragment, segment)].append(pysam.AlignedSegment.fromstring(sam, self.header))
        return groups


def _spool(bam, db, collector):
    db.execute('CREATE TABLE records (rid TEXT PRIMARY KEY, rg TEXT, name TEXT, segment INTEGER, sam TEXT, excluded TEXT)')
    scanned, excluded = 0, Counter()
    bam.reset()
    for read in bam.fetch(until_eof=True):
        scanned += 1
        reason = 'unnamed_record' if read.query_name is None else _ineligible(read, collector)
        if reason:
            excluded[reason] += 1
        identity = segment_identity(read)
        sam, rid = read.to_string(), _record_id(read)
        try:
            db.execute('INSERT INTO records VALUES (?, ?, ?, ?, ?, ?)', (rid, *identity, sam, reason))
        except sqlite3.IntegrityError:
            previous = db.execute('SELECT sam FROM records WHERE rid=?', (rid,)).fetchone()[0]
            if previous != sam:
                raise ValueError('Conflicting original SAM records share SV record identity %s' % rid)
        if scanned % 10000 == 0:
            db.commit()
    db.execute('CREATE INDEX segments ON records (rg, name, segment, rid)')
    db.commit()
    return scanned, dict(sorted(excluded.items()))


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
        any(op == 4 for op, _ in read.cigartuples)
        and any(read.reference_name == b.contig and abs(edge - b.position) <= 1
                for b in (event.donor, event.acceptor) for edge in (read.reference_start, read.reference_end))
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


def _select(groups, db, options, collector, spans, event, annotated):
    db.execute('CREATE TABLE discovery (rg TEXT, name TEXT, segment INTEGER, priority INTEGER, signature TEXT, size INTEGER)')
    for template in groups.templates():
        members = defaultdict(list)
        for read in template:
            members[segment_identity(read)].append(read)
        for identity, records in members.items():
            if not any(a < end and start < b for read in records for a, b in read.get_blocks()
                       for start, end in spans.get(read.reference_name, ())):
                continue
            # Ignore read name/labels, but retain library, sequence, placements,
            # orientation, mate, quality and available path structure.
            signature = sha256(repr((identity[0], identity[2], sorted(
                (r.reference_name, r.reference_start, r.flag, r.cigarstring, r.query_sequence,
                 r.qual, r.mapping_quality, r.get_tag('SA') if r.has_tag('SA') else '')
                for r in records))).encode()).hexdigest()
            db.execute('INSERT INTO discovery VALUES (?, ?, ?, ?, ?, ?)',
                       (*identity, _priority(records, collector, event, annotated), signature, len(records)))
    db.commit()
    # Reserve a threshold-sized bundle per sequence/placement class before
    # spending the remaining budget on abundant copies of any one class.
    query = '''SELECT rg, name, segment, size FROM (
        SELECT *, row_number() OVER (PARTITION BY priority, signature ORDER BY rg, name, segment) AS copy
        FROM discovery) ORDER BY priority, (copy - 1) / ?, signature, copy'''
    selected, count = [], 0
    for rg, name, segment, size in db.execute(query, (options['min_alternative_fragments'],)):
        if count + size <= options['max_records']:
            selected.extend(groups[rg, name, segment])
            count += size
    total = db.execute('SELECT coalesce(sum(size), 0) FROM discovery').fetchone()[0]
    return selected, total


def _positions(path):
    positions = [None] * len(path['sequence'])
    for block in path['blocks']:
        bases = range(block['reference_start'], block['reference_end'])
        if block['strand'] == '-':
            bases = reversed(bases)
        for q, base in zip(range(block['query_start'], block['query_end']), bases):
            positions[q] = (block['contig'], base, block['strand'])
    return tuple(positions)


def _support(groups, result, collector, event, annotated, options):
    """Stream complete templates, keeping only candidate junction witnesses."""
    header = groups.header.to_dict()
    joins = {(tuple(j['left']), tuple(j['right'])) for path in result['paths'] for j in path['junctions']
             if j['left'] is not None and j['right'] is not None}
    placements = {(b['contig'], b['reference_start'], b['reference_end'])
                  for path in result['paths'] for b in path['blocks']}
    records, reasons = [], Counter()
    for template in groups.templates() if result['paths'] else ():
        # No per-base path building for unrelated genomic background.
        if not any(read.reference_name == contig and a < high and low < b
                   for read in template for a, b in read.get_blocks() for contig, low, high in placements):
            continue
        store = _Observations(template, collector, options['max_extension_segments'], header,
                              options['sample_id'], options['source'])
        observations = [o for identity in store.groups for o in store.build(identity)]
        reasons.update(store.reasons)
        hit = any((o.positions[i], o.positions[j]) in joins
                  or event.relation(o.positions[i], o.positions[j], j - i - 1)[0] == 'breakpoint_junction'
                  for o in observations for i, j, _ in o.breaks)
        hit |= any(event.clips(o.positions, len(o.sequence)) for o in observations)
        if hit:
            records.extend(template)  # Include visible mates when resolving labels.
    store = _Observations(records, collector, options['max_extension_segments'], header,
                          options['sample_id'], options['source'])
    for identity in store.groups:
        store.build(identity)
    # Parents can be outside the candidates' intervals; keep original producer
    # and ancestry resolution over the complete eligible input, without loading it.
    store.lineage = ReadLineage(groups, header)
    adjacency = [entry for (left, right), entries in store.joins.items() for entry in entries
                 if event.relation(left, right, entry[2] - entry[1] - 1)[0] == 'breakpoint_junction']
    clips = defaultdict(list)
    for o in store.observations.values():
        for i, j, assignment in event.clips(o.positions, len(o.sequence)):
            clips['donor' if assignment[0] == 0 else 'acceptor'].append((o, i, j, 'clip'))
    used = {}

    def support(witnesses):
        identities = {store.observations[w['observation']].identity for w in witnesses}
        return dict(**store.cell_umi.support(identities), read_lineage=store.lineage.support(identities),
                    missing_quality_reads=len({store.observations[w['observation']].identity for w in witnesses
                                               if store.observations[w['observation']].missing_qualities}),
                    scope='complete_input_direct_junction_observations', witnesses=witnesses)

    for path in result['paths']:
        positions = _positions(path)
        reported, evidence = _path_result(path['sequence'], positions, set(), store, event, annotated, adjacency, clips,
                                          lambda *args: {}, lambda *args: None)
        path['junctions'] = reported['junctions']
        used.update(evidence)
        path['sequence_evidence']['voting_support_scope'] = 'bounded_discovery'
        path['full_interval_support'] = support(_witnesses(path['sequence'], positions, 0, len(positions),
                                                         path['junctions'], evidence))
        for candidate in path['exploratory_orfs']['candidates']:
            candidate['start_evidence']['observation_scope'] = 'bounded_discovery'
            start, end = candidate['query_interval']
            crossed = [path['junctions'][i] for i in candidate['crossed_junctions']]
            candidate['full_interval_support'] = support(_witnesses(path['sequence'], positions, start, end,
                                                                  crossed, evidence))
    # Discovery-only context also needs original visible mates/parents. Its
    # voting counts remain explicitly scoped to the bounded discovery set.
    consulted = {tuple(row['identity']) for row in result['cell_umi_evidence']['reads']}
    labels = {}
    for fragment in sorted({identity[:2] for identity in consulted}):
        original = CellUmiEvidence(groups.template(fragment), header, options['sample_id'], options['source'])
        for identity in original.groups:
            if identity in consulted:
                original.segment(identity)
        labels.update((tuple(row['identity']), row) for row in original.evidence()['reads'])
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
        for read in store.groups[o.identity]:
            rid = _record_id(read)
            if rid in o.records:
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
        with closing(sqlite3.connect(str(Path(directory) / 'records.sqlite'))) as db:
            db.execute('PRAGMA temp_store=FILE')
            db.execute('PRAGMA cache_size=-8192')
            scanned, excluded = _spool(bam, db, collector)
            groups = _DiskGroups(db, bam.header)
            selected, available = _select(groups, db, options, collector, spans, event, annotated)
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
                                                  records_scanned=scanned, unique_records=db.execute('SELECT count(*) FROM records').fetchone()[0],
                                                  eligible_segments=len(groups),
                                                  excluded_records=excluded, read_path_notes=support_notes)
            result['discovery'] = dict(policy='event_first_sequence_placement_bundles.v1',
                                       eligible_regional_records=available, selected_records=len(selected),
                                       omitted_records=available - len(selected),
                                       complete=available == len(selected) and not result['acquisition']['limitations'])
            if available > len(selected):
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
