"""Grouped computation keeps original identity, ordering, and interval evidence."""

from dataclasses import replace
import sqlite3

import pysam
import pytest

from isovar.read_collector import ReadCollector
from isovar.read_identity import segment_identity
from isovar.read_lineage import ReadLineage
from isovar.sv_rna import _GeometryCache, _Observations, _segment_path
from isovar.sv_rna_dense import _DiskGroups, _spool
from isovar.sv_rna_orfs import _witnesses
from tests.test_sv_rna import HEADER, record, write_bam
from tests.test_sv_rna_orfs import inputs


@pytest.mark.parametrize('input_format', ['bam', 'sam'])
def test_binary_spool_keeps_all_payloads_mates_and_name_spellings(tmp_path, input_format):
    reads = []
    for name in ('read1', 'read01', 'read001'):
        for flag in (65, 129, 2113, 321):
            read = record(name, '1', 100, '40M', 'A' * 40, flag=flag)
            read.set_tag('CB', name)
            reads.append(read)
    other_group = pysam.AlignedSegment.fromstring(reads[0].to_string(), HEADER)
    other_group.set_tag('RG', 'different-library')
    reads.append(other_group)
    excluded = record('excluded', '1', 100, '40M', 'A' * 40, flag=1024)
    reads += [excluded, reads[0], excluded]
    bam_path = write_bam(tmp_path / 'original.bam', reads)
    if input_format == 'sam':
        bam_path = tmp_path / 'original.sam'
        with pysam.AlignmentFile(bam_path, 'wh', header=HEADER) as output:
            for read in reads:
                output.write(read)
    with pysam.AlignmentFile(bam_path) as bam, sqlite3.connect(':memory:') as db:
        scanned, reasons, unique, segments, path = _spool(bam, db, ReadCollector(), tmp_path)
        assert scanned == len(reads)
        assert reasons == {'duplicate': 2}
        assert unique == len({r.to_string() for r in reads})
        with pysam.AlignmentFile(path) as reader:
            groups = _DiskGroups(db, reader, HEADER, segments)
            wanted = {r.to_string() for r in reads if not r.is_duplicate}
            assert {r.to_string() for group in groups.templates() for r in group} == wanted
            assert len(groups) == len({segment_identity(r) for r in reads if not r.is_duplicate})
            for identity in groups:
                assert identity in groups
                assert {r.to_string() for r in groups[identity]} == {
                    r.to_string() for r in reads if segment_identity(r) == identity and not r.is_duplicate}


def test_grouped_geometry_keeps_the_original_overlap_cap_and_ties():
    reads = [record(name, '1', start, cigar, 'A' * length) for name, start, cigar, length in [
        ('a-short', 100, '10M', 10), ('z-long', 90, '40M', 40),
        ('b-long', 90, '40M', 40), ('a-long', 90, '40M', 40),
        ('partial', 105, '60M', 60), ('terminal-gap', 90, '40M20D', 40)]]
    store = _Observations(reads, ReadCollector(), 3, HEADER.to_dict(), 'sample', 'source')
    window = tuple(('1', q, '+') for q in range(100, 110))
    found = store.overlapping(window, 0, len(window))
    # Equal full coverage, then furthest reference extent, then original ID.
    assert {store.observations[key].identity[1] for key, _ in found} == {
        'terminal-gap', 'a-long', 'b-long'}
    assert store.limited
    reverse = tuple(('1', q, '-') for q in range(109, 99, -1))
    found = store.overlapping(reverse, 0, len(reverse))
    assert {store.observations[key].identity[1] for key, _ in found} == {'a-long', 'b-long', 'terminal-gap'}


def test_shared_interval_comparison_keeps_each_original_orientation_and_query_offset():
    sequence, positions, observations, junction = inputs()
    original = observations.pop('full')
    for key, interval, reverse in [('first', (10, 25), False), ('second', (30, 45), False),
                                   ('reverse', (40, 55), True)]:
        observations[key] = replace(original, key=key, identity=('library', key, 0),
                                    query_interval=interval, reverse_complement=reverse)
    observations['elsewhere'] = replace(original, key='elsewhere',
                                        positions=tuple(('3', q, '+') for q in range(len(sequence))))
    junction['direct_observations'] = sorted(observations)
    cache = {}
    witnesses = _witnesses(sequence, positions, 3, 9, [junction], observations, group_cache=cache)
    assert {w['observation']: w['original_query_interval'] for w in witnesses} == {
        'first': [13, 19], 'second': [33, 39], 'reverse': [46, 52]}
    assert all(w['observation_interval'] == [3, 9] for w in witnesses)
    assert len(next(iter(cache.values()))) == 2  # Same sequence at two distinct placements.
    assert witnesses == _witnesses(sequence, positions, 3, 9, [junction], observations, group_cache=cache)


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('keep_clips', [False, True])
def test_run_geometry_preserves_hard_clips_insertions_gap_kinds_and_qualities(reverse, keep_clips):
    read = record('complex', '1', 100, '5H3S6M2I5M4D4M30N5M2S7H', 'A' * 27, flag=16 if reverse else 0)
    read.query_qualities = list(range(27))
    sequence, positions, breaks, interval, qualities = _segment_path(
        [read], ReadCollector(use_soft_clipped_bases=keep_clips), {}, _GeometryCache())
    bases = list(range(100, 106)) + [None] * 2 + list(range(106, 111)) + list(range(115, 119)) + list(range(149, 154))
    if keep_clips:
        bases = [None] * 3 + bases + [None] * 2
    expected_qualities = list(range(27)) if keep_clips else list(range(3, 25))
    if reverse:
        bases.reverse()
        expected_qualities.reverse()
    assert positions == tuple(None if b is None else ('1', b, '-' if reverse else '+') for b in bases)
    assert sequence == ('T' if reverse else 'A') * len(bases)
    assert qualities == tuple(expected_qualities)
    expected_interval = ((7, 34) if reverse else (5, 32)) if keep_clips else (
        (9, 31) if reverse else (8, 30))
    assert interval == expected_interval
    if keep_clips:
        assert breaks == (((6, 7, 'N'), (10, 11, 'D')) if reverse else ((15, 16, 'D'), (19, 20, 'N')))
    else:
        assert breaks == (((4, 5, 'N'), (8, 9, 'D')) if reverse else ((12, 13, 'D'), (16, 17, 'N')))


def test_inapplicable_ancestry_does_not_fetch_original_records():
    class Unfetched(dict):
        def __getitem__(self, key):
            raise AssertionError('Illumina ancestry does not require a BAM seek')

    lineage = ReadLineage(Unfetched({('library', 'first', mate): None for mate in (64, 128)}),
                          {'RG': [{'ID': 'library', 'PL': 'ILLUMINA'}]})
    support = lineage.support([('library', 'first', 64), ('library', 'first', 128)])
    assert support['reads'] == support['unresolved_reads'] == 2
    assert support['statuses'] == {'unknown_producer': 2}
    assert support['signal_groups'] == 0
