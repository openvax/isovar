"""Nominated breakpoint evidence takes precedence over background depth (#314)."""

from collections import Counter

import pysam
import pytest

from isovar.fusion import FusionBreakpoint
from isovar.sv_rna import collect_sv_records
from .test_sv_rna import record, write_bam


def collect(path, breakpoints=(), regions=(), **kwargs):
    with pysam.AlignmentFile(str(path)) as bam:
        return collect_sv_records(bam, regions, [], breakpoints=breakpoints,
                                  breakpoint_window=10, **kwargs)


def test_breakpoints_precede_deep_background(tmp_path):
    priority = [record('%s%d' % (contig, i), contig, 10000, '10M', 'A' * 10)
                for contig in ('1', '2') for i in range(2)]
    background = [record('background%d' % i, '1', 1000, '10M', 'C' * 10) for i in range(20)]
    path = write_bam(tmp_path / 'deep.bam', background + priority)
    breakpoints = [FusionBreakpoint(c, 10000, '+') for c in ('1', '2')]
    reads, acquisition = collect(path, breakpoints, [('1', 1000, 1100)], max_records=6)
    assert {r.query_name for r in priority} <= {r.query_name for r in reads}
    assert len(reads) == 6 and acquisition['limitations'] == ['record_limit']
    # Neither records nor qualities/flags/tags are rewritten by prioritization.
    assert {r.to_string() for r in priority} <= {r.to_string() for r in reads}


@pytest.mark.parametrize('second_contig,second_position', [('2', 1000), ('1', 1010)])
def test_deep_breakpoint_windows_share_the_budget(tmp_path, second_contig, second_position):
    # In the same-contig case the windows overlap: merging them would starve
    # the second breakpoint behind the deep, coordinate-sorted first window.
    reads = [record('left%d' % i, '1', 991, '1M', 'A') for i in range(20)]
    reads += [record('right%d' % i, second_contig, second_position + 1, '1M', 'C') for i in range(20)]
    path = write_bam(tmp_path / 'balanced.bam', reads)
    breakpoints = [FusionBreakpoint('1', 1000, '+'), FusionBreakpoint(second_contig, second_position, '-')]
    selected, acquisition = collect(path, breakpoints, max_records=5)
    assert Counter(r.query_name.rstrip('0123456789') for r in selected) == dict(left=3, right=2)
    assert acquisition['queries'] == 2 and acquisition['limitations'] == ['record_limit']
    again, _ = collect(path, reversed(breakpoints), max_records=5)
    assert [r.to_string() for r in selected] == [r.to_string() for r in again]


@pytest.mark.parametrize('hint', ['SA', 'mate'])
def test_priority_links_inside_unvisited_background_are_recovered_first(tmp_path, hint):
    first = record('linked', '1', 1000, '10M', 'A' * 10, flag=65 if hint == 'mate' else 0)
    second = record('linked', '2', 8000, '10M', 'C' * 10, flag=129 if hint == 'mate' else 2048)
    third = record('linked', '1', 15000, '10M', 'G' * 10, flag=2177 if hint == 'mate' else 2048)
    if hint == 'SA':
        first.set_tag('SA', '2,8001,+,10M,60,0;')
    else:
        first.next_reference_id, first.next_reference_start = second.reference_id, second.reference_start
        for read in (second, third):
            read.next_reference_id, read.next_reference_start = first.reference_id, first.reference_start
    second.set_tag('SA', '1,15001,+,10M,60,0;')
    other_side = record('other-side', '2', 12000, '10M', 'T' * 10)
    background = [record('background%d' % i, '1', 5000, '10M', 'C' * 10) for i in range(20)]
    # A same-name record from another library must not enter via a link.
    collision = record('linked', '2', 8000, '10M', 'C' * 10, group='b')
    expected = [first, second, third, other_side]
    path = write_bam(tmp_path / 'links.bam', background + expected + [collision])
    breakpoints = [FusionBreakpoint('1', 1000, '+'), FusionBreakpoint('2', 12000, '+')]
    reads, acquisition = collect(path, breakpoints, [('1', 5000, 5100), ('2', 8000, 8100)], max_records=4)
    assert {r.to_string() for r in reads} == {r.to_string() for r in expected}
    assert [(h['contig'], h['start']) for h in acquisition['hop_queries'] if h['fetched']] == [
        ('2', 8000), ('1', 15000)]
    assert 'record_limit' in acquisition['limitations']


def test_overlapping_windows_deduplicate_without_losing_records(tmp_path):
    reads = [record('shared', '1', 1000, '20M', 'A' * 20),
             record('right', '1', 1015, '10M', 'C' * 10),
             record('background', '2', 5000, '10M', 'G' * 10)]
    path = write_bam(tmp_path / 'overlap.bam', reads + reads[:1])
    breakpoints = [FusionBreakpoint('1', 1000, '+'), FusionBreakpoint('1', 1010, '-')]
    selected, acquisition = collect(path, breakpoints, [('1', 990, 1030), ('2', 5000, 5100)], max_records=3)
    assert Counter(r.to_string() for r in selected) == Counter(r.to_string() for r in reads)
    assert acquisition['limitations'] == []  # Exactly full, with nothing omitted.


def test_query_limit_reports_unsearched_breakpoint(tmp_path):
    reads = [record(c, c, 1000, '10M', 'A' * 10) for c in ('1', '2')]
    path = write_bam(tmp_path / 'limited.bam', reads)
    selected, acquisition = collect(path, [FusionBreakpoint(c, 1000, '+') for c in ('1', '2')], max_queries=1)
    assert [r.query_name for r in selected] == ['1']
    assert acquisition['queries'] == 1 and acquisition['limitations'] == ['query_limit']
    assert acquisition['region_queries'] == [
        dict(contig='1', start=990, end=1010, priority=True, fetched=True, complete=True),
        dict(contig='2', start=990, end=1010, priority=True, fetched=False, complete=False)]


def test_exhausted_window_returns_its_remaining_budget_to_the_other(tmp_path):
    reads = [record('sparse', '1', 1000, '10M', 'A' * 10)]
    reads += [record('deep%d' % i, '2', 1000, '10M', 'C' * 10) for i in range(20)]
    path = write_bam(tmp_path / 'asymmetric.bam', reads)
    selected, acquisition = collect(path, [FusionBreakpoint(c, 1000, '+') for c in ('1', '2')], max_records=6)
    assert len(selected) == 6 and sum(r.reference_name == '2' for r in selected) == 5
    assert [q['complete'] for q in acquisition['region_queries']] == [True, False]
    assert acquisition['limitations'] == ['record_limit']


def test_query_cap_prevents_link_and_background_fetches(tmp_path):
    first = record('linked', '1', 1000, '10M', 'A' * 10)
    first.set_tag('SA', '2,8001,+,10M,60,0;')
    other_side = record('other', '2', 1000, '10M', 'C' * 10)
    linked = record('linked', '2', 8000, '10M', 'G' * 10, flag=2048)
    path = write_bam(tmp_path / 'query-cap.bam', [first, other_side, linked])
    selected, acquisition = collect(path, [FusionBreakpoint(c, 1000, '+') for c in ('1', '2')],
                                    [('2', 7900, 8100)], max_queries=2)
    assert {r.to_string() for r in selected} == {first.to_string(), other_side.to_string()}
    assert acquisition['queries'] == 2 and acquisition['limitations'] == ['query_limit']
    assert acquisition['hop_queries'] == [
        dict(contig='2', start=8000, end=8001, fragments=1, fetched=False, complete=False)]
    assert len(acquisition['region_queries']) == 2 and all(q['complete'] for q in acquisition['region_queries'])
