"""Complete-input support remains separate from bounded hypothesis discovery."""

from dataclasses import asdict
import json
from pathlib import Path
import weakref

import pysam
import pytest

from isovar import export_sv_rna_orfs
from isovar.cli.isovar_sv_rna import run as cli_run
from isovar.read_collector import ReadCollector
from isovar.fusion import FusionBreakpoint, FusionReference
from isovar.sv_rna import reconstruct_sv_rna
from tests.test_sv_rna import Scenario, aligned, record, write_bam


@pytest.mark.parametrize('strand', ['+', '-'])
def test_frozen_span_index_keeps_nested_and_repeated_genomic_blocks(strand):
    from isovar.sv_rna import _PathSpans, _Positions

    positions = []
    for low, high in [(1000, 1010), (100, 150), (105, 110), (100, 150)]:
        bases = range(low, high) if strand == '+' else range(high - 1, low - 1, -1)
        positions.extend(('1', base, strand) for base in bases)
    span = _PathSpans('A' * len(positions), _Positions(positions))
    assert {start for start, _, _ in span.overlapping_runs('1', strand, 108, 109)} == {10, 60, 65}
    assert span.overlapping_runs('1', strand, 150, 1000) == []
    assert not span.overlapping_runs('2', strand, 100, 1001)
    assert not span.overlapping_runs('1', '+' if strand == '-' else '-', 100, 1001)


def test_fragmented_candidate_span_queries_skip_disjoint_blocks():
    from types import SimpleNamespace
    from isovar.sv_rna import _PathSpans, _Positions

    positions = _Positions([('1', q + q // 10, '+') for q in range(4000)])
    span = _PathSpans('A' * 4000, positions)
    assert len(span.overlapping_runs('1', '+', 110, 120)) == 1
    # Placed-base alleles do not alter the co-observed interval; exact sequence
    # agreement remains the separate full-interval witness calculation.
    observation = SimpleNamespace(key='original', positions=positions, sequence='C' * 4000)
    assert span(observation) == (0, 3999)


@pytest.mark.parametrize('strand', ['+', '-'])
def test_span_rejects_incompatible_query_offsets_before_remaining_blocks(strand):
    from types import SimpleNamespace
    from isovar.sv_rna import _PathSpans, _Positions

    bases = [q + q // 10 for q in range(4000)]
    if strand == '-':
        bases.reverse()
    positions = _Positions([('1', q, strand) for q in bases])
    span = _PathSpans('A' * 4000, positions)
    original = span.overlapping_runs
    visited = []

    def record(*args):
        visited.append(args)
        return original(*args)

    span.overlapping_runs = record
    # One inserted query base changes the offset after the first block. The
    # old equal-gap-length rule rejects it regardless of all later blocks.
    observed = _Positions(positions[:10] + (None,) + positions[10:])
    observation = SimpleNamespace(key='indel', positions=observed, sequence='A' * 4001)
    assert span(observation) == ()
    assert len(visited) == 2


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('batch_bases', [1, 400, 1700])
def test_stream_batches_preserve_whole_results_and_release_observations(tmp_path, monkeypatch, strand, batch_bases):
    from isovar import sv_rna_dense, sv_rna_stream

    scenario = Scenario(strand)
    reads = dense_reads(scenario, 23)
    # Shared cell labels across batches must count once, while individual
    # qualities and an allele difference still belong to their original reads.
    for read in reads:
        if read.query_name.startswith('event-'):
            read.set_tag('CB', 'shared-cell')
            read.set_tag('UB', 'shared-umi')
    for read in reads[-4:-2]:
        quality = read.query_qualities
        read.query_sequence = 'A' + read.query_sequence[1:]
        read.query_qualities = quality
    bam = write_bam(tmp_path / 'batches.bam', reads)
    expected = scenario.run(bam, dense_support=True, max_records=4, assemble=False)
    monkeypatch.setattr(sv_rna_stream, '_SUPPORT_BATCH_BASES', batch_bases)
    live, sizes = weakref.WeakSet(), []
    build, consume = sv_rna_dense._build_segment, sv_rna_stream._SupportStream.consume

    def tracked(*args):
        built = build(*args)
        live.update(built)
        # One indivisible template may exceed the base budget; the next
        # template may be built before the previous batch is consumed.
        assert sum(len(o.sequence) for o in live) <= max(batch_bases, 400) + 400
        return built

    def counted(self, observations):
        sizes.append(sum(len(o.sequence) for o in observations))
        return consume(self, observations)

    monkeypatch.setattr(sv_rna_dense, '_build_segment', tracked)
    monkeypatch.setattr(sv_rna_stream._SupportStream, 'consume', counted)
    observed = scenario.run(bam, dense_support=True, max_records=4, assemble=False)
    assert observed == expected
    assert not live
    assert len(sizes) > 1
    assert max(sizes) <= max(batch_bases, 400)
    assert candidate(observed)['junctions'][0]['direct_support']['cells'] == 1


def test_splice_ambiguous_event_joins_survive_abundant_breakpoint_clips(tmp_path):
    sequence = "ATG" + "GCT" * 29
    reads = [record("clip-%02d" % i, "1", 10, "90M10S", "C" * 90 + "A" * 10)
             for i in range(20)]
    for i in range(7):
        read = record("ambiguous-%02d" % i, "1", 60, "40M200N50M", sequence, qualities=False)
        read.set_tag("CB", "cell-%d" % i)
        reads.append(read)
    with pysam.AlignmentFile(str(write_bam(tmp_path / "ambiguous.bam", reads))) as bam:
        result = reconstruct_sv_rna(
            bam, event_id="ordinary-splice-compatible", reference_name="test",
            donor=FusionBreakpoint("1", 100, "+"), acceptor=FusionBreakpoint("1", 200, "+"),
            regions=[], references=[], sample_id="sample", source="pinned-original",
            event_provenance={"test": "dense-splice-priority"}, dense_support=True, assemble=False, max_records=2,
            read_collector=ReadCollector(use_soft_clipped_bases=True))
    path, = [p for p in result["paths"] if p["sequence"] == sequence]
    assert path["junctions"][0]["relation"] == "splice_ambiguous_event_junction"
    assert path["full_interval_support"]["fragments"] == 7
    assert path["full_interval_support"]["cells"] == 7
    assert path["full_interval_support"]["missing_quality_reads"] == 7
    assert result["support_acquisition"]["records_scanned"] == len(reads)
    assert result["discovery"]["selected_records"] == 2


def dense_reads(scenario, n=13):
    reads = [record('background-%04d' % i, '1', 2000, '90M', scenario.donor_ref.sequence[100:190])
             for i in range(250)]
    for i in range(n):
        pieces = aligned('event-%03d' % i, scenario.sequence[200:300], scenario.positions[200:300],
                         qualities=i != n - 1)
        for read in pieces:
            read.set_tag('CB', 'cell-%d' % i)
            read.set_tag('UB', 'umi-%d' % i)
        reads.extend(pieces)
    return reads


def candidate(result):
    path, = result['paths']
    return path


def test_late_rare_event_and_complete_support_do_not_depend_on_discovery_cap_or_order(tmp_path):
    scenario = Scenario()
    reads = dense_reads(scenario)
    bam = tmp_path / 'dense.bam'
    write_bam(bam, reads)
    assert scenario.run(bam, max_records=2, assemble=False)['status'] == 'no_candidate_paths'
    results = []
    for cap, ordering in [(2, reads), (4, list(reversed(reads))), (10, reads)]:
        write_bam(bam, ordering)
        result = scenario.run(bam, max_records=cap, assemble=False, dense_support=True, scratch_dir=tmp_path)
        path = candidate(result)
        assert path['sequence'] == scenario.sequence[200:300]
        assert path['full_interval_support']['reads'] == path['full_interval_support']['fragments'] == 13
        assert path['full_interval_support']['cells'] == path['full_interval_support']['umis'] == 13
        assert path['full_interval_support']['missing_quality_reads'] == 1
        assert path['junctions'][0]['direct_support']['fragments'] == 13
        assert len(path['full_interval_support']['witnesses']) == 13
        assert path['sequence_evidence']['voting_support']['fragments'] <= cap // 2
        assert path['sequence_evidence']['voting_support_scope'] == 'bounded_discovery'
        assert 'discovery_record_limit' in result['limitations']
        assert result['discovery']['omitted_records'] > 200
        assert result['support_acquisition']['complete']
        assert result['support_acquisition']['records_scanned'] == len(reads)
        assert all('background' not in sam for sam in result['original_records'].values())
        assert {sam for sam in result['original_records'].values()} == {
            read.to_string() for read in reads if read.query_name.startswith('event-')}
        assert sorted(row['cell_barcode'] for row in result['cell_umi_evidence']['reads']) == [
            'cell-%d' % i for i in sorted(range(13), key=lambda i: 'cell-%d' % i)]
        results.append(path['full_interval_support'])
    assert results[0] == results[1] == results[2]
    # All spool/BAM files are cleaned on success.
    assert not list(tmp_path.glob('isovar-sv-*'))


def test_quality_and_label_variation_share_discovery_work_but_keep_all_support(tmp_path, monkeypatch):
    from isovar import sv_rna_dense

    scenario = Scenario()
    reads = dense_reads(scenario, 30)
    for i, read in enumerate(reads):
        if read.query_qualities is not None:
            read.query_qualities = [1 + i % 40] * len(read.query_sequence)
    calls = []
    priority = sv_rna_dense._priority

    def counted(records, *args):
        calls.append(records[0].query_name)
        return priority(records, *args)

    monkeypatch.setattr(sv_rna_dense, '_priority', counted)
    result = scenario.run(write_bam(tmp_path / 'quality.bam', reads), dense_support=True, assemble=False)
    assert len(calls) == 2  # Background and event, independent of quality/name/labels.
    assert result['discovery']['eligible_regional_groups'] == 2
    assert result['discovery']['selected_records'] == len(reads)
    assert result['discovery']['represented_records'] == len(reads)
    assert result['discovery']['groups_complete']
    assert result['discovery']['complete']
    path = candidate(result)
    assert path['full_interval_support']['fragments'] == 30
    assert path['full_interval_support']['cells'] == path['full_interval_support']['umis'] == 30
    assert path['full_interval_support']['missing_quality_reads'] == 1
    assert set(result['original_records'].values()) == {
        read.to_string() for read in reads if read.query_name.startswith('event-')}


def test_grouping_preserves_major_minor_allele_proportions(tmp_path):
    scenario = Scenario()
    major = scenario.sequence[200:300]
    minor = major[:20] + ('A' if major[20] != 'A' else 'C') + major[21:]
    reads = [r for name, count, sequence in [('major', 20, major), ('minor', 2, minor)]
             for i in range(count) for r in aligned('%s-%d' % (name, i), sequence, scenario.positions[200:300])]
    result = scenario.run(write_bam(tmp_path / 'proportions.bam', reads), dense_support=True, assemble=False)
    path = candidate(result)
    assert path['sequence'] == major
    assert path['full_interval_support']['fragments'] == 20
    assert path['junctions'][0]['direct_support']['fragments'] == 22


@pytest.mark.parametrize('reverse', [False, True])
def test_structural_cache_shares_positions_but_preserves_bases_and_quality(reverse, monkeypatch):
    from collections import Counter
    from isovar import sv_rna
    from isovar.read_identity import segment_identity

    collector, cache = ReadCollector(use_soft_clipped_bases=True), {}
    reads = [record('read-%d' % i, '1', 100, '5S40M100N50M5S', 'A' * 100,
                    flag=16 if reverse else 0, qualities=i != 2) for i in range(3)]
    reads[1].query_sequence = 'A' * 60 + 'C' + 'A' * 39
    reads[1].query_qualities = [7] * 100
    expected = [sv_rna._build_segment(segment_identity(r), [r], collector, Counter(), {}) for r in reads]
    calls, placements = [], sv_rna._single_record_geometry

    def counted(*args):
        calls.append(1)
        return placements(*args)

    monkeypatch.setattr(sv_rna, '_single_record_geometry', counted)
    observed = [sv_rna._build_segment(segment_identity(r), [r], collector, Counter(), {}, cache) for r in reads]
    assert [[asdict(o) for o in group] for group in observed] == [[asdict(o) for o in group] for group in expected]
    assert len(calls) == 1
    assert observed[0][0].positions is observed[1][0].positions
    assert observed[0][1].positions is observed[2][1].positions


def test_support_screen_skips_overlapping_background_without_changing_evidence(tmp_path, monkeypatch):
    from isovar import sv_rna_dense

    scenario = Scenario()
    reads = dense_reads(scenario, 5)
    bam = write_bam(tmp_path / 'background.bam', reads)
    screened = scenario.run(bam, dense_support=True, assemble=False)
    monkeypatch.setattr(sv_rna_dense, '_may_support', lambda *args: True)
    exhaustive = scenario.run(bam, dense_support=True, assemble=False)
    assert screened == exhaustive


@pytest.mark.parametrize('dense', [False, True])
@pytest.mark.parametrize('leading', [False, True])
@pytest.mark.parametrize('operation', ['D', 'N'])
def test_clip_discovery_and_support_use_actual_aligned_edges(tmp_path, dense, leading, operation):
    scenario = Scenario()
    contig, start, cigar = ('2', 5995, '50S5%s50M' % operation) if leading else (
        '1', 2100, '50M5%s50S' % operation)
    reads = [record('terminal-%d' % i, contig, start, cigar, scenario.sequence[200:300]) for i in range(3)]
    options = dict(dense_support=dense, assemble=False, read_collector=ReadCollector(use_soft_clipped_bases=True))
    result = scenario.run(write_bam(tmp_path / 'terminal.bam', reads), *scenario.exact, **options)
    path = candidate(result)
    assert path['junctions'][0]['relation'] == 'breakpoint_clip_partner_unplaced'
    assert path['junctions'][0]['direct_support']['fragments'] == 3
    if dense:
        assert path['full_interval_support']['fragments'] == 3
    reads.append(record('ordinary', contig, 6000 if leading else 2100,
                        '50S50M' if leading else '50M50S', scenario.sequence[200:300]))
    result = scenario.run(write_bam(tmp_path / 'mixed.bam', reads), *scenario.exact, **options)
    assert candidate(result)['junctions'][0]['direct_support']['fragments'] == 4


@pytest.mark.parametrize('change', ['base', 'insertion', 'coverage', 'placement', 'supplementary'])
def test_grouping_keeps_full_alleles_and_actual_alignment_paths(change):
    from isovar.sv_rna_dense import _group_signatures

    original = record('a', '1', 100, '90M', 'A' * 90)
    other = record('b', '1', 100, '90M', 'A' * 90)
    records = [other]
    if change == 'base':
        # Even a distal difference, outside a breakpoint fingerprint, stays distinct.
        other.query_sequence = 'A' * 89 + 'C'
    elif change == 'insertion':
        other.cigarstring = '45M1I45M'
        other.query_sequence = 'A' * 45 + 'C' + 'A' * 45
    elif change == 'coverage':
        other.cigarstring = '80M'
        other.query_sequence = 'A' * 80
    elif change == 'placement':
        other.reference_start += 1
    else:
        records.append(record('b', '2', 100, '90M', 'A' * 90, flag=2048))
    structure, allele = _group_signatures(('a', 'a', 0), [original])
    other_structure, other_allele = _group_signatures(('a', 'b', 0), records)
    assert (structure == other_structure) == (change == 'base')
    assert allele != other_allele


@pytest.mark.parametrize('strand', ['+', '-'])
def test_original_query_witnesses_and_orf_export_use_all_records(tmp_path, strand):
    scenario = Scenario(strand)
    reads = []
    for i in range(9):
        reads += aligned('long-%d' % i, scenario.sequence, scenario.positions, qualities=False)
    result = scenario.run(write_bam(tmp_path / 'long.bam', reads), dense_support=True,
                          max_records=4, assemble=False)
    path = candidate(result)
    assert path['full_interval_support']['fragments'] == 9
    assert path['full_interval_support']['missing_quality_reads'] == 9
    assert path['exploratory_orfs']['candidates']
    for orf in path['exploratory_orfs']['candidates']:
        support = orf['full_interval_support']
        assert support['fragments'] == support['missing_quality_reads'] == 9
        assert len(support['witnesses']) == 9
        start, end = orf['query_interval']
        assert all(w['original_query_interval'] == [start, end] for w in support['witnesses'])
    for orf in export_sv_rna_orfs(result)['candidates']:
        assert orf['rna_support']['fragments'] == 9
        assert orf['rna_support']['missing_quality_reads'] == 9


def test_filters_and_conflicting_mate_labels_are_preserved(tmp_path):
    scenario = Scenario()
    reads = []
    for i in range(5):
        pieces = aligned('template-%d' % i, scenario.sequence[200:300], scenario.positions[200:300], flag=65)
        for read in pieces:
            read.set_tag('CB', 'cell-A')
            read.set_tag('UB', 'umi-%d' % i)
        reads += pieces
        mate = record('template-%d' % i, '2', 15000, '80M', 'A' * 80, flag=129)
        mate.set_tag('CB', 'cell-B')
        mate.set_tag('UB', 'umi-%d' % i)
        reads.append(mate)
    # Duplicate and low-MAPQ event observations remain excluded in both passes.
    reads += aligned('duplicate', scenario.sequence[200:300], scenario.positions[200:300], flag=1024)
    low = aligned('low', scenario.sequence[200:300], scenario.positions[200:300])
    for read in low:
        read.mapping_quality = 0
    reads += low
    result = scenario.run(write_bam(tmp_path / 'mates.bam', reads), dense_support=True,
                          max_records=4, assemble=False)
    support = candidate(result)['full_interval_support']
    assert support['reads'] == support['fragments'] == 5
    assert support['cells'] == support['umis'] == 0
    assert support['label_statuses'] == {'conflicting_template_labels': 5}
    assert result['support_acquisition']['excluded_records']
    assert all('duplicate' not in sam and not sam.startswith('low\t') for sam in result['original_records'].values())


def test_distinct_sequences_are_not_combined_into_one_supported_hypothesis(tmp_path):
    scenario = Scenario()
    original = scenario.sequence[200:300]
    changed = original[:20] + ('A' if original[20] != 'A' else 'C') + original[21:]
    reads = [read for i in range(5) for read in aligned('a-%d' % i, original, scenario.positions[200:300])]
    reads += [read for i in range(4) for read in aligned('b-%d' % i, changed, scenario.positions[200:300])]
    result = scenario.run(write_bam(tmp_path / 'alternatives.bam', reads), dense_support=True,
                          max_records=8, assemble=False)
    assert {path['sequence']: path['full_interval_support']['fragments'] for path in result['paths']} == {
        original: 5, changed: 4}
    assert all(path['junctions'][0]['direct_support']['fragments'] == 9 for path in result['paths'])
    assert {path['sequence']: path['event_local_support']['fragments'] for path in result['paths']} == {
        original: 5, changed: 4}
    assert {path['event_junction_window_support']['fragments'] for path in result['paths']} == {9}
    assert {path['event_junction_window_support']['specificity'] for path in result['paths']} == {
        'shared_or_unresolved'}
    assert all(path['event_local_support']['specificity'] == 'unique_among_retained_paths'
               for path in result['paths'])


def test_local_window_requires_the_branch_beyond_a_shared_junction():
    from isovar.sv_rna_stream import _local_windows

    shared = 'A' * 120
    first = shared[:90] + 'C' + shared[91:]
    second = shared[:90] + 'G' + shared[91:]
    paths = [{'sequence': first}, {'sequence': second}]
    positions = [(None,) * len(first), (None,) * len(second)]
    queries = [[{'key': 'event', 'relation': 'breakpoint_junction',
                 'query_interval': [60, 60]}] for _ in paths]
    windows = _local_windows(paths, positions, queries)
    assert {row['flank_bases'] for row in windows.values()} == {48}
    assert {row['specificity'] for row in windows.values()} == {'unique_among_retained_paths'}
    assert {tuple(row['query_interval']) for row in windows.values()} == {(12, 108)}


def test_local_window_marks_distant_branch_unresolved():
    from isovar.sv_rna_stream import _local_windows

    first = 'A' * 450
    second = first[:400] + 'C' + first[401:]
    paths = [{'sequence': first}, {'sequence': second}]
    positions = [(None,) * len(first), (None,) * len(second)]
    queries = [[{'key': 'event', 'relation': 'breakpoint_junction',
                 'query_interval': [200, 200]}] for _ in paths]
    windows = _local_windows(paths, positions, queries)
    assert {row['flank_bases'] for row in windows.values()} == {150}
    assert {row['specificity'] for row in windows.values()} == {'shared_or_unresolved'}


def test_pinned_hcc_and_colo_sv_local_evidence():
    from reports.hcc1395_colo829.build_sv_transcript_evidence import build

    rows, _, _ = build()
    assert len(rows) == 7
    chr9 = [row for row in rows if row['call_id'].startswith('chr9_')]
    assert len(chr9) == 3
    assert {row['direct_junction_fragments'] for row in chr9} == {7}
    assert sorted(row['local_exact_fragments'] for row in chr9) == [1, 1, 3]
    assert sorted(row['junction_12_flank_exact_fragments'] for row in chr9) == [2, 2, 3]
    assert {row['junction_12_flank_specificity'] for row in chr9} == {'shared_or_unresolved'}
    assert all(row['local_specificity'] == 'unique_among_retained_paths' for row in chr9)
    sv3, = [row for row in rows if row['sample'] == 'COLO829']
    assert (sv3['direct_junction_fragments'], sv3['direct_junction_fragments_mapq20'],
            sv3['local_exact_fragments'], sv3['full_path_exact_fragments']) == (3, 0, 3, 0)
    assert sv3['junction_inserted'] == 'CTATAGTAGAACTTCTTTT'
    assert sv3['reference_donor_after'] != sv3['junction_inserted']


def test_no_candidates_and_failed_scans_never_claim_negative_evidence(tmp_path):
    scenario = Scenario()
    reads = [record('ordinary', '1', 2000, '90M', scenario.donor_ref.sequence[100:190])]
    result = scenario.run(write_bam(tmp_path / 'ordinary.bam', reads), dense_support=True, max_records=2)
    assert result['status'] == 'no_candidate_paths'
    assert result['discovery']['complete'] and result['support_acquisition']['complete']
    with pytest.raises(ValueError, match='scratch_dir requires'):
        scenario.run(tmp_path / 'ordinary.bam', scratch_dir=tmp_path)
    with pytest.raises(ValueError, match='positive integers'):
        scenario.run(tmp_path / 'ordinary.bam', dense_support=True, max_records=0)
    with pytest.raises(OSError):
        scenario.run(tmp_path / 'ordinary.bam', dense_support=True, scratch_dir=tmp_path / 'absent')


def test_cli_dense_support_and_scratch_directory(tmp_path):
    scenario = Scenario()
    bam = write_bam(tmp_path / 'cli.bam', dense_reads(scenario, 5))
    inputs = dict(event_id='D--A', reference_name='test', sample_id='sample',
                  donor=asdict(scenario.donor), acceptor=asdict(scenario.acceptor), regions=[],
                  references=[asdict(scenario.donor_ref), asdict(scenario.acceptor_ref)],
                  event_provenance=dict(caller='synthetic DNA call'))
    event, output = tmp_path / 'event.json', tmp_path / 'output.json'
    event.write_text(json.dumps(inputs))
    cli_run(['--input', str(event), '--bam', str(bam), '--output', str(output), '--dense-support',
             '--sv-scratch-dir', str(tmp_path), '--max-records', '2', '--no-assembly'])
    result = json.loads(output.read_text())
    assert candidate(result)['full_interval_support']['fragments'] == 5
    assert result['parameters']['dense_support'] is True


def test_duplicate_fetches_and_alternative_placements_do_not_create_support(tmp_path):
    scenario = Scenario()
    reads = dense_reads(scenario, 4)
    original = scenario.run(write_bam(tmp_path / 'original.bam', reads), max_records=4,
                            assemble=False, dense_support=True)
    # Repeated identical records are one original alignment, not new molecules.
    duplicate = scenario.run(write_bam(tmp_path / 'duplicate.bam', reads + reads), max_records=4,
                             assemble=False, dense_support=True)
    assert candidate(original)['full_interval_support'] == candidate(duplicate)['full_interval_support']
    assert duplicate['support_acquisition']['records_scanned'] == 2 * len(reads)
    assert duplicate['support_acquisition']['unique_records'] == len(reads)

    alternative = record('event-000', '1', 12000, '100M', scenario.sequence[200:300], flag=256)
    ambiguous = scenario.run(write_bam(tmp_path / 'ambiguous.bam', reads + [alternative]), max_records=8,
                             assemble=False, dense_support=True)
    assert candidate(ambiguous)['full_interval_support']['fragments'] == 4


def test_interrupted_scan_leaves_no_partial_result(tmp_path, monkeypatch):
    from isovar import sv_rna_dense
    scenario = Scenario()
    reads = dense_reads(scenario, 3)
    bam = write_bam(tmp_path / 'interrupted.bam', reads)
    original = sv_rna_dense._spool

    def interrupted(*args):
        original(*args)
        raise OSError('injected scan interruption')

    monkeypatch.setattr(sv_rna_dense, '_spool', interrupted)
    with pytest.raises(OSError, match='injected scan interruption'):
        scenario.run(bam, dense_support=True, scratch_dir=tmp_path)
    assert not list(tmp_path.glob('isovar-sv-*'))


@pytest.mark.parametrize('dense', [False, True])
def test_distinct_original_payloads_preserve_ambiguous_labels_without_extra_fragments(tmp_path, dense):
    scenario = Scenario()
    reads = dense_reads(scenario, 3)
    conflict = pysam.AlignedSegment.fromstring(reads[-1].to_string(), reads[-1].header)
    conflict.set_tag('CB', 'another-cell')
    bam = write_bam(tmp_path / 'conflicting.bam', reads + [conflict])
    result = scenario.run(bam, dense_support=dense, assemble=False, **({'scratch_dir': tmp_path} if dense else {}))
    path = candidate(result)
    assert path['junctions'][0]['direct_support']['fragments'] == 3
    assert path['junctions'][0]['direct_support']['cells'] == 2
    assert path['junctions'][0]['direct_support']['label_statuses']['conflicting_tags'] == 1
    assert {read.to_string() for read in reads if read.query_name.startswith('event-')} | {
        conflict.to_string()} <= set(result['original_records'].values())
    if dense:
        assert path['full_interval_support']['fragments'] == 3
        assert path['full_interval_support']['cells'] == 2
        assert result['support_acquisition']['unique_records'] == len(reads) + 1
    assert not list(tmp_path.glob('isovar-sv-*'))


@pytest.mark.parametrize('field', ['mapping_quality', 'next_reference_id', 'next_reference_start',
                                 'template_length', 'query_sequence', 'query_qualities', 'CB', 'UB', 'custom'])
def test_indexed_collection_keeps_distinct_sam_payloads_and_deduplicates_repeat_fetches(tmp_path, field):
    from isovar.sv_rna import _record_id, collect_sv_records
    original = record('same-segment', '1', 100, '90M', 'A' * 90)
    different = pysam.AlignedSegment.fromstring(original.to_string(), original.header)
    if field == 'query_sequence':
        different.query_sequence = 'C' + 'A' * 89
    elif field == 'query_qualities':
        different.query_qualities = [17] * 90
    elif field in ('CB', 'UB', 'custom'):
        different.set_tag('ZZ' if field == 'custom' else field, 'changed')
    else:
        setattr(different, field, 1 if field == 'next_reference_id' else 123)
    assert _record_id(original) != _record_id(different)
    bam_path = write_bam(tmp_path / 'distinct.bam', [original, different, original])
    with pysam.AlignmentFile(bam_path) as bam:
        records, _ = collect_sv_records(bam, [('1', 100, 150), ('1', 110, 190)], [])
    assert {read.to_string() for read in records} == {original.to_string(), different.to_string()}
    assert len(records) == 2


def test_original_short_read_fixture_keeps_pinned_translation_and_support(tmp_path):
    from tests.test_sv_rna import corpus_bam, run_corpus
    from isovar.read_end_inference import Adapter, ReadEndProfile
    _, bam, inputs = corpus_bam(tmp_path, 'ATP5MG--KMT2A')
    original = run_corpus(bam, inputs)
    result = run_corpus(bam, inputs, dense_support=True, max_records=10)
    path = candidate(result)
    assert path['sequence'] == candidate(original)['sequence']
    assert [t['amino_acids'] for t in path['translations']] == ['MAQFVRNLVEKTPALVNG']
    assert path['junctions'][0]['direct_support']['fragments'] == 2
    assert path['full_interval_support']['fragments'] <= 2
    assert result['support_acquisition']['records_complete']
    assert export_sv_rna_orfs(result)['support_acquisition'] == result['support_acquisition']
    profile = ReadEndProfile('TruSeq', '1', (Adapter('TruSeq', 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'),))
    trimmed = run_corpus(bam, inputs, dense_support=True, max_records=10, read_collector=ReadCollector(
        use_soft_clipped_bases=True, read_end_profile=profile, trim_adapters=True))
    assert [p['sequence'] for p in trimmed['paths']] == [path['sequence']]


def test_iterable_parameters_and_input_scoped_signal_parents(tmp_path):
    from tests.test_sv_rna import HEADER
    scenario = Scenario()
    header = HEADER.to_dict()
    header['RG'][0].update(PL='ONT', PG='dorado', LB='library', SM='sample')
    header['PG'] = [dict(ID='dorado', PN='dorado', CL='dorado basecaller model input')]
    header = pysam.AlignmentHeader.from_dict(header)
    reads = []
    for i in range(3):
        pieces = aligned('child-%d' % i, scenario.sequence[200:300], scenario.positions[200:300])
        for read in pieces:
            read.set_tag('pi', 'parent')
            read.set_tag('dx', 0)
            read.set_tag('sp', i * 100)
        reads += pieces
    # The visible parent has invalid ancestry, beyond every candidate interval.
    parent = record('parent', '2', 18000, '50M', 'A' * 50)
    parent.set_tag('dx', 1)
    reads.append(parent)
    reads = [pysam.AlignedSegment.fromstring(r.to_string(), header) for r in reads]
    bam = write_bam(tmp_path / 'parent.bam', reads, header=header)
    result = scenario.run(bam, dense_support=True, max_records=2, assemble=False,
                          peptide_lengths=(n for n in (8, 9)))
    support = candidate(result)['full_interval_support']
    assert support['fragments'] == 3
    assert support['read_lineage']['signal_groups'] == 0
    assert support['read_lineage']['statuses'] == {'unresolved_parent': 3}
    assert any(r['identity'][1] == 'parent' for r in result['read_lineage']['reads'])


def test_discovery_context_labels_consult_original_mates_outside_search_regions(tmp_path):
    scenario = Scenario()
    reads = aligned('event', scenario.sequence[200:300], scenario.positions[200:300])
    context = record('context', '1', 2090, '60M', scenario.donor_ref.sequence[190:250], flag=65)
    mate = record('context', '2', 18000, '60M', 'A' * 60, flag=129)
    for read, cell in [(context, 'cell-A'), (mate, 'cell-B')]:
        read.set_tag('CB', cell)
        read.set_tag('UB', 'umi')
    reads.extend([context, mate])
    result = scenario.run(write_bam(tmp_path / 'context.bam', reads), dense_support=True,
                          max_records=3, min_alternative_fragments=1, min_overlap=20)
    rows = [row for row in result['cell_umi_evidence']['reads'] if row['identity'][1] == 'context']
    assert len(rows) == 2
    assert {row['status'] for row in rows} == {'conflicting_template_labels'}
    assert all(path['sequence_evidence']['voting_support']['cells'] == 0 for path in result['paths'])


def test_placed_rare_event_precedes_unlinked_breakpoint_clips(tmp_path):
    scenario = Scenario()
    reads = []
    for i in range(40):
        pieces = aligned('unlinked-%03d' % i, scenario.sequence[200:300], scenario.positions[200:300])
        for read in pieces:
            read.set_tag('SA', '1,17000,+,100M,60,0;')
        reads += pieces
    reads += aligned('rare-linked', scenario.sequence[200:300], scenario.positions[200:300])
    result = scenario.run(write_bam(tmp_path / 'clips.bam', reads), *scenario.exact,
                          dense_support=True, max_records=2, assemble=False,
                          read_collector=ReadCollector(use_soft_clipped_bases=True))
    path = candidate(result)
    assert path['sequence'] == scenario.sequence[200:300]
    assert path['event_linkage']['status'] == 'breakpoint_junction'
    assert path['full_interval_support']['fragments'] == 1


def test_dense_exon_clip_priority_and_complete_original_support(tmp_path):
    scenario = Scenario()
    reads = [record('clip-%d' % i, '1', 2120, '30M30S', scenario.sequence[220:280], mapq=60)
             for i in range(16)]
    reads += aligned('donor', scenario.sequence[160:250], scenario.positions[160:250])
    partner = aligned('partner', scenario.sequence[250:350], scenario.positions[250:350])
    partner[0].mapping_quality = 1
    reads += partner + [partner[0]]  # Repeated original SAM payload is not another fragment.
    reads += aligned('ordinary', scenario.donor_ref.sequence[160:350], scenario.donor_positions[160:350])
    reads += [record('background-%d' % i, '1', 2010, '50M', scenario.donor_ref.sequence[10:60])
              for i in range(30)]
    bam = write_bam(tmp_path / 'dense-exon-clip.bam', reads)
    options = dict(read_collector=ReadCollector(min_mapping_quality=0, use_soft_clipped_bases=True),
                   include_regional_candidates=False, exon_anchored_clips=True)

    limited = scenario.run(bam, dense_support=True, max_records=4, **options)
    clip, = [seed for seed in limited['seeds']
             if seed['relation'] == 'exon_anchored_clip_partner_unplaced']
    assert clip['fragments'] == 4
    assert clip['complete_input_support']['fragments'] == 16
    assert len(clip['complete_input_support']['original_read_ids']) == 16
    assert clip['complete_input_support']['raw_mapping_qualities'] == [60]
    assert limited['discovery']['policy'] == 'event_first_full_allele_groups_with_exon_clips.v1'
    assert limited['support_acquisition']['unique_records'] == len(reads) - 1
    assert 'discovery_record_limit' in limited['limitations']

    bounded = scenario.run(bam, max_records=80, **options)
    dense = scenario.run(bam, max_records=80, dense_support=True, **options)
    expected = scenario.sequence[160:350]
    assert {path['sequence'] for path in bounded['paths']} == {expected}
    assert {path['sequence'] for path in dense['paths']} == {expected}
    assert dense['paths'][0]['frame_status'] == 'translated'
    assert dense['paths'][0]['junctions'][0]['direct_support']['fragments'] == 0
    clip, = [seed for seed in dense['seeds']
             if seed['relation'] == 'exon_anchored_clip_partner_unplaced']
    assert clip['complete_input_support']['fragments'] == 16
    assert all(row['identity'][1].startswith('clip-') for row in
               (dense['observations'][key] for key in clip['complete_input_support']['observations']))


def test_dense_exon_clips_replay_original_hcc1395_reads():
    import gzip

    root = Path(__file__).parent / 'data/hcc1395_sv_junction'
    models = json.loads(gzip.open(root / 'ensembl87_models.json.gz', 'rt').read())['models']['SUGT1']
    references = [FusionReference(**row) for row in models]
    bam = root / 'SRR23456709.sv_candidates.grch38.bam'
    options = dict(event_id='chr13_52688329_chr13_60966095_+-', reference_name='GRCh38',
                   sample_id='HCC1395', source='SRA:SRR23456709',
                   donor=FusionBreakpoint('chr13', 52688329, '+'),
                   acceptor=FusionBreakpoint('chr13', 60966095, '+'),
                   regions=[], references=references, event_provenance={'dna_call_id': 'chr13'},
                   read_collector=ReadCollector(min_mapping_quality=20, use_soft_clipped_bases=True),
                   max_records=5000, max_queries=100, max_paths=20,
                   include_regional_candidates=False, exon_anchored_clips=True)
    with pysam.AlignmentFile(str(bam)) as alignment:
        bounded = reconstruct_sv_rna(alignment, **options)
    with pysam.AlignmentFile(str(bam)) as alignment:
        dense = reconstruct_sv_rna(alignment, dense_support=True, **options)
    assert {path['path_id'] for path in dense['paths']} == {path['path_id'] for path in bounded['paths']}
    assert max(j['direct_support']['fragments'] for path in dense['paths'] for j in path['junctions']
               if j['relation'] == 'breakpoint_junction') == 3
    clips = [seed for seed in dense['seeds'] if seed['relation'] == 'exon_anchored_clip_partner_unplaced']
    assert clips and all(seed['complete_input_support']['fragments'] >= 1 for seed in clips if seed['paths'])
