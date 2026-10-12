"""Complete-input support remains separate from bounded hypothesis discovery."""

from dataclasses import asdict
import json

import pysam
import pytest

from isovar import export_sv_rna_orfs
from isovar.cli.isovar_sv_rna import run as cli_run
from isovar.read_collector import ReadCollector
from isovar.fusion import FusionBreakpoint
from isovar.sv_rna import reconstruct_sv_rna
from tests.test_sv_rna import Scenario, aligned, record, write_bam


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
