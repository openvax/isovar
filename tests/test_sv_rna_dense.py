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


def test_conflicting_records_and_interrupted_scan_leave_no_partial_result(tmp_path, monkeypatch):
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
    monkeypatch.setattr(sv_rna_dense, '_spool', original)
    conflict = pysam.AlignedSegment.fromstring(reads[-1].to_string(), reads[-1].header)
    conflict.set_tag('CB', 'another-cell')
    bam = write_bam(tmp_path / 'conflicting.bam', reads + [conflict])
    with pytest.raises(ValueError, match='Conflicting original SAM'):
        scenario.run(bam, dense_support=True, scratch_dir=tmp_path)
    assert not list(tmp_path.glob('isovar-sv-*'))


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
