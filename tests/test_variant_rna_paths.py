"""Original intronic RNA survives missing CDS models and keeps its evidence scope."""

from dataclasses import asdict

import pysam
import pytest

from isovar.fusion import FusionReference
from isovar.rna_path_linkage import OriginalPathMatch, PathFeature, template_path_linkage
from isovar.variant_rna_paths import reconstruct_variant_rna_paths, variant_rna_path_input_from_dict


EXON1 = 'ATG' + 'GCC' * 5
EXON2 = 'AAG' * 6 + 'TAA'


def inputs(*, stop_before=False, contig='chr1'):
    intron = ('TAA' if stop_before else 'GGA') + 'GGA' * 13
    reference = EXON1 + intron + EXON2
    model = FusionReference('coding.1', 'GRCh38', 'synthetic.1', '1', '+', ((100, 118), (160, 181)),
                            EXON1 + EXON2, 0, 39)
    data = dict(reference_name='GRCh38', sample_id='sample',
                variant=dict(contig=contig, start=125, ref='G', alt=''), references=[asdict(model)])
    return data, reference


def bam_file(path, sequence, *, contig='chr1'):
    header = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': contig, 'LN': 1000}], 'RG': [{'ID': 'rg', 'SM': 'sample'}]})
    with pysam.AlignmentFile(path, 'wb', header=header) as bam:
        for name, seq, cigar in [('alt', sequence[:24] + sequence[25:], '24M1D56M'),
                                  ('ref', sequence, '81M'), ('ordinary', EXON1 + EXON2, '18M42N21M')]:
            read = pysam.AlignedSegment(header)
            read.query_name, read.query_sequence = name, seq
            read.query_qualities = [30] * len(seq)
            read.reference_id, read.reference_start, read.mapping_quality = 0, 100, 60
            read.cigarstring = cigar
            read.set_tag('RG', 'rg')
            bam.write(read)
    return path


def run_paths(tmp_path, *, stop_before=False, models=True):
    data, reference = inputs(stop_before=stop_before)
    if not models:
        data['references'] = []
    path = bam_file(tmp_path / 'original.bam', reference)
    with pysam.AlignmentFile(path) as bam:
        return reconstruct_variant_rna_paths(bam, **variant_rna_path_input_from_dict(data), source='original-library',
                                             min_overlap=12, min_anchor_bases=12, flank=100)


def test_no_coding_models_does_not_prevent_intronic_alt_or_ref_assembly(tmp_path):
    result = run_paths(tmp_path, models=False)
    assert result['allele_counts'] == {'ALT': {'fragments': 1, 'segments': 1},
                                      'REF': {'fragments': 1, 'segments': 1},
                                      'OTHER': {'fragments': 0, 'segments': 0}}
    assert result['paths']
    assert all(p['protein_interpretation'] == 'not_assessed' and not p['coding'] for p in result['paths'])
    assert {f['state'] for p in result['paths'] for f in p['focal_alleles']} == {'ALT', 'REF'}
    assert all(p['linkage']['single_read_whole_path_fragments'] >= 1 for p in result['paths'])
    assert all(p['linkage']['paired_template_whole_path_fragments'] == 0 for p in result['paths'])


@pytest.mark.parametrize('stop_before', [False, True])
def test_retention_frame_keeps_stops_before_the_allele_explicit(tmp_path, stop_before):
    result = run_paths(tmp_path, stop_before=stop_before)
    alt_frames = [c for p in result['paths'] for c in p['coding'] if c['focal_state'] == 'ALT']
    assert alt_frames
    upstream = [c for c in alt_frames if c['complete_5prime']]
    assert {c['focal_coverage'] for c in upstream} == ({'after_stop'} if stop_before else {'covered'})
    assert all(c['focal_coverage'] == 'upstream_of_frame_anchor' for c in alt_frames if not c['complete_5prime'])
    for path in result['paths']:
        if path['coding']:
            assert path['linkage']['single_read_whole_path_fragments'] == 1
            assert any(row['links'] for row in path['linkage']['witnesses'])
    assert all(c['protein_status'] == 'conditional_RNA_frame' for c in alt_frames)
    assert all(c['reference_protein'] == 'M' + 'A' * 5 + 'K' * 6 for c in alt_frames)


def test_input_decoder_rejects_typos_and_reference_mismatch(tmp_path):
    data, reference = inputs()
    with pytest.raises(ValueError, match='Unexpected literal RNA'):
        variant_rna_path_input_from_dict(dict(data, refereces=[]))
    data['references'][0]['reference_name'] = 'GRCh37'
    with pysam.AlignmentFile(bam_file(tmp_path / 'original.bam', reference)) as bam:
        with pytest.raises(ValueError, match='reference declarations differ'):
            reconstruct_variant_rna_paths(bam, **variant_rna_path_input_from_dict(data), source='original')


def match(name, segment, low, high, placement='a'):
    return OriginalPathMatch(('rg', name, segment), name + placement, low, high)


def test_alternative_placements_cannot_supply_mate_or_single_read_links():
    features = [PathFeature('allele', 2, 3), PathFeature('exon', 15, 18)]
    alternatives = [match('one', 64, 0, 10), match('one', 64, 10, 20, 'b')]
    result = template_path_linkage(alternatives, 20, features)
    assert result == dict(whole_path_single_read_segments=[], whole_path_paired_template=False, links=[])
    paired = [alternatives[0], match('one', 128, 9, 20)]
    result = template_path_linkage(paired, 20, features)
    assert result['whole_path_paired_template']
    assert result['links'][0]['paired_template']
    assert result['links'][0]['single_read_segments'] == []
    separated = [alternatives[0], match('one', 128, 12, 20)]
    result = template_path_linkage(separated, 20, features)
    assert not result['whole_path_paired_template']
    assert result['links'][0]['paired_template']


def test_allele_links_need_quality_approved_original_and_deletion_flanks():
    left = match('one', 64, 0, 10)
    features = [PathFeature('allele', 5, 5, ('other',)), PathFeature('exon', 1, 4)]
    assert not template_path_linkage([left], 10, features)['links']
    features[0] = PathFeature('allele', 5, 5, (left.observation_key,))
    assert template_path_linkage([left], 10, features)['links'][0]['single_read_segments'] == [('rg', 'one', 64)]
    assert not match('one', 64, 5, 10).covers(features[0])


def test_cli_api_and_cohort_share_one_input_and_reuse_library(tmp_path, monkeypatch):
    from isovar.cli.isovar_variant_rna_paths import run as cli
    from isovar.cohort import run_cohort
    from isovar.sv_rna_audit import file_sha256, read_json, write_json
    import isovar.rna_library as rna_library

    data, reference = inputs()
    bam = bam_file(tmp_path / 'original.bam', reference)
    supplied = tmp_path / 'input.json'
    write_json(supplied, data)
    result_path = tmp_path / 'cli.json.gz'
    parameters = dict(min_overlap=12, min_anchor_bases=12, flank=100)
    cli(['--bam', str(bam), '--input', str(supplied), '--source', 'original', '--output', str(result_path),
         '--min-overlap', '12', '--min-anchor-bases', '12', '--flank', '100'])
    cli_result = read_json(result_path)
    with pysam.AlignmentFile(bam) as reader:
        api = reconstruct_variant_rna_paths(reader, **variant_rna_path_input_from_dict(data), source='original', **parameters)
    # JSON normalizes tuples; no schema-specific translation is needed.
    import json
    assert json.loads(json.dumps(api)) == cli_result
    manifest = dict(schema='isovar.cohort_manifest.v1', samples=[dict(sample_id='sample')],
        files={key: dict(path=str(path), sha256=file_sha256(path)) for key, path in [('bam', bam), ('input', supplied)]},
        sources=[dict(source_id='original', sample_id='sample', library_id='one', bam='bam',
                      scope='Synthetic original RNA', provenance={'method': 'regression'})],
        variant_rna_paths=[dict(id=identifier, source_id='original', sample_id='sample', input='input', parameters=parameters)
                           for identifier in ['event', 'second-view']])
    manifest_path = tmp_path / 'manifest.json'
    write_json(manifest_path, manifest)
    preparations = []
    prepare = rna_library.prepare_rna_library

    def tracked_prepare(*args, **kwargs):
        preparations.append(args[0])
        return prepare(*args, **kwargs)

    monkeypatch.setattr(rna_library, 'prepare_rna_library', tracked_prepare)
    output = tmp_path / 'cohort'
    run_cohort(manifest_path, output, render=False)
    assert len(preparations) == 1
    assert read_json(output / 'variant_rna_paths/event/result.json.gz') == cli_result
    run_cohort(manifest_path, output, render_only=True, render=False)
    run_cohort(manifest_path, output, resume=True, render=False)
    assert len(preparations) == 1
    manifest['variant_rna_paths'][0]['parameters']['min_overlap'] = 13
    write_json(manifest_path, manifest)
    with pytest.raises(ValueError, match='identity differs'):
        run_cohort(manifest_path, output, render_only=True, render=False)


@pytest.mark.parametrize('contig', ['MT', 'chrUn_KI270442v1', 'viral_custom', 'HLA-A*01:01:01:01'])
def test_exact_custom_and_alias_contigs_without_coding_models(tmp_path, contig):
    from isovar.contig_resolution import UnresolvedContigError

    data, reference = inputs(contig=contig)
    data['references'] = []
    bam = bam_file(tmp_path / 'original.bam', reference, contig='chrM' if contig == 'MT' else contig)
    with pysam.AlignmentFile(bam) as reader:
        result = reconstruct_variant_rna_paths(reader, **variant_rna_path_input_from_dict(data), source='original', min_overlap=12)
    assert result['allele_counts']['ALT']['fragments'] == 1
    assert result['resolved_contig'] == ('chrM' if contig == 'MT' else contig)
    data['variant']['contig'] = 'not_in_header'
    with pysam.AlignmentFile(bam) as reader, pytest.raises(UnresolvedContigError):
        reconstruct_variant_rna_paths(reader, **variant_rna_path_input_from_dict(data), source='original')


def test_seed_budget_does_not_truncate_allele_or_support_counts(tmp_path):
    data, reference = inputs()
    path = bam_file(tmp_path / 'initial.bam', reference)
    expanded = tmp_path / 'expanded.bam'
    with pysam.AlignmentFile(path) as original, pysam.AlignmentFile(expanded, 'wb', header=original.header) as bam:
        for read in original:
            bam.write(read)
            if read.query_name == 'alt':
                read.query_name = 'independent-alt'
                bam.write(read)
    with pysam.AlignmentFile(expanded) as bam:
        result = reconstruct_variant_rna_paths(bam, **variant_rna_path_input_from_dict(data), source='original',
                                               max_seed_segments=1, min_overlap=12)
    assert result['allele_counts']['ALT']['fragments'] == 2
    assert result['discovery']['ALT']['seed_limit_reached']
    assert result['discovery']['ALT']['selected_seed_segments'] == 1
    paths = [p for p in result['paths'] if any(f['state'] == 'ALT' for f in p['focal_alleles'])]
    assert paths and all(p['linkage']['single_read_whole_path_fragments'] == 2 for p in paths)


def test_cohort_literal_path_pdf_keeps_conditional_frames_separate(tmp_path):
    from isovar.cohort_pdf import write_cohort_reports
    from isovar.sv_rna_audit import file_sha256, write_json
    from pypdf import PdfReader

    result = run_paths(tmp_path, stop_before=True)
    result_path = tmp_path / 'paths.json.gz'
    write_json(result_path, result)
    source = dict(source_id='original-library', sample_id='sample', library_id='one',
                  scope='Synthetic original RNA', provenance={'method': 'regression'})
    manifest = dict(samples=[dict(sample_id='sample')], sources=[source], files={})
    summary = dict(manifest=str(tmp_path / 'manifest.json'), sources=[source], small_variants=[], additional=[dict(
        family='variant_rna_paths', task=dict(id='retention', sample_id='sample', source_id='original-library'),
        status='recorded', path=str(result_path), sha256=file_sha256(result_path), identity=dict(inputs={}))])
    paths = write_cohort_reports(manifest, summary, tmp_path)
    text = ' '.join(' '.join(page.extract_text() for page in PdfReader(paths[0]).pages).split())
    assert 'after stop' in text
    assert 'paired-only whole-path templates' in text
    assert 'Frames are conditional' in text
    assert 'Prior assessment retained' not in text


def test_ordinary_splice_competitor_is_counted_without_becoming_alt_rna(tmp_path):
    result = run_paths(tmp_path)
    ordinary, = result['reference_paths']
    assert ordinary['origin'] == 'reference_model'
    assert ordinary['reference_transcript_ids'] == ['coding.1']
    assert ordinary['sequence'] == EXON1 + EXON2
    assert len(ordinary['whole_path_segments']) == 1
    assert ordinary['whole_path_segments'][0][1] == 'ordinary'
    assert not ordinary['coding']  # the declared model is not an observed protein
    assert all(p['origin'] == 'observed_assembly' for p in result['paths'])


def test_report_exemplars_include_ref_controls_and_hide_mirror_duplicates(tmp_path):
    from isovar.variant_path_report import _exemplars, _allele_fragments

    data = run_paths(tmp_path)
    selected = _exemplars(data)
    assert len(selected) == 2
    assert {f['state'] for candidate in selected for f in candidate['focal_alleles']} == {'ALT', 'REF'}
    assert all(block[-1] == '+' for candidate in selected for block in candidate['mapping'])
    assert sorted(_allele_fragments(candidate, 'ALT') for candidate in selected) == [0, 1]
