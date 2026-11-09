"""Mixed-source evidence remains separate and viral support uses fragments."""

import csv
import json

import pytest

from isovar.antigen_report import compose
from reports.hcc1395_colo829.build_sv_exon_annotations import sv_footprints


def write_json(path, payload):
    path.write_text(json.dumps(payload))
    return path


def small_variant_tsv(path, *, ncorf_assessed=True):
    rows = [
        dict(sample_id="tumor", call_id="SNV", ref="A", alt="G",
             exon_overlap="True", alt_reads="2", protein_count="1",
             orfdata_coding_block_overlap="False"),
        dict(sample_id="tumor", call_id="INS", ref="A", alt="AT",
             exon_overlap="True", alt_reads="0", protein_count="0",
             orfdata_coding_block_overlap="True"),
        dict(sample_id="tumor", call_id="OUT", ref="C", alt="T",
             exon_overlap="False", alt_reads="0", protein_count="0",
             orfdata_coding_block_overlap="False"),
    ]
    if not ncorf_assessed:
        for row in rows:
            row["orfdata_coding_block_overlap"] = ""
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def test_sv_exon_footprint_uses_span_flanks_and_unknown_mate():
    base = dict(chrom1="chr1", pos1="100", chrom2="chr1", pos2="200")
    assert sv_footprints(dict(base, **{"class": "DEL"})) == (
        [("1", 100, 200)], "intrachromosomal_span")
    assert sv_footprints(dict(base, **{"class": "INS"})) == (
        [("1", 100, 100), ("1", 101, 101), ("1", 200, 200)],
        "insertion_flanks")
    assert sv_footprints(dict(base, **{"class": "BND"}, pos2="NA")) == (
        [("1", 100, 100)], "one_nominated_breakpoint_only")
    assert sv_footprints(dict(base, **{"class": "TRA"}, chrom2="chr2")) == (
        [("1", 100, 100), ("2", 200, 200)], "nominated_breakpoints")


def test_mixed_source_report_preserves_denominators_and_viral_zero_bins(tmp_path):
    pytest.importorskip("reportlab")
    pypdf = pytest.importorskip("pypdf")
    from isovar.cli.isovar_antigen_report import run
    from isovar.evidence_pdf import viral_support_figures

    small = small_variant_tsv(tmp_path / "small.tsv")
    sv = write_json(tmp_path / "sv.json", dict(
        schema="isovar.sv_rna_candidates.v5", sample_id="tumor", source="dna-sv",
        event_id="sv1", status="event_linked_candidates", paths=[dict(
            event_linkage=dict(status="breakpoint_junction"),
            junctions=[dict(relation="breakpoint_junction",
                            direct_support=dict(fragments=3))],
            translations=[dict(amino_acids="MUT",
                               departure_relations=["breakpoint_junction"])],
            exploratory_orfs=dict(candidates=[dict(amino_acids="ALT")]))]))
    fusion = write_json(tmp_path / "fusion.json", dict(
        schema="isovar.fusion_rna.v3", sample_id="tumor", event_id="fusion1",
        provenance=dict(source="rna-fusion"), status="translated",
        paths=[dict(junctions=[dict(direct_support=dict(fragments=2))],
                    translations=[dict(amino_acids="FUS")])]))
    contigs = [dict(contig="HPV16", aligned_fragments=10,
                    sequence_assigned_fragments=8, status="sequences_reconstructed"),
               dict(contig="EBV", aligned_fragments=2,
                    sequence_assigned_fragments=0, status="aligned_but_not_reconstructed"),
               dict(contig="HBV", aligned_fragments=0,
                    sequence_assigned_fragments=0, status="no_eligible_alignments")]
    viral = write_json(tmp_path / "viral.json", dict(
        schema="isovar.viral_sequences.v1", sample_id="tumor", source="rna-viral",
        contigs=contigs, sequences=[dict(id="path1", coding=[dict(
            amino_acids="VIRAL", status="complete_CDS")])],
        discovery=dict(truncated=False)))
    isoforms = write_json(tmp_path / "isoforms.json", dict(
        schema="isovar.antigen_isoforms.v1", sample_id="tumor", source="rna-cta",
        targets=[dict(gene_id="G1"), dict(gene_id="G2")],
        isoforms=[dict(id="observed", gene_ids=["G1"], compatible_fragments=5)],
        partial_sequences=[dict(id="reference-only", gene_ids=["G2"],
                                compatible_fragments=0)],
        protein_groups=[dict(fully_covered_isoform_ids=["observed"],
                             compatible_fragments=5)],
        discovery=dict(truncated=False)))
    output = tmp_path / "all-categories.pdf"
    run(["--sample-id", "tumor", "--output", str(output),
         "--small-variant-tsv", str(small), "--sv-result", str(sv),
         "--fusion-result", str(fusion), "--viral-result", str(viral),
         "--isoforms-result", str(isoforms)])
    summary = json.loads(output.with_suffix(".json").read_text())
    categories = summary["categories"]
    assert categories["small_variants"]["records"][0]["exon_involved"] == 2
    assert categories["small_variants"]["records"][0]["exonic_reads_by_class"] == {
        "SNV": [2], "insertion": [0]}
    assert categories["sv"]["records"][0]["event_linked_translated_paths"] == 1
    assert categories["sv"]["records"][0]["max_direct_junction_fragments"] == 3
    assert categories["fusion"]["records"][0]["translated_candidates"] == 1
    assert categories["viral"]["records"][0]["assigned_contigs"] == 1
    assert categories["isoforms"]["records"][0]["rna_supported_genes"] == 1
    figures = viral_support_figures(contigs)
    assert figures[0].values["0"] == 1
    assert "0" not in figures[1].bins
    assert figures[2].values["0"] == 2
    content = " ".join(" ".join(page.extract_text() for page in pypdf.PdfReader(output).pages).split())
    assert "Small variants: named calls" in content
    assert "HPV16" in content and "EBV" in content and "HBV" in content
    assert "Sequence-assigned fragments" in content
    assert "Targets: G1, G2" in content
    assert "SNV" in content and "INS" in content and "OUT" in content
    assert "Structural variants and supplied fusions" in content
    assert "Targeted antigen-gene isoforms" in content


def test_missing_categories_and_wrong_sample_or_schema_are_explicit(tmp_path):
    small = small_variant_tsv(tmp_path / "small.tsv")
    inventory = compose("tumor", small_variants=small)
    assert inventory["categories"]["small_variants"]["records"][0][
        "ncorf_reference_candidates"] == 1
    assert inventory["categories"]["viral"]["status"] == "not_assessed"
    assert inventory["categories"]["sv"]["status"] == "not_assessed"
    wrong = write_json(tmp_path / "wrong.json", dict(
        schema="isovar.viral_sequences.v1", sample_id="other", source="rna",
        contigs=[], sequences=[]))
    with pytest.raises(ValueError, match="sample_id"):
        compose("tumor", viral=[wrong])
    wrong.write_text(json.dumps({"schema": "unsupported", "sample_id": "tumor"}))
    with pytest.raises(ValueError, match="expected schema"):
        compose("tumor", viral=[wrong])
    with pytest.raises(ValueError, match="at least one"):
        compose("tumor")

    small_variant_tsv(small, ncorf_assessed=False)
    inventory = compose("tumor", small_variants=small)
    assert inventory["categories"]["small_variants"]["records"][0][
        "ncorf_reference_candidates"] is None
    pytest.importorskip("reportlab")
    pypdf = pytest.importorskip("pypdf")
    from isovar.antigen_evidence_pdf import write_antigen_report

    output = tmp_path / "unassessed-ncorf.pdf"
    write_antigen_report(output, inventory)
    content = " ".join(" ".join(page.extract_text() for page in pypdf.PdfReader(output).pages).split())
    assert "Reference ncORF coding-block overlap: not assessed" in content


def test_sv_translation_must_cross_supported_event_junction():
    from isovar.antigen_report import sv_summary

    result = dict(event_id="event", status="event_linked_candidates", paths=[dict(
        event_linkage=dict(status="breakpoint_junction"),
        junctions=[dict(relation="breakpoint_junction",
                        direct_support=dict(fragments=2)),
                   dict(relation="regional_novel_junction",
                        direct_support=dict(fragments=9))],
        translations=[dict(amino_acids="REFERENCE",
                           departure_relations=["regional_novel_junction"])],
        exploratory_orfs=dict(candidates=[]))])
    summary = sv_summary(result)
    assert summary["max_direct_junction_fragments"] == 2
    assert summary["rna_linked_events"]
    assert summary["event_linked_translated_paths"] == 0


def test_mixed_sample_cli_preserves_separate_original_sources(tmp_path):
    from isovar.cli.isovar_antigen_report import run
    from pypdf import PdfReader
    paths = []
    for sample in ('T0', 'T1'):
        path = small_variant_tsv(tmp_path / (sample + '.tsv'))
        with path.open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        for row in rows:
            row['sample_id'] = sample
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        paths.append(path)
    output = tmp_path / 'Sid.pdf'
    args = ['--sample-id', 'Sid', '--output', str(output)]
    for sample, path in zip(('T0', 'T1'), paths):
        args += ['--small-variant-tsv', str(path), '--input-sample-id', str(path), sample]
    run(args)
    summary = json.loads(output.with_suffix('.json').read_text())
    records = summary['categories']['small_variants']['records']
    assert summary['sample_id'] == 'Sid'
    assert [row['input_sample_id'] for row in records] == ['T0', 'T1']
    assert [row['input_records'] for row in records] == [3, 3]
    assert [row['calls'][0]['sample_id'] for row in records] == ['T0', 'T1']
    assert [row['calls'][0]['alt_reads'] for row in records] == ['2', '2']
    pages = PdfReader(output).pages
    text = '\n'.join(page.extract_text() for page in pages)
    assert 'T0' in text and 'T1' in text and 'Sid' in text
    assert len(pages) <= 2
    first = output.read_bytes()
    run(args)
    assert output.read_bytes() == first


def test_input_mapping_rejects_unknown_paths_and_preserves_sample_checks(tmp_path):
    path = small_variant_tsv(tmp_path / 'input.tsv')
    with pytest.raises(ValueError, match='unsupplied'):
        compose('Sid', small_variants=path, input_sample_ids={tmp_path / 'other.tsv': 'tumor'})
    with pytest.raises(ValueError, match='sample_id'):
        compose('Sid', small_variants=path, input_sample_ids={path: 'wrong'})
    with pytest.raises(ValueError, match='nonempty'):
        compose('Sid', small_variants=path, input_sample_ids={path: ''})
    with pytest.raises(ValueError, match='Duplicate input export'):
        compose('tumor', small_variants=[path, path])


def test_isoform_rows_use_producer_counts_and_distinct_original_span_witnesses(tmp_path):
    from isovar.antigen_report import isoform_summary
    first, second = ['rg', 'fragment', 64], ['rg', 'fragment', 128]
    result = dict(targets=[dict(gene_id='G', gene_name='CTA1')], source='rna',
        isoforms=[dict(id='isoform-1', gene_ids=['G'], rna_support=dict(reads=7, fragments=5),
            supporting_segments=7, supporting_fragments=5, compatible_fragments=4, unique_fragments=2,
            whole_path_segments=[first, first, second], sequence_status='whole_path_observed',
            coding=[dict(status='complete_CDS', whole_cds_segments=[first, second]),
                    dict(status='partial_CDS', whole_cds_segments=[])])],
        partial_sequences=[], protein_groups=[], support=dict(scope='complete_declared_input'))
    summary = isoform_summary(result)
    rows = {row['subject']: row for row in summary['measurements']}
    assert rows['RNA path contributors']['reads'] == 7
    assert rows['RNA path contributors']['fragments'] == 5
    assert rows['Compatible with RNA path']['fragments'] == 4
    assert rows['Unique among retained paths']['fragments'] == 2
    assert rows['Whole retained RNA path']['reads'] == 2
    assert rows['Whole retained RNA path']['fragments'] == 1
    assert rows['Whole complete CDS 1']['reads'] == 2
    assert rows['Whole complete CDS 2']['reads'] is None
    assert all('CTA1' in row['name'] for row in rows.values())


def test_fusion_reason_and_missing_read_count_remain_explicit(tmp_path):
    from isovar.antigen_report import fusion_summary
    result = dict(event_id='fusion', status='insufficient_support',
        reasons=['insufficient_direct_junction_fragments'],
        paths=[dict(junctions=[dict(direct_support=dict(fragments=0))])])
    summary = fusion_summary(result)
    assert summary['reasons'] == result['reasons']
    assert summary['measurements'][0]['fragments'] == 0
    assert summary['measurements'][0]['reads'] is None


def test_long_named_tables_paginate_without_dropping_rows(tmp_path):
    from pypdf import PdfReader
    from isovar.antigen_evidence_pdf import write_antigen_report

    path = small_variant_tsv(tmp_path / 'small.tsv')
    inventory = compose('tumor', small_variants=path)
    record = inventory['categories']['small_variants']['records'][0]
    template = record['calls'][0]
    record['calls'] = [dict(template, call_id=f'call-{index:03d}', gene=f'GENE{index:03d}',
                          top_effect='p.A1G') for index in range(60)]
    record.update(input_records=60, exon_involved=60, rna_alt_observed=60,
                  exonic_reads_by_class={'SNV': [2] * 60})
    output = tmp_path / 'long-table.pdf'
    write_antigen_report(output, inventory)
    pages = PdfReader(output).pages
    content = ' '.join(page.extract_text() for page in pages)
    assert len(pages) > 1
    assert all(f'GENE{index:03d}' in content for index in range(60))
    assert content.count('ALT reads') >= 2  # Repeated headers on split tables.


def test_native_isoform_and_viral_exports_keep_original_support_units(tmp_path):
    from pypdf import PdfReader
    from isovar.antigen_evidence_pdf import write_antigen_report
    from .test_antigen_isoforms import target_models, aligned, reconstruct as isoforms
    from .test_viral_sequences import inputs, read, reconstruct as viral

    targets, positions = target_models()
    cta_dir, viral_dir = tmp_path / 'cta', tmp_path / 'viral'
    cta_dir.mkdir()
    viral_dir.mkdir()
    cta = isoforms(cta_dir, [r for i in range(2)
        for r in aligned(f'full-{i}', targets[0].reference.sequence, positions[0])], targets[:1])
    viral_inputs = inputs.__wrapped__(viral_dir)
    virus = viral(viral_dir, viral_inputs,
                  [read(viral_inputs, f'virus-{i}', 'HPV16', 100) for i in range(3)])
    cta_path = write_json(cta_dir / 'result.json', cta)
    viral_path = write_json(viral_dir / 'result.json', virus)
    inventory = compose('report-group', isoforms=[cta_path], viral=[viral_path],
                        input_sample_ids={cta_path: 'sample', viral_path: 'tumor'})
    cta_rows = inventory['categories']['isoforms']['records'][0]['measurements']
    viral_rows = inventory['categories']['viral']['records'][0]['measurements']
    assert {row['reads'] for row in cta_rows if row['subject'] in
            ('RNA path contributors', 'Whole retained RNA path', 'Whole complete CDS 1')} == {2}
    assert {row['reads'] for row in viral_rows if row['subject'] in
            ('RNA path contributors', 'Whole retained RNA path')} == {3}
    assert next(row for row in cta_rows if row['subject'] == 'Compatible with RNA path')['fragments'] == 2
    output = tmp_path / 'native-paths.pdf'
    write_antigen_report(output, inventory)
    content = ' '.join(page.extract_text() for page in PdfReader(output).pages)
    assert 'HPV16' in content and 'CTA' in content and 'Whole complete CDS 1' in content
