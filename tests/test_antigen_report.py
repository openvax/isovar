"""Mixed-source evidence remains separate and viral support uses fragments."""

import csv
import json

import pytest

from isovar.antigen_report import compose
from reports.hcc1395_colo829.build_sv_exon_annotations import sv_footprints


def write_json(path, payload):
    path.write_text(json.dumps(payload))
    return path


def small_variant_tsv(path):
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
    content = "\n".join(page.extract_text() for page in pypdf.PdfReader(output).pages)
    assert "Exon-involved small variants" in content
    assert "Viral contigs by aligned fragments" in content
    assert "Viral contigs by sequence-assigned fragments" in content
    assert "Structural variants and supplied fusions" in content
    assert "Targeted antigen-gene isoforms" in content


def test_missing_categories_and_wrong_sample_or_schema_are_explicit(tmp_path):
    small = small_variant_tsv(tmp_path / "small.tsv")
    inventory = compose("tumor", small_variants=small)
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
