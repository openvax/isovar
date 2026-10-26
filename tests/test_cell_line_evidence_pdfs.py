"""Check the published cell-line report counts and reference comparisons."""

from pathlib import Path

import pytest


@pytest.mark.parametrize(
    ("sample", "sv_counts", "example_calls"),
    [
        ("HCC1395", (1788, 354, 172, 6, 3, 0), ("COSM32759", "COSM32783")),
        ("COLO829", (188, 97, 64, 6, 0, 0),
         ("chr9:21971154:CCG>C", "chr2:105361307:G>GCCCC",
          "chr1:156265990:G>A")),
    ],
)
def test_cell_line_sv_funnels_and_protein_comparisons(tmp_path, sample, sv_counts,
                                                      example_calls):
    pypdf = pytest.importorskip("pypdf")
    from reports.hcc1395_colo829.build_evidence_pdfs import build_one, collect

    data = collect(sample)
    assert tuple(data["counts"][name] for name in (
        "sv_total", "sv_exonic", "sv_exonic_covered", "sv_exonic_candidate",
        "sv_exonic_junction", "sv_exonic_product")) == sv_counts
    assert set(example_calls) <= {row['call_id'] for row in data['example_references']}
    for row in data["example_references"]:
        mutant = data["summary"]["top"][row["call_id"]]
        start = int(mutant["mutation_start_idx"])
        if row['placement_method'] == 'unique_shared_protein_prefix':
            assert row["reference_context"][:start] == mutant["protein_sequence"][:start]
        else:
            assert row['placement_method'] == 'genomic_codon_anchor_RNA_flanks_differ'
        if row['call_id'] in example_calls:
            assert row["reference_context"] != mutant["protein_sequence"]

    build_one(sample, tmp_path)
    pdf = pypdf.PdfReader(Path(tmp_path) / f"{sample}_somatic_RNA_evidence.pdf")
    cover = pdf.pages[0].extract_text()
    assert "Small variants: DNA to RNA and local protein" in cover
    assert "Structural variants: complementary routes" in cover
    assert "Soft clips / unmapped" in cover
    assert "not an eligibility gate" in cover
    full_text = "\n".join(page.extract_text() for page in pdf.pages)
    assert "RNA link candidates per exon-involved SV" in full_text
    assert "Original-read replay: junction versus RNA path" in full_text
    assert "one read need not span the assembled context" in full_text
    if sample == "COLO829":
        assert "1,213" in full_text
    assert 'Variant landscape' in full_text
    assert 'Highest-support SNVs' in full_text
    assert 'Top RNA-supported frameshift nominations' in full_text
    assert 'ncORFs and other antigen-source categories' in full_text
    assert 'IGV (GRCh38)' in full_text
    assert 'Focal change' in full_text
    if sample == 'HCC1395':
        assert 'Unclassified' in full_text or 'unclassified provenance' in full_text
        assert '531-aa' in full_text and '149 exact Illumina junction' in full_text
    else:
        assert '5,248 fragments support' in full_text and 'CNIH3' in full_text
