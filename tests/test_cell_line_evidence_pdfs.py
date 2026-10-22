"""Check the published cell-line report counts and reference comparisons."""

from pathlib import Path

import pytest


@pytest.mark.parametrize(
    ("sample", "sv_counts", "example_calls"),
    [
        ("HCC1395", (1788, 354, 172, 6, 0, 0), ("COSM32759", "COSM32783")),
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
    assert tuple(row["call_id"] for row in data["example_references"]) == example_calls
    for row in data["example_references"]:
        mutant = data["summary"]["top"][row["call_id"]]
        start = int(mutant["mutation_start_idx"])
        assert row["reference_context"][:start] == mutant["protein_sequence"][:start]
        assert row["reference_context"] != mutant["protein_sequence"]

    build_one(sample, tmp_path)
    pdf = pypdf.PdfReader(Path(tmp_path) / f"{sample}_somatic_RNA_evidence.pdf")
    cover = pdf.pages[0].extract_text()
    assert "Small variants: DNA to RNA and local protein" in cover
    assert "Structural variants: DNA to RNA breakpoint evidence" in cover
    assert "Split-read or cross-end mate candidate" in cover
    assert "RNA junction sequence resolved" in cover
    full_text = "\n".join(page.extract_text() for page in pdf.pages)
    assert "RNA link candidates per exon-involved SV" in full_text
    assert "Reference and mutant protein windows" in full_text
    assert all(call_id in full_text for call_id in example_calls)
