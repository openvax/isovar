"""Original matched DNA/RNA evidence; no assumed origin for the RNA indel."""

import pytest

from tests.data.osteosarc.figure_comparisons import nr2f2


@pytest.mark.parametrize("sample", ["normal", "tumor", "rna"])
@pytest.mark.parametrize("mode,options", [("q20", {}), ("q30", {"min_quality": 30}), ("mapq1", {"min_mapq": 1})])
def test_original_records_reproduce_anchored_counts(sample, mode, options):
    data = nr2f2.load()
    source = data["sources"][sample]
    assert nr2f2.recount(source, data["reference"], **options) == source["audits"][mode]


def test_rna_deletion_is_not_confirmed_by_matched_dna_or_independent_endpoints():
    sources = nr2f2.load()["sources"]
    for key, expected in [("normal", {"reference": 137, "other": 2}),
                          ("tumor", {"reference": 323, "other": 3}),
                          ("rna", {"reference": 23, "deletion": 5})]:
        assert sources[key]["audits"]["q20"]["counts"]["deletion"] == expected
    rna = sources["rna"]["audits"]["q20"]
    assert rna["deletion_alignment_paths"] == rna["deletion_fragment_endpoints"] == 1
    assert rna["counts"]["joint"] == {"alternate/deletion": 4, "reference/reference": 22}
    assert sources["tumor"]["audits"]["q20"]["counts"]["joint"]["alternate/reference"] == 16
    # Stricter qualities reduce support; they never create DNA confirmation.
    assert sources["rna"]["audits"]["q30"]["counts"]["deletion"]["deletion"] == 3
    assert all("deletion" not in sources[s]["audits"]["q30"]["counts"]["deletion"] for s in ["normal", "tumor"])
