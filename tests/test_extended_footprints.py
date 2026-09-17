"""Real Sid evidence: complete source bookkeeping and no inferred SV protein."""

import pysam

from examples.osteosarc_extended_figures import load_dlg5
from examples.osteosarc_footprint_figures import load as baseline
from isovar.protein_comparison import comparison_rows
from tests.data.osteosarc.figure_comparisons.extended_footprints import load, source_products
from tests.data.osteosarc.figure_comparisons.dlg5 import local_evidence
from tests.data.osteosarc.figure_comparisons.footprints import fusion_footprint, split_deletion_paths
from tests.test_chimeric_phasing import record
from tests.data.fusions.build_osteosarc import query_mapping, extract


def test_all_nine_candidates_and_fourteen_products_are_retained_separately():
    data = load()
    assert len(data["sources"]) == len(source_products()) == 14
    assert len(data["indels"]) + len(data["deletions"]) + len(data["fusions"]) == 9
    sources = {s["id"]: s for s in data["sources"]}
    assert sources["T3-ONT-dedup"]["processing_family"] == sources["T3-ONT-tagged"]["processing_family"]
    assert all(s["receipt"]["assembly"] == "GRCh38" for s in sources.values())
    for old in baseline()["indels"]:
        new = next(e for e in data["indels"] if e["gene"] == old["gene"])
        assert len(new["products"]) == 11  # Tagged predecessors only used for path checks.
        for previous in old["products"]:
            current = next(p for p in new["products"] if p["source"] == previous["source"])
            assert current["default_counts"] == previous["default_counts"]
            for now, before in zip(current["visualization"]["modes"], previous["visualization"]["modes"]):
                assert (now["protein"]["amino_acids"] if now["protein"] else None) == (before["protein"]["amino_acids"] if before["protein"] else None)


def test_every_uncapped_protein_is_validated_and_represented_without_reference_filling():
    for entry in load()["indels"]:
        products = [dict(p, label=p["source"]) for p in entry["products"]]
        rows = comparison_rows(products)
        for product in products:
            for mode in product["visualization"]["modes"]:
                displayed = [r for r in rows if r["source"] == product["source"] and r["assembly"] == mode["assembly"]]
                assert {i for r in displayed for i in r["covered_ranks"]} == set(range(1, len(mode["proteins"]) + 1))
                validation = next(v for v in product["visualization"]["provenance"]["independent_validation"] if v["assembly"] == mode["assembly"])
                assert len(validation["all_results"]) == len(mode["proteins"])
                assert all(c["validation_status"] == "ok" for c in validation["all_results"])
    gtf = next(e for e in load()["indels"] if e["gene"] == "GTF3C5")
    ont = next(p for p in gtf["products"] if p["source"] == "T1-ONT-dedup")
    assert any(not c["matches_expected"] for v in ont["visualization"]["provenance"]["independent_validation"] for c in v["all_results"])
    assert all(m["protein"] for m in ont["visualization"]["modes"])


def test_t3_gabbr1_has_a_path_but_no_quality_passing_window():
    data = load()
    entry = next(e for e in data["fusions"] if e["name"] == "GABBR1--SLC29A1")
    product = next(p for p in entry["products"] if p["source"] == "T3-ONT-tagged")
    header = pysam.AlignmentHeader.from_references(["chr6"], [170805979])
    records = [pysam.AlignedSegment.fromstring(s, header) for p in product["paths"] for s in p["original_records"]]
    result = fusion_footprint(records, (entry["name"], *entry["breakpoints"]))
    assert result["complete_paths"] == product["complete_paths"] == 1
    assert all(not p["windows"] for p in result["paths"])
    primary = next(r for r in records if not r.is_supplementary)
    assert ("chr6", 29612905, "-") not in query_mapping(primary).values()
    assert ("chr6", 29612919, "-") in query_mapping(primary).values()
    assert min(primary.get_forward_qualities()[1013:1068]) == 3
    assert extract(primary, "chr6", 29612920, "chr6", 44218908, 20, records, ("-", "+")) is None
    # Counterfactual diagnostic only: distinguish endpoint failure from QUAL.
    # Originals remain pinned; these in-memory records are disposable copies.
    for r in records:
        r.query_qualities = [60] * len(r.query_sequence)
    assert extract(primary, "chr6", 29612906, "chr6", 44218908, 20, records, ("-", "+")) is None
    assert extract(primary, "chr6", 29612920, "chr6", 44218908, 20, records, ("-", "+")) is not None


def test_dlg5_caller_fields_are_not_silently_repaired_into_assembled_haplotype():
    data = load_dlg5()
    assert data["purple_signature_in_all_dragen_contigs"]
    assert all(not c["signature_in_contig"] for c in data["dragen"])
    assert all(c["insertion"] == "GAAATGATGC" for c in data["dragen"])
    assert data["signatures"]["dragen_fields"] != data["signatures"]["purple_bnd"]
    assert all(p["signatures"]["dragen_fields"]["sequence_only"] == 0 for p in data["products"])


def test_dlg5_start_and_exact_junction_counts_recount_from_original_sam():
    data = load_dlg5()
    for product in data["products"]:
        header = pysam.AlignmentHeader.from_text(product["header"])
        records = [pysam.AlignedSegment.fromstring(s, header) for s in product["original_records"]]
        result = local_evidence(records, data["signatures"])
        assert result["signatures"] == product["signatures"]
        assert result["start_codon"] == product["start_codon"]
        originals = {s for group in product["split_paths"]["observed_cigar"] for s in group}
        split_records = [pysam.AlignedSegment.fromstring(s, header) for s in originals]
        assert len(split_deletion_paths(split_records, 77850921, 77930460)) == len(product["split_paths"]["observed_cigar"])
    by_id = {p["source"]["id"]: p for p in data["products"]}
    assert {s: by_id[s]["signatures"]["purple_bnd"]["q20"] for s in ("T0-DNA", "T1-DNA", "T2-DNA", "T1-organoid-DNA")} == {
        "T0-DNA": 0, "T1-DNA": 2, "T2-DNA": 4, "T1-organoid-DNA": 4}
    assert by_id["T3-scRNA"]["start_codon"]["reference_start"] == 24
    assert by_id["T3-CD45neg"]["start_codon"]["reference_start"] == 22
    assert all(p["signatures"]["purple_bnd"]["q20"] == 0 for p in data["products"] if p["source"]["product"] != "DNA")


def test_low_quality_mate_does_not_cancel_an_observed_start_or_invent_a_conflict():
    good, poor = record(cigar="3M", flag=83), record(cigar="3M", flag=163)
    for r in (good, poor):
        r.query_sequence = "CAT"
    good.query_qualities, poor.query_qualities = [30] * 3, [5] * 3
    assert local_evidence([good, poor], {}, (100, 103))["start_codon"] == {"reference_start": 1}
    poor.query_sequence, poor.query_qualities = "CAA", [30] * 3
    assert local_evidence([good, poor], {}, (100, 103))["start_codon"] == {"conflicting_template": 1}
    skip = record(cigar="1M3N3M")
    assert local_evidence([skip], {}, (101, 104))["start_codon"] == {}


def test_missing_pacbio_qualities_remain_an_input_limitation():
    for entry in load()["indels"]:
        p = next(p for p in entry["products"] if p["source"] == "T1-PacBio")
        quality = p["visualization"]["provenance"]["input_quality"]
        assert quality["primary_records_missing_qualities"] > 0
        assert p["default_counts"]["alt"] == 0
