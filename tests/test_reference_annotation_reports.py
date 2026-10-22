"""Offline consistency checks for the HCC1395/COLO829 reference audit."""

import csv
import gzip
from collections import Counter
from pathlib import Path


REPORT = Path(__file__).parents[1] / "reports/hcc1395_colo829"
ANNOTATIONS = REPORT / "reference_coding_orfdata_annotations.tsv.gz"


def rows(path):
    with gzip.open(path, "rt", newline="") as source:
        return list(csv.DictReader(source, delimiter="\t"))


def test_reference_annotations_match_each_source_atlas():
    annotations = rows(ANNOTATIONS)
    assert len(annotations) == 87136
    assert not any(row["annotation_error"] for row in annotations)
    cases = (
        ("HCC1395", "hcc1395_direct_rna", 41072, 1728, 292, 491, 17, 505, 135, 5),
        ("COLO829", "colo829", 46064, 1395, 502, 373, 5, 378, 180, 2),
    )
    for (sample, prefix, expected, exonic, exonic_alt, coding, ncorf,
         union, coding_alt, ncorf_alt) in cases:
        calls = rows(REPORT / f"{prefix}_atlas_calls.tsv.gz")
        per_call = {r["call_id"]: r for r in annotations if r["sample"] == sample}
        assert len(calls) == len(per_call) == expected
        assert {r["call_id"] for r in calls} == set(per_call)
        counts = Counter()
        for call in calls:
            annotation = per_call[call["call_id"]]
            assert (annotation["chrom"], annotation["pos"], annotation["ref"],
                    annotation["alt"]) == (call["chrom"], call["pos"], call["ref"],
                                            call["alt"])
            is_coding = annotation["reference_coding_effect"] == "True"
            is_ncorf = annotation["orfdata_coding_block_overlap"] == "True"
            is_exonic = annotation["exon_overlap"] == "True"
            counts["exonic"] += is_exonic
            counts["coding"] += is_coding
            counts["ncorf"] += is_ncorf
            counts["union"] += is_coding or is_ncorf
            if call["rna_status"] == "alternate_observed":
                counts["exonic_alt"] += is_exonic
                counts["coding_alt"] += is_coding
                counts["ncorf_alt"] += is_ncorf
        assert tuple(counts[key] for key in (
            "exonic", "exonic_alt", "coding", "ncorf", "union",
            "coding_alt", "ncorf_alt")) == (
                exonic, exonic_alt, coding, ncorf, union, coding_alt, ncorf_alt)


def test_known_orfdata_hits_are_reference_candidates():
    by_key = {(r["sample"], r["call_id"]): r for r in rows(ANNOTATIONS)}
    hcc = by_key["HCC1395", "chr1:11805214:C>G"]
    assert hcc["orfdata_orf_ids"] == "c1norep35"
    assert hcc["orfdata_overlap_rule"] == "changed_reference_bases"
    colo = by_key["COLO829", "chr17:58089866:G>A"]
    assert colo["orfdata_orf_ids"] == "c17norep133;c17riboseqorf123"
    assert colo["orfdata_orf_count"] == "2"
    crossing = by_key["HCC1395", "chr14:20256549:CCTGAATGGGGATCATTCTACACAATACTCCT>C"]
    assert crossing["genomic_transcript_overlap"] == "False"
    assert crossing["exon_overlap"] == "True"
