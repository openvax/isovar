"""Original-read regression for selected HCC1395 and COLO829 SV evidence."""

import csv

from reports.hcc1395_colo829.build_sv_junction_audit import write_audit


def test_candidate_screen_is_separate_from_validated_junctions_and_paths(tmp_path):
    output = tmp_path / "sv_junction_audit.tsv"
    write_audit(output)
    with output.open() as stream:
        rows = {row["call_id"]: row for row in csv.DictReader(stream, delimiter="\t")}

    chr9 = rows["chr9_83969153_chr9_84015113_-+"]
    assert (chr9["screen_split_candidates"], chr9["mapq20_junction_fragments"],
            chr9["mapq20_max_full_path_fragments"]) == ("7", "7", "1")
    assert chr9["mapq20_event_paths_retained"] == "2"
    assert len(chr9["mapq20_junction_read_ids"].split(";")) == 7
    assert "path_limit" in chr9["mapq20_limitations"]

    chr13 = rows["chr13_52688329_chr13_60966095_+-"]
    assert (chr13["screen_split_candidates"], chr13["mapq1_junction_fragments"],
            chr13["mapq20_junction_fragments"], chr13["mapq20_max_full_path_fragments"]) == (
                "5", "4", "3", "1")
    chr21 = rows["chr21_42975947_chr21_43035721_-+"]
    assert (chr21["mapq20_junction_fragments"], chr21["mapq20_event_paths_retained"],
            chr21["mapq20_max_full_path_fragments"]) == ("2", "2", "1")

    # COLO829's six mate links at SV3 are not six junctions. Three additional
    # low-MAPQ split records appear only with the permissive collection gate.
    sv3 = rows["COLO829_SOMATIC_SV3_DEL"]
    assert (sv3["screen_mate_candidates"], sv3["mapq1_junction_fragments"],
            sv3["mapq20_junction_fragments"], sv3["mapq20_max_full_path_fragments"]) == (
                "6", "3", "0", "0")
    for call_id in ("COLO829_SOMATIC_SV42_DEL", "COLO829_SOMATIC_SV49_DEL"):
        assert rows[call_id]["screen_mate_candidates"] == "1"
        assert rows[call_id]["mapq20_junction_fragments"] == "0"
