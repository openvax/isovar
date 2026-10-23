"""Original-read regression for selected HCC1395 and COLO829 SV evidence."""

import csv

from reports.hcc1395_colo829.build_sv_junction_audit import write_audit


def test_low_mapq_junctions_are_retained_and_assembly_support_is_distinct(tmp_path):
    output = tmp_path / "sv_junction_audit.tsv"
    write_audit(output)
    with output.open() as stream:
        rows = {row["call_id"]: row for row in csv.DictReader(stream, delimiter="\t")}

    chr9 = rows["chr9_83969153_chr9_84015113_-+"]
    assert (chr9["screen_split_candidates"], chr9["junction_fragments_all_mapq"]) == ("7", "7")
    assert chr9["assembled_event_paths_retained"] == "3"
    assert len(chr9["junction_read_ids_all_mapq"].split(";")) == 7
    assert "path_limit" in chr9["assembly_limitations"]

    chr13 = rows["chr13_52688329_chr13_60966095_+-"]
    assert (chr13["screen_split_candidates"], chr13["junction_fragments_all_mapq"],
            chr13["junction_fragments_mapq20_subset"]) == ("5", "4", "3")
    chr21 = rows["chr21_42975947_chr21_43035721_-+"]
    assert (chr21["junction_fragments_all_mapq"], chr21["assembled_event_paths_retained"]) == ("2", "2")

    # The six mate links are not six junctions. Actual split alignments from
    # three distinct fragments survive when low-MAPQ pieces are retained.
    sv3 = rows["COLO829_SOMATIC_SV3_DEL"]
    assert (sv3["screen_mate_candidates"], sv3["junction_fragments_all_mapq"],
            sv3["junction_fragments_mapq20_subset"]) == ("6", "3", "0")
    assert sv3["junction_piece_mapqs_all_mapq"] == "1;3;6;52"
    assert sv3["junction_unplaced_bases"] == "CTATAGTAGAACTTCTTTT"
    assert int(sv3["max_assembled_context_bases"]) > 1000
    assert "extension_segment_limit" in sv3["assembly_limitations"]
    for call_id in ("COLO829_SOMATIC_SV42_DEL", "COLO829_SOMATIC_SV49_DEL"):
        assert rows[call_id]["screen_mate_candidates"] == "1"
        assert rows[call_id]["junction_fragments_all_mapq"] == "0"
