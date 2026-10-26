"""Source-derived SV3 exon links remain distinct from an event-linked frame."""

from reports.hcc1395_colo829.build_sv3_frame_anchor_audit import (
    _aligned_exon_clip_contains, build)


def test_source_exon_clips_do_not_phase_the_sv3_junction():
    result = build()
    assert result["schema"] == "isovar.sv3_frame_anchor_audit.v1"
    assert result["bam_alignment_records_scanned"] == 52212
    assert result["bam_spots_scanned"] == 26976
    assert result["bam_spots_with_coding_exon_alignment"] == 673
    assert result["exon_and_path_spot_count"] == 9
    assert result["exon_and_path_spots_with_30_base_match"] == 6
    assert result["exon_and_path_spots_with_cds_alignment"] == 0
    assert result["exact_sv3_join_spots_with_exon_alignment"] == []
    assert result["raw_exact_junction_spots_without_exact_bam_window"] == [
        "SRR8615617.110737497"]
    assert all(link["exact_sv3_43mer"] is False for link in result["exon_path_links"])
    assert all(match["path_region"] == "downstream_reference"
               and match["distance_after_join_bases"] == 280
               and match["same_end_exon_clip_contains_match"]
               for link in result["exon_path_links"] for match in link["matches"])
    strong = next(link for link in result["exon_path_links"]
                  if link["spot"] == "SRR8615617.7306533")
    assert strong["matches"][0]["genomic_interval_0based_half_open"] == [
        "chr1", 224598611, 224598642]
    assert any(exon["exon_end_0based"] == 224587420
               for alignment in strong["exon_alignments"] for exon in alignment["exons"])
    assert result["coding_reference_starts"] == [
        {"transcript_id": "ENST00000272133", "start_grch38_1based": 224617175}]
    assert [row["same_end_start_to_core_join_spots"] for row in result["exploratory_orfs"]] == [1, 2]
    assert [row["same_end_start_to_junction_spots"] for row in result["exploratory_orfs"]] == [0, 0]
    assert result["conclusion"] == "no_event_linked_frame_anchor_in_scanned_inputs"


def test_exon_clip_link_requires_the_exon_alignments_own_clipped_bases():
    record = dict(query_sequence="AACCGGTTGATCAGTA", query_aligned_interval=[0, 8])
    assert _aligned_exon_clip_contains(record, "TTGATCAG")
    assert not _aligned_exon_clip_contains(record, "AACCGGTT")
    assert not _aligned_exon_clip_contains(
        dict(query_sequence="AACCGGTTGATCAGTA", query_aligned_interval=[0, 16]),
        "TTGATCAG")
