"""COLO829 SV3 original-read evidence does not assign an unsupported frame."""

from reports.hcc1395_colo829.build_sv3_path_fusion import build


def test_source_reads_tile_both_exploratory_orfs_but_do_not_resolve_protein():
    result = build()
    assert result["schema"] == "isovar.sv3_path_fusion.v1"
    assert result["path_id"] == "653d5d9c015d3e590fb149b3"
    assert len(result["exact_junction_window_spots"]) == 5
    assert len(result["placed_split_join_spots"]) == 3
    assert sum(spot["opposite_mate_exact_path"] for spot in result["source_spots"]) == 4
    assert not any(spot["opposite_mate_exact_cnih3_exon"] for spot in result["source_spots"])
    assert result["exact_exon_anchor_bases"] == 0
    assert result["frame_status"] == "no_gene_anchor"
    assert result["event_linked_translation_candidates"] == 0
    assert "path_limit" in result["limitations"]
    short, long = result["candidates"]
    assert [short["amino_acids"], long["amino_acids"]] == [37, 162]
    assert [len(c["start_to_junction_overlap_chain"]) for c in (short, long)] == [2, 2]
    assert [len(c["full_orf_sequence_overlap_chain"]) for c in (short, long)] == [2, 7]
    assert not short["annotated_start"] and not long["annotated_start"]
    assert not short["initiation_observed"] and not long["initiation_observed"]
    assert result["conclusion"] == "mature_path_and_frame_unresolved"
