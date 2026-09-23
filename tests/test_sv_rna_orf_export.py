"""Export fidelity and evidence denominators, using only synthetic sequences."""

from copy import deepcopy
import csv
import json

import pytest

from isovar import export_sv_rna_orfs, write_sv_rna_orfs
from isovar.cli.isovar_sv_rna import run as cli_run
from tests.test_sv_rna_orfs import inputs, insertion_inputs, run


def reconstruction(args=None):
    seq, pos, reads, junction = args or inputs()
    orfs = run(seq, pos, reads, junction)
    path = dict(path_id="path1", sequence=seq, junctions=[junction], exploratory_orfs=orfs,
                frame_status="unresolved", unresolved_models=[], reconstruction_scopes=["seed_spanning"])
    observations = {key: dict(identity=list(o.identity), reverse_complement=o.reverse,
                             missing_qualities=o.missing_qualities) for key, o in reads.items()}
    return dict(schema="isovar.sv_rna_candidates.v2", event_id="D--A", reference_name="synthetic",
                sample_id="sample", source="source", event_provenance={}, donor={}, acceptor={},
                reference_models=[], parameters={}, limitations=[], paths=[path], observations=observations)


def only(result):
    candidate, = export_sv_rna_orfs(result)["candidates"]
    return candidate


def test_export_unions_witnesses_across_paths_and_does_not_mutate_input():
    result = reconstruction()
    path = deepcopy(result["paths"][0])
    path["path_id"] = "path2"
    witnesses = path["exploratory_orfs"]["candidates"][0]["full_interval_support"]["witnesses"]
    witnesses *= 3
    result["paths"] += [path, deepcopy(path)]
    before = deepcopy(result)
    candidate = only(result)
    assert result == before
    assert len(candidate["occurrences"]) == 2
    assert all(len(o["witnesses"]) == 1 for o in candidate["occurrences"])
    support = candidate["rna_support"]
    assert support["fragments"] == support["segments"] == 1
    assert support["cell_umi_support"]["complete_label_count"] is None
    assert support["cell_umi_support"]["status_counts"] == {"metadata_unavailable": 1}
    assert "segment_ids" not in support["cell_umi_support"]
    assert "segment_ids" not in support["read_lineage"]
    assert support["independent_molecules"] is None


def test_sequence_ids_keep_synonymous_alternatives_and_scope_evidence_independently():
    result = reconstruction()
    first = only(result)
    alternative = deepcopy(result["paths"][0])
    alternative["path_id"] = "synonymous"
    alternative["sequence"] = alternative["sequence"].replace("AAA", "AAG")
    # This synthetic path has no full witnesses for its changed sequence.
    alternative["exploratory_orfs"]["candidates"][0]["full_interval_support"]["witnesses"] = []
    result["paths"].append(alternative)
    exported = export_sv_rna_orfs(result)
    assert len(exported["candidates"]) == 2
    assert len({c["protein_sequence_id"] for c in exported["candidates"]}) == 1
    assert len({c["nucleotide_sequence_id"] for c in exported["candidates"]}) == 2
    result["paths"].reverse()
    assert export_sv_rna_orfs(result) == exported
    result["paths"] = [result["paths"][1]]
    result["sample_id"] = "another_sample"
    changed = only(result)
    assert changed["candidate_id"] == first["candidate_id"]
    assert changed["rna_support"]["segment_ids"] != first["rna_support"]["segment_ids"]
    result["sample_id"] = "sample"
    result["event_id"] = "another_event"
    changed = only(result)
    assert changed["candidate_id"] != first["candidate_id"]
    assert changed["rna_support"]["segment_ids"] == first["rna_support"]["segment_ids"]
    result["source"] = "another_source"
    assert only(result)["rna_support"]["fragment_ids"] != first["rna_support"]["fragment_ids"]


def add_witness(result, key, identity, reverse=False, missing=False):
    result["observations"][key] = dict(identity=identity, reverse_complement=reverse, missing_qualities=missing)
    witnesses = result["paths"][0]["exploratory_orfs"]["candidates"][0]["full_interval_support"]["witnesses"]
    witness = deepcopy(witnesses[0])
    witness["observation"] = key
    witnesses.append(witness)


def test_mates_and_alternative_orientations_count_segments_and_fragments_separately():
    result = reconstruction()
    result["observations"]["full"]["identity"] = ["library", "read", 64]
    add_witness(result, "mate", ["library", "read", 128], reverse=True, missing=True)
    add_witness(result, "placement", ["library", "read", 64])
    candidate = only(result)
    support = candidate["rna_support"]
    assert support["segments"] == 2 and support["fragments"] == 1
    assert support["fragment_query_orientations"] == dict(original_query=0, reverse_complement=0, mixed=1)
    assert support["missing_quality_segments"] == 1
    assert {"mixed_query_orientations", "missing_base_qualities"} <= set(candidate["uncertainty_flags"])


def test_known_shared_library_labels_and_signal_descendants_use_the_shared_ledgers():
    result = reconstruction()
    add_witness(result, "other", ["other_rg", "other_read", 0])
    identities = [o["identity"] for o in result["observations"].values()]
    result["cell_umi_evidence"] = dict(segments=[dict(identity=i, library_scope_known=True,
        label=["sample", "source", "shared_library", "cell", "umi"], status="resolved_label") for i in identities])
    result["read_lineage"] = dict(segments=[dict(identity=i, signal_group=[i[0], "parent"], status="split_read")
                                           for i in identities])
    support = only(result)["rna_support"]
    assert support["cell_umi_support"]["complete_label_count"] == 1
    assert support["read_lineage"]["resolved_signal_groups"] == 2  # RG still partitions signal ancestry.
    result["read_lineage"]["segments"][1]["signal_group"] = ["library", "parent"]
    candidate = only(result)
    assert candidate["rna_support"]["fragments"] == 2
    assert candidate["rna_support"]["read_lineage"]["resolved_signal_groups"] == 1
    assert "shared_signal_ancestry" in candidate["uncertainty_flags"]
    del result["cell_umi_evidence"]["segments"][1]
    assert only(result)["rna_support"]["cell_umi_support"]["complete_label_count"] is None


def test_partial_unsupported_and_limited_hypotheses_are_exported_with_warnings():
    result = reconstruction(inputs("ATGAAAATGAAA"))
    result["paths"][0]["exploratory_orfs"]["candidates"][0]["full_interval_support"]["witnesses"] = []
    result["paths"][0]["exploratory_orfs"]["candidate_limit_reached"] = True
    result["limitations"] = ["max_paths", "repeated_genomic_position"]
    result["parameters"]["min_alternative_fragments"] = 1
    candidate = only(result)
    assert candidate["nucleotide_sequence"] == "ATGAAAATGAAA"
    assert candidate["amino_acids"] == "MKMK" and not candidate["ends_with_stop_codon"]
    assert {"partial_orf", "no_full_fragment_witness", "orf_candidate_limit_reached",
            "reconstruction_limit:max_paths", "reconstruction_limit:repeated_genomic_position",
            "nondefault_branch_thresholds"} <= set(candidate["uncertainty_flags"])
    assert candidate["rna_support"]["cell_umi_support"]["complete_label_count"] is None
    assert not candidate["initiation_observed"] and not candidate["translation_observed"]


def test_stop_only_and_reverse_complement_only_flags_preserve_junction_provenance():
    result = reconstruction(insertion_inputs("ATG" + "GCC" * 14, "TGA" + "CCC" * 3, "CCC" * 8))
    result["paths"][0]["junctions"][0]["unplaced_bases"] = "TGA" + "CCC" * 3
    result["observations"]["full"]["reverse_complement"] = True
    candidate = only(result)
    assert {"termination_only_junction_crossing", "unplaced_junction_sequence",
            "reverse_complement_query_witnesses_only"} <= set(candidate["uncertainty_flags"])
    junction, = candidate["occurrences"][0]["junctions"]
    assert junction["boundaries"] == ["donor_to_unplaced"]
    assert junction["query_interval"] == [44, 57]
    assert "direct_fragments" not in junction


def test_older_v2_boundary_metadata_is_flagged_not_invented():
    result = reconstruction()
    del result["paths"][0]["exploratory_orfs"]["candidates"][0]["junction_crossings"]
    del result["paths"][0]["reconstruction_scopes"]
    assert "junction_boundary_classification_unavailable" in only(result)["uncertainty_flags"]


def test_rejects_wrong_schema_and_inconsistent_translation():
    result = reconstruction()
    result["paths"][0]["sequence"] = "A" * 15
    with pytest.raises(ValueError, match="disagree"):
        export_sv_rna_orfs(result)
    with pytest.raises(ValueError, match="requires"):
        export_sv_rna_orfs({})


def test_writer_round_trip_and_empty_outputs(tmp_path):
    exported = export_sv_rna_orfs(reconstruction())
    paths = write_sv_rna_orfs(exported, tmp_path / "nested" / "sample.orfs")
    assert paths["json"].name == "sample.orfs.json"
    assert json.loads(paths["json"].read_text()) == exported
    with paths["tsv"].open() as handle:
        row, = csv.DictReader(handle, delimiter="\t")
    assert row["amino_acids"] == "MKMK" and row["nucleotide_sequence"] == "ATGAAAATGAAATAA"
    assert row["complete_cell_umi_labels"] == row["independent_molecules"] == ""
    assert paths["protein"].read_text().splitlines()[1:] == ["MKMK"]
    assert paths["nucleotide"].read_text().splitlines()[1:] == ["ATGAAAATGAAATAA"]
    exported["candidates"] = []
    write_sv_rna_orfs(exported, tmp_path / "nested" / "sample.orfs")
    assert paths["protein"].read_text() == paths["nucleotide"].read_text() == ""
    assert len(paths["tsv"].read_text().splitlines()) == 1


@pytest.mark.parametrize("suffix", [".json", ".tsv", ".protein.fasta", ".nucleotide.fasta"])
def test_cli_rejects_output_collision_before_opening_inputs(tmp_path, capsys, suffix):
    prefix = tmp_path / "out"
    with pytest.raises(SystemExit):
        cli_run(["--input", "nonexistent", "--bam", "nonexistent", "--output", str(prefix) + suffix,
                 "--orf-output-prefix", str(prefix)])
    assert "collides with --output" in capsys.readouterr().err
