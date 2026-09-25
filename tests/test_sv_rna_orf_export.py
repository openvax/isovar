"""Export fidelity and evidence denominators, using only synthetic sequences."""

from copy import deepcopy
from tests.testing_helpers import complete_umis
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
                frame_status="unresolved", unresolved_models=[], reconstruction_scopes=["seed_spanning"],
                end_reasons={"5prime": ["observations_end"], "3prime": ["observations_end"]})
    observations = {key: dict(identity=list(o.identity), reverse_complement=o.reverse_complement,
                             missing_qualities=o.missing_qualities) for key, o in reads.items()}
    return dict(schema="isovar.sv_rna_candidates.v5", event_id="D--A", reference_name="synthetic",
                sample_id="sample", source="source", event_provenance={}, donor={}, acceptor={},
                reference_models=[], parameters={}, limitations=[], paths=[path], observations=observations,
                cell_umi_evidence=dict(policy="isovar.cell_umi_labels.v1", reads=[]), read_lineage=dict(reads=[]))


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
    assert support["fragments"] == support["reads"] == 1
    assert complete_umis(support) is None
    assert support["label_statuses"] == {"metadata_unavailable": 1}
    assert "segment_ids" not in support["read_lineage"]


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
    assert changed["rna_support"]["read_ids"] != first["rna_support"]["read_ids"]
    result["sample_id"] = "sample"
    result["event_id"] = "another_event"
    changed = only(result)
    assert changed["candidate_id"] != first["candidate_id"]
    assert changed["rna_support"]["read_ids"] == first["rna_support"]["read_ids"]
    result["source"] = "another_source"
    assert only(result)["rna_support"]["fragment_ids"] != first["rna_support"]["fragment_ids"]


def add_witness(result, key, identity, reverse_complement=False, missing=False):
    result["observations"][key] = dict(identity=identity, reverse_complement=reverse_complement, missing_qualities=missing)
    witnesses = result["paths"][0]["exploratory_orfs"]["candidates"][0]["full_interval_support"]["witnesses"]
    witness = deepcopy(witnesses[0])
    witness["observation"] = key
    witnesses.append(witness)


def test_mates_and_alternative_orientations_count_segments_and_fragments_separately():
    result = reconstruction()
    result["observations"]["full"]["identity"] = ["library", "read", 64]
    add_witness(result, "mate", ["library", "read", 128], reverse_complement=True, missing=True)
    add_witness(result, "placement", ["library", "read", 64])
    candidate = only(result)
    support = candidate["rna_support"]
    assert support["reads"] == 2 and support["fragments"] == 1
    assert support["fragment_query_orientations"] == dict(original_query=0, reverse_complement=0, mixed=1)
    assert support["missing_quality_reads"] == 1
    assert {"mixed_read_orientations", "missing_base_qualities"} <= set(candidate["uncertainty_flags"])


def test_known_shared_library_labels_and_signal_descendants_use_the_shared_ledgers():
    result = reconstruction()
    add_witness(result, "other", ["other_rg", "other_read", 0])
    identities = [o["identity"] for o in result["observations"].values()]
    scope = dict(source="source", sample_id="sample", library="shared_library")
    result["cell_umi_evidence"] = dict(reads=[dict(identity=i, library_scope_known=True, scope=scope, cell_barcode="cell",
        label=["sample", "source", "shared_library", "cell", "umi"], status="resolved_label") for i in identities])
    result["read_lineage"] = dict(reads=[dict(identity=i, signal_group=[i[0], "parent"], status="split_read")
                                        for i in identities])
    support = only(result)["rna_support"]
    assert complete_umis(support) == 1
    assert support["read_lineage"]["signal_groups"] == 2  # RG still partitions signal ancestry.
    result["read_lineage"]["reads"][1]["signal_group"] = ["library", "parent"]
    candidate = only(result)
    assert candidate["rna_support"]["fragments"] == 2
    assert candidate["rna_support"]["read_lineage"]["signal_groups"] == 1
    assert "shared_signal_ancestry" in candidate["uncertainty_flags"]
    del result["cell_umi_evidence"]["reads"][1]
    assert complete_umis(only(result)["rna_support"]) is None


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
    # No witnesses: zero reads, so zero UMIs, exactly.
    assert candidate["rna_support"]["reads"] == complete_umis(candidate["rna_support"]) == 0
    assert not candidate["initiation_observed"] and not candidate["translation_observed"]


def test_stop_only_and_reverse_complement_only_flags_preserve_junction_provenance():
    result = reconstruction(insertion_inputs("ATG" + "GCC" * 14, "TGA" + "CCC" * 3, "CCC" * 8))
    result["paths"][0]["junctions"][0]["unplaced_bases"] = "TGA" + "CCC" * 3
    result["observations"]["full"]["reverse_complement"] = True
    candidate = only(result)
    assert {"termination_only_junction_crossing", "unplaced_junction_sequence",
            "reverse_complement_support_only"} <= set(candidate["uncertainty_flags"])
    junction, = candidate["occurrences"][0]["junctions"]
    assert junction["boundaries"] == ["donor_to_unplaced"]
    assert junction["query_interval"] == [45, 57]
    assert "direct_support" not in junction


def test_missing_or_disagreeing_start_annotations_do_not_choose_the_best_occurrence(tmp_path):
    from tests.test_orf_start import annotation, reference

    result = reconstruction()
    # Two paths with the same ORF can disagree about the transcript origin.
    # This test controls their annotation to exercise export's reconciliation.
    original_id = only(result)["candidate_id"]
    other = deepcopy(result["paths"][0])
    other["path_id"] = "other"
    other["exploratory_orfs"]["candidates"][0]["start_evidence"] = annotation(reference())
    result["paths"].append(other)
    candidate = only(result)
    assert candidate["candidate_id"] == original_id
    assert candidate["start_evidence_summary"] == dict(status="ambiguous", tier=None, priority=None)
    assert "start_tier_ambiguous" in candidate["uncertainty_flags"]
    paths = write_sv_rna_orfs(export_sv_rna_orfs(result), tmp_path / "ambiguous")
    with paths["tsv"].open() as handle:
        row, = csv.DictReader(handle, delimiter="\t")
    assert row["start_tier_status"] == "ambiguous" and row["start_priority"] == row["start_tier"] == ""


def test_rejects_wrong_schema_and_inconsistent_translation():
    result = reconstruction()
    result["paths"][0]["sequence"] = "A" * 15
    with pytest.raises(ValueError, match="disagree"):
        export_sv_rna_orfs(result)
    with pytest.raises(ValueError, match="requires"):
        export_sv_rna_orfs({})
    with pytest.raises(ValueError, match="requires"):
        export_sv_rna_orfs(dict(result, schema="isovar.sv_rna_candidates.v3"))
    with pytest.raises(ValueError, match="Expected"):
        write_sv_rna_orfs(dict(schema="isovar.sv_rna_orfs.v2"), "unused")


def test_writer_round_trip_and_empty_outputs(tmp_path):
    exported = export_sv_rna_orfs(reconstruction())
    paths = write_sv_rna_orfs(exported, tmp_path / "nested" / "sample.orfs")
    assert paths["json"].name == "sample.orfs.json"
    assert json.loads(paths["json"].read_text()) == exported
    with paths["tsv"].open() as handle:
        row, = csv.DictReader(handle, delimiter="\t")
    assert row["amino_acids"] == "MKMK" and row["nucleotide_sequence"] == "ATGAAAATGAAATAA"
    # The one witness has no label row, so it is unlabelled and the UMI count is not exact.
    assert (row["reads"], row["umis"], row["umis_complete"], row["unlabeled_reads"]) == ("1", "0", "False", "1")
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


def test_path_end_reasons_are_exported_and_truncation_is_flagged():
    result = reconstruction()
    natural = only(result)
    assert natural["occurrences"][0]["end_reasons"] == result["paths"][0]["end_reasons"]
    assert "path_end_truncated" not in natural["uncertainty_flags"]
    result["paths"][0]["end_reasons"]["3prime"].append("path_limit")
    assert "path_end_truncated" in only(result)["uncertainty_flags"]
