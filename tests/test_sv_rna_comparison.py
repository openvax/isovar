"""Comparison inputs are constructed inline; partial RNA cannot verify a protein."""
from copy import deepcopy

import pytest

from isovar.sv_rna_comparison import compare_sv_rna_predictions


def inputs(sequence="MAK", stop=False, complete=False):
    candidate = dict(amino_acids=sequence, translation_start=0, translation_end=9,
                     ends_with_stop_codon=stop, frame_evidence=[dict(complete_5prime=complete)],
                     departure_relations=["breakpoint_junction"])
    result = dict(event_id="event", reference_name="test", sample_id="sample", source="reads.bam",
                  limitations=[], observations={}, paths=[dict(path_id="p", sequence="ATGGCTAAA",
                  event_linkage=dict(status="breakpoint_junction"), translations=[candidate],
                  junctions=[dict(relation="breakpoint_junction")], exploratory_orfs=dict(candidates=[]))])
    predictions = dict(event_id="event", reference_name="test", source=dict(caller="varcode", version="9.4.2"),
                       predictions=[dict(prediction_id="model", amino_acids="MMAKM", complete=True)])
    return result, predictions


def test_partial_sequence_match_does_not_claim_complete_protein():
    result, predictions = inputs()
    row, = compare_sv_rna_predictions(result, predictions)["hypotheses"]
    match, = row["comparisons"]
    assert match == dict(prediction_id="model", status="rna_fragment_matches_prediction",
                         intervals=[[1, 4]], scope="partial_or_unknown")
    assert row["full_interval_fragments"] is None  # Junction support is not full-span support.
    assert row["translation_observed"] is False


def test_complete_shorter_orf_is_not_called_a_partial_match():
    result, predictions = inputs(stop=True, complete=True)
    row, = compare_sv_rna_predictions(result, predictions)["hypotheses"]
    assert row["comparisons"][0]["status"] == "no_exact_sequence_match"
    predictions["predictions"][0]["amino_acids"] = "MAK"
    row, = compare_sv_rna_predictions(result, predictions)["hypotheses"]
    assert row["comparisons"][0]["status"] == "same_amino_acid_sequence"
    assert row["comparisons"][0]["scope"] == "complete_candidate"


def test_unresolved_prediction_is_preserved_without_invented_comparison():
    result, predictions = inputs()
    predictions["predictions"][0]["amino_acids"] = None
    comparison = compare_sv_rna_predictions(result, predictions)
    assert comparison["predictions"] == predictions
    assert comparison["hypotheses"][0]["comparisons"][0]["status"] == "prediction_has_no_protein"
    result["paths"] = []
    assert compare_sv_rna_predictions(result, predictions)["status"] == "no_event_linked_protein_hypotheses"


def test_shared_orf_witnesses_are_not_added_across_path_contexts():
    from tests.test_sv_rna_orf_export import reconstruction
    result = reconstruction()
    _, predictions = inputs()
    predictions.update(event_id=result["event_id"], reference_name=result["reference_name"])
    p = result["paths"][0]
    p.update(translations=[], event_linkage=dict(status="breakpoint_junction"))
    other = deepcopy(p)
    other["path_id"] = "other"
    other["sequence"] += "CCC"
    result["paths"].append(other)
    row, = compare_sv_rna_predictions(result, predictions)["hypotheses"]
    assert {o["path_id"] for o in row["orf"]["occurrences"]} == {"path1", "other"}
    assert row["full_interval_fragments"] == 1
    assert len(row["orf"]["rna_support"]["fragment_ids"]) == 1
    # A same-sequence regional-only occurrence cannot contribute witnesses.
    regional = deepcopy(p)
    regional["path_id"] = "regional"
    regional["junctions"][0]["relation"] = "regional_novel_junction"
    result["paths"].append(regional)
    again, = compare_sv_rna_predictions(result, predictions)["hypotheses"]
    assert again == row


@pytest.mark.parametrize("key,value", [("event_id", "other"), ("reference_name", "other"),
                                       ("sample_id", "other"), ("source", None)])
def test_mismatched_context_is_rejected(key, value):
    result, predictions = inputs()
    predictions[key] = value
    with pytest.raises(ValueError):
        compare_sv_rna_predictions(result, predictions)


def test_regional_translation_cannot_validate_an_event():
    result, predictions = inputs()
    result["paths"][0]["translations"][0]["departure_relations"] = ["regional_novel_junction"]
    assert not compare_sv_rna_predictions(result, predictions)["hypotheses"]


def test_cli_comparison_and_shared_orf_export_agree(tmp_path, monkeypatch):
    import json
    from contextlib import nullcontext
    from isovar.cli import isovar_sv_rna
    from tests.test_sv_rna_orf_export import reconstruction

    result = reconstruction()
    result["paths"][0].update(translations=[], event_linkage=dict(status="breakpoint_junction"))
    _, predictions = inputs()
    predictions.update(event_id=result["event_id"], reference_name=result["reference_name"])
    predictions["predictions"][0].update(amino_acids="MKMK", complete=True)
    event, predicted, output = (tmp_path / name for name in ("input.json", "predictions.json", "result.json"))
    event.write_text("{}")
    predicted.write_text(json.dumps(predictions))
    monkeypatch.setattr(isovar_sv_rna, "sv_rna_input_from_dict", lambda value: {})
    monkeypatch.setattr(isovar_sv_rna, "alignment_file_from_args", lambda options: nullcontext(None))
    monkeypatch.setattr(isovar_sv_rna, "reconstruct_sv_rna", lambda *args, **kwargs: result)
    isovar_sv_rna.run(["--input", str(event), "--bam", "original.bam", "--output", str(output),
                       "--predictions", str(predicted), "--orf-output-prefix", str(tmp_path / "orfs")])
    row, = json.loads(output.read_text())["prediction_comparison"]["hypotheses"]
    exported, = json.loads((tmp_path / "orfs.json").read_text())["candidates"]
    assert row["orf"] == exported
    assert row["comparisons"][0]["status"] == "same_amino_acid_sequence"
    assert row["full_interval_fragments"] == 1
