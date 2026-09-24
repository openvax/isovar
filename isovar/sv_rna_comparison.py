"""Compare RNA-derived protein hypotheses with source-identified predictions.

Exact sequence comparisons do not establish initiation, translation, event
causation, or identity of isoforms. Partial RNA is never called a full protein.
"""
from copy import deepcopy
from hashlib import sha256

from .sv_rna_orf_export import export_sv_rna_orfs
from .sv_rna_relations import EVENT_CROSSING_RELATIONS

EVENT_RELATIONS = frozenset(EVENT_CROSSING_RELATIONS)


def _matches(observed, complete, prediction):
    expected = prediction.get("amino_acids")
    if not expected:
        return dict(status="prediction_has_no_protein", prediction_id=prediction["prediction_id"], intervals=[])
    positions, offset = [], expected.find(observed)
    while offset >= 0:
        positions.append([offset, offset + len(observed)])
        offset = expected.find(observed, offset + 1)
    if observed == expected:
        status = "same_amino_acid_sequence"
    elif positions and not complete:
        status = "rna_fragment_matches_prediction"
    else:
        status = "no_exact_sequence_match"
    return dict(prediction_id=prediction["prediction_id"], status=status, intervals=positions,
                scope="complete_candidate" if complete and prediction.get("complete", False) else "partial_or_unknown")


def compare_sv_rna_predictions(result, predictions):
    """Compare every event-linked RNA translation/ORF to supplied predictions.

    ``predictions`` is a mapping with matching ``event_id`` and
    ``reference_name``, nonempty ``source`` provenance, and a ``predictions``
    list. Each entry has a unique ``prediction_id``, ``amino_acids`` (or None
    for unresolved), optional ``complete`` (default False), and any additional
    metadata, preserved verbatim. Optional sample_id must match the RNA result.

    Returns a comparison mapping suitable for inclusion in the reconstruction
    JSON. Identical nucleotide/protein hypotheses share path IDs and original
    witness identities. Fragment counts are deduplicated within this source;
    junction support never substitutes for a full ORF witness.
    """
    if not predictions.get("source"):
        raise ValueError("Prediction source provenance is required")
    for key in ("event_id", "reference_name"):
        if not result.get(key) or result[key] != predictions.get(key):
            raise ValueError("RNA/prediction %s mismatch" % key)
    if predictions.get("sample_id") not in (None, result.get("sample_id")):
        raise ValueError("RNA/prediction sample_id mismatch")
    ids = set()
    for prediction in predictions["predictions"]:
        identifier = prediction.get("prediction_id")
        if not identifier or identifier in ids:
            raise ValueError("Prediction IDs must be nonempty and unique")
        ids.add(identifier)
        sequence = prediction.get("amino_acids")
        if sequence is not None and (not isinstance(sequence, str) or not sequence
                                     or set(sequence) - set("ACDEFGHIKLMNPQRSTVWYOU")):
            raise ValueError("Predicted proteins must be explicit uppercase amino acids or None")
        if type(prediction.get("complete", False)) is not bool:
            raise ValueError("Prediction completeness must be boolean")
    # Filter occurrences before the shared exporter unions their witnesses.
    # Regional-only occurrences must not inflate event-linked support.
    # exploratory_orfs already selects these relations; this guards results
    # read from files, which may be older or edited.
    linked = dict(result, paths=[])
    for path in result["paths"]:
        candidates = [c for c in path["exploratory_orfs"]["candidates"]
                      if any(path["junctions"][i]["relation"] in EVENT_RELATIONS
                             for i in c["crossed_junctions"])]
        if candidates:
            linked["paths"].append(dict(path, exploratory_orfs=dict(path["exploratory_orfs"], candidates=candidates)))
    rows = []
    if linked["paths"]:
        for candidate in export_sv_rna_orfs(linked)["candidates"]:
            complete = candidate["ends_with_stop_codon"]
            rows.append(dict(hypothesis_id=candidate["candidate_id"], kind="exploratory_orf",
                amino_acids=candidate["amino_acids"], nucleotide_sequence=candidate["nucleotide_sequence"],
                complete_candidate=complete, initiation_observed=False, translation_observed=False,
                full_interval_fragments=candidate["rna_support"]["fragments"], orf=candidate,
                comparisons=[_matches(candidate["amino_acids"], complete, p) for p in predictions["predictions"]]))
    annotated = {}
    for path in result["paths"]:
        for candidate in path["translations"]:
            if not EVENT_RELATIONS.intersection(candidate["departure_relations"]):
                continue
            start, end = candidate["translation_start"], candidate["translation_end"]
            complete = candidate["ends_with_stop_codon"] and any(
                e["complete_5prime"] for e in candidate["frame_evidence"])
            sequence, nucleotide = candidate["amino_acids"], path["sequence"][start:end]
            key = nucleotide, sequence, complete
            item = annotated.setdefault(key, dict(kind="annotated_frame", amino_acids=sequence,
                nucleotide_sequence=nucleotide, complete_candidate=complete,
                initiation_observed=False, translation_observed=False, paths={},
                full_interval_fragments=None,
                comparisons=[_matches(sequence, complete, p) for p in predictions["predictions"]]))
            item["paths"][path["path_id"]] = dict(linkage=path["event_linkage"]["status"],
                                                 query_interval=[start, end], candidate=deepcopy(candidate))
    for key, item in sorted(annotated.items()):
        item["hypothesis_id"] = sha256(repr((result["event_id"], result["reference_name"], key)).encode()).hexdigest()
        rows.append(item)
    return dict(schema="isovar.sv_rna_prediction_comparison.v2", event_id=result["event_id"],
                reference_name=result["reference_name"], sample_id=result.get("sample_id"),
                rna_source=result["source"], predictions=deepcopy(predictions), hypotheses=rows,
                reconstruction_limitations=deepcopy(result["limitations"]),
                status="sequence_comparisons" if rows else "no_event_linked_protein_hypotheses",
                translation_observed=False, somatic_causation_proven=False)
