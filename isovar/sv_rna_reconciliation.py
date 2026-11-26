"""Competitive, coverage-aware reconciliation of nominated RNA hypotheses.

Comparisons concern observed nucleotide intervals, not absence of a transcript
from a sample. Conflicting reads can reflect coexisting isoforms. Predictions
without sequence and unexamined alternatives remain explicit.
"""

from collections.abc import Mapping
from copy import deepcopy
from hashlib import sha256


def _normalize_junction(entry):
    contig, start, end, strand = entry
    start_i, end_i = int(start), int(end)
    if start_i > end_i:
        start_i, end_i = end_i, start_i
    return (str(contig), start_i, end_i, str(strand))


def _materialize_open_template(model, sequence_provider):
    if model is None or sequence_provider is None:
        return None
    template_fn = getattr(model, "template_sequence", None)
    if not callable(template_fn):
        return None
    model_evidence = getattr(model, "evidence", None)
    if isinstance(model_evidence, Mapping) and not model_evidence.get("template_segments"):
        return None
    if not isinstance(model_evidence, Mapping) and hasattr(model, "evidence"):
        return None
    template = template_fn(sequence_provider)
    return None if template is None else str(template).upper()


def hypotheses_from_varcode(effects, sequence_provider=None):
    """Preserve all Varcode EffectCandidates, including unresolved products.

    Varcode owns consequence prediction. This adapter reads its materialized
    ``MutantTranscript`` and candidate metadata (``mechanisms``, ``phase``,
    ``uncertainty``, 1-based inclusive ascending ``(contig, start, end, strand)``
    junctions, ``mate_region``, and ``exit_frame``); it neither chooses the
    first effect nor reconstructs missing cDNA from protein or partial
    reference segments.

    When ``sequence_provider`` is supplied (a callable
    ``(contig, start_1based, end_1based) -> str``), open rearrangement
    templates (``not model.is_complete`` with template segments) materialize
    their oriented 5'-to-3' pre-mRNA sequence via
    ``model.template_sequence(sequence_provider)`` while retaining
    ``is_complete=False``. Without ``sequence_provider``, open templates keep
    ``sequence=None``.
    """
    import varcode

    rows, seen = [], set()
    for effect in effects:
        candidates = getattr(effect, "candidates", ())
        if not candidates:
            candidates = (varcode.EffectCandidate(effect=effect),)
        for candidate in candidates:
            inner = candidate.effect
            model = getattr(inner, "mutant_transcript", None)
            transcript = getattr(inner, "transcript", None)
            cdna_sequence = getattr(model, "cdna_sequence", None)
            is_complete = bool(getattr(model, "is_complete", cdna_sequence is not None))
            sequence = cdna_sequence
            if sequence is None and not is_complete and sequence_provider is not None:
                sequence = _materialize_open_template(model, sequence_provider)
            protein = getattr(inner, "mutant_protein_sequence", None)
            evidence = deepcopy(dict(getattr(candidate, "evidence", None) or {}))
            mechanisms = tuple(
                str(m) for m in (getattr(candidate, "mechanisms", None) or ())
            )
            phase = tuple(
                (str(variant), str(state))
                for variant, state in (getattr(candidate, "phase", None) or ())
            )
            uncertainty = tuple(
                str(axis) for axis in (getattr(candidate, "uncertainty", None) or ())
            )
            novel_junctions = tuple(
                _normalize_junction(j)
                for j in (getattr(model, "novel_junctions", None) or ())
            )
            back_junctions = tuple(
                _normalize_junction(j)
                for j in (getattr(model, "back_junctions", None) or ())
            )
            wildtype_junctions = tuple(
                _normalize_junction(j)
                for j in (getattr(model, "wildtype_junctions", None) or ())
            )
            raw_mate_region = getattr(model, "mate_region", None)
            mate_region = (
                _normalize_junction(raw_mate_region)
                if raw_mate_region is not None
                else None
            )
            raw_exit_frame = getattr(model, "exit_frame", None)
            exit_frame = None if raw_exit_frame is None else int(raw_exit_frame)
            row = dict(
                kind="varcode_prediction",
                sequence=sequence,
                amino_acids=protein,
                transcript_id=getattr(transcript, "id", None),
                effect_type=type(inner).__name__,
                description=inner.short_description,
                source=dict(
                    producer=getattr(candidate, "source", "varcode"),
                    varcode_version=varcode.__version__,
                ),
                evidence=evidence,
                mechanisms=mechanisms,
                phase=phase,
                uncertainty=uncertainty,
                is_complete=is_complete,
                novel_junctions=novel_junctions,
                back_junctions=back_junctions,
                wildtype_junctions=wildtype_junctions,
                mate_region=mate_region,
                exit_frame=exit_frame,
            )
            # Ordinals retain genuinely separate hypotheses even when proteins
            # coincide; do not discard their phase/splice evidence.
            identifier = sha256(repr(row).encode()).hexdigest()
            if identifier not in seen:
                rows.append(dict(hypothesis_id="varcode:" + identifier, **row))
                seen.add(identifier)
    return rows


def _unique_offset(sequence, fragment):
    start = sequence.find(fragment)
    return start if start >= 0 and sequence.find(fragment, start + 1) < 0 else None


def _compare_window(sequence, window, flank):
    if not sequence:
        return "unassessed", None
    start = sequence.find(window)
    if start >= 0:
        return "compatible", [start, start + len(window)]
    # Only two unique, collinear exact flanks make the intervening difference
    # informative. A missing match, one flank or repeat placement is unknown.
    left = _unique_offset(sequence, window[:flank])
    right = _unique_offset(sequence, window[-flank:])
    if left is not None and right is not None and right >= left + flank:
        return "contradicted_by_informative_evidence", [left, right + flank]
    return "unassessed", None


def reconcile_sv_rna_hypotheses(assessment, hypotheses=(), *, include_references=True):
    """Return every candidate and the original fragments discriminating it.

    A candidate can have both support and contradiction. Neither eliminates it:
    multiple transcripts can coexist. ``supported`` only applies to the listed
    local windows, never to initiation, a full isoform or somatic causation.
    """
    candidates = [deepcopy(h) for h in hypotheses]
    for path in assessment["paths"]:
        candidates.append(dict(hypothesis_id="rna:" + path["path_id"], kind="reconstructed_rna",
                               sequence=path["sequence"], path_id=path["path_id"],
                               source=assessment["source"]))
    if include_references:
        for reference in assessment["reference_models"]:
            candidates.append(dict(hypothesis_id="reference:" + reference["annotation"] + ":" + reference["transcript_id"],
                                   kind="ordinary_transcript", sequence=reference["sequence"],
                                   transcript_id=reference["transcript_id"],
                                   source=dict(reference_name=reference["reference_name"],
                                               annotation=reference["annotation"])))
    for placement in assessment.get("alternative_placements", []):
        candidates.append(dict(hypothesis_id="placement:" + placement["placement_id"],
                               kind="alternative_original_alignment", sequence=placement["sequence"],
                               source=dict(source_id=assessment["source"], alignment=placement["original_alignment"]),
                               placement=deepcopy(placement), scope="local_alignment_only"))
    ids = [h["hypothesis_id"] for h in candidates]
    if len(set(ids)) != len(ids):
        raise ValueError("Hypothesis IDs must be unique across supplied, RNA and reference candidates")
    for candidate in candidates:
        sequence = candidate.get("sequence")
        if not candidate.get("source") or not candidate["hypothesis_id"]:
            raise ValueError("Every hypothesis needs source provenance and an identity")
        if sequence is not None and set(sequence) - set("ACGT"):
            raise ValueError("Hypotheses require uppercase DNA, a predicted empty transcript, or unresolved null")
        if candidate.get("reference_name", assessment["reference_name"]) != assessment["reference_name"]:
            raise ValueError("Hypothesis reference assembly mismatch")
    flank = assessment["parameters"]["junction_flank"]
    windows = {}
    for path in assessment["paths"]:
        for junction in path["junctions"]:
            window = junction["window_sequence"]
            if len(window) < 2 * flank:
                continue
            entry = windows.setdefault(window, dict(path_ids=set(), fragments=set(), relations=set()))
            entry["path_ids"].add(path["path_id"])
            entry["relations"].add(junction["relation"])
            entry["fragments"].update(tuple(f) for f in junction["exact_sequence_support"]["complete_fragments"])
        # Local branch windows can distinguish RNA paths that share a join.
        branch = path.get("local_branch_support") or {}
        interval = branch.get("query_interval")
        if interval and interval[1] - interval[0] >= 2 * flank:
            window = path["sequence"][interval[0]:interval[1]]
            entry = windows.setdefault(window, dict(path_ids=set(), fragments=set(), relations=set()))
            entry["path_ids"].add(path["path_id"])
            entry["relations"].add("local_branch")
            entry["fragments"].update(tuple(m["fragment_id"]) for m in path["exact_matches"]
                                       if m["path_interval"][0] <= interval[0]
                                       and m["path_interval"][1] >= interval[1])
    observations, matrix = [], {identifier: [] for identifier in ids}
    for window, evidence in sorted(windows.items()):
        identifier = sha256(window.encode()).hexdigest()
        comparisons = [_compare_window(h.get("sequence"), window, flank) for h in candidates]
        compatible = [h["hypothesis_id"] for h, (status, _) in zip(candidates, comparisons)
                      if status == "compatible"]
        contradicted = [h["hypothesis_id"] for h, (status, _) in zip(candidates, comparisons)
                        if status == "contradicted_by_informative_evidence" and evidence["fragments"]]
        observations.append(dict(evidence_id=identifier, sequence=window,
                                 fragment_ids=[list(f) for f in sorted(evidence["fragments"])],
                                 path_ids=sorted(evidence["path_ids"]), relations=sorted(evidence["relations"]),
                                 compatible_hypothesis_ids=compatible, contradicted_hypothesis_ids=contradicted))
        for hypothesis, (status, interval) in zip(candidates, comparisons):
            if not evidence["fragments"]:
                status = "insufficient_RNA"
            elif status == "compatible":
                status = "supported" if contradicted else "indistinguishable_in_observed_RNA"
            matrix[hypothesis["hypothesis_id"]].append(dict(evidence_id=identifier, status=status,
                                                         hypothesis_interval=interval))
    rows = []
    for hypothesis in candidates:
        comparisons = matrix[hypothesis["hypothesis_id"]]
        states = {row["status"] for row in comparisons}
        if not hypothesis.get("sequence"):
            status = "unassessed"
        elif "supported" in states and "contradicted_by_informative_evidence" in states:
            status = "mixed_support_and_contradiction"
        elif "supported" in states:
            status = "supported"
        elif "contradicted_by_informative_evidence" in states:
            status = "contradicted_by_informative_evidence"
        elif "indistinguishable_in_observed_RNA" in states:
            status = "indistinguishable_in_observed_RNA"
        else:
            status = "insufficient_RNA"
        rows.append(dict(hypothesis=hypothesis, status=status, comparisons=comparisons,
                         sequence_state=("unavailable" if hypothesis.get("sequence") is None else
                                         "predicted_empty_transcript" if hypothesis["sequence"] == "" else "materialized")))
    return dict(schema="isovar.sv_rna_reconciliation.v1", event_id=assessment["event_id"],
                reference_name=assessment["reference_name"], sample_id=assessment["sample_id"],
                source=assessment["source"], hypotheses=rows, evidence=observations,
                limitations=deepcopy(assessment["limitations"]),
                scope="local_observed_intervals; coexisting_hypotheses_retained",
                translation_observed=False, somatic_causation_proven=False)
