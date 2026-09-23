"""Portable exploratory SV ORFs with source-scoped, deduplicated RNA evidence."""

from collections import defaultdict
from copy import deepcopy
import csv
from hashlib import sha256
import json
from pathlib import Path
from urllib.parse import quote

from .cell_umi import summarize_cell_umi_rows
from .default_parameters import (
    SV_MIN_ALTERNATIVE_FRAGMENTS, SV_MIN_ALTERNATIVE_FRACTION, SV_MIN_LOCAL_VARIANT_FRACTION,
)
from .genetic_code import standard_genetic_code
from .orf_start import summarize_orf_start_evidence
from .read_lineage import summarize_lineage_rows


def normalize_sv_rna_orf_export(export):
    """Copy a v1 or v2 ORF export into the v2 warning vocabulary.

    Parameters
    ----------
    export : dict
        An ``isovar.sv_rna_orfs.v1`` or ``isovar.sv_rna_orfs.v2`` mapping.

    Returns
    -------
    dict
        A deep copy with schema v2 and each legacy orientation flag renamed
        once. Unknown flags and every other field are preserved. Evidence,
        sequence IDs and counts are not recomputed. Unrecognized schemas
        raise ``ValueError``. This is also the writer's migration policy.
    """
    if export.get("schema") not in ("isovar.sv_rna_orfs.v1", "isovar.sv_rna_orfs.v2"):
        raise ValueError("Expected isovar.sv_rna_orfs.v1 or isovar.sv_rna_orfs.v2")
    aliases = dict(reverse_complement_query_witnesses_only="reverse_complement_support_only",
                   rna_polarity_unresolved="rna_strand_unresolved",
                   mixed_query_orientations="mixed_read_orientations")
    result = deepcopy(export)
    result["schema"] = "isovar.sv_rna_orfs.v2"
    for candidate in result["candidates"]:
        candidate["uncertainty_flags"] = sorted({aliases.get(flag, flag)
                                                 for flag in candidate["uncertainty_flags"]})
    return result


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _identifier(kind, value):
    return kind + "_" + sha256(_canonical([kind, value]).encode()).hexdigest()


def _unique(rows):
    return [json.loads(row) for row in sorted({_canonical(row) for row in rows})]


def _support(result, witnesses):
    scope = [result["sample_id"], result["source"]]
    observations = [result["observations"][w["observation"]] for w in witnesses]
    identities = sorted({tuple(o["identity"]) for o in observations})
    fragments = sorted({identity[:2] for identity in identities})
    segment_ids = [_identifier("segment", [scope, identity]) for identity in identities]
    fragment_ids = [_identifier("fragment", [scope, identity]) for identity in fragments]
    labels = {tuple(row["identity"]): row for row in result.get("cell_umi_evidence", {}).get("segments", [])}
    lineage = {tuple(row["identity"]): row for row in result.get("read_lineage", {}).get("segments", [])}
    label_rows = [labels.get(identity, dict(identity=list(identity), label=None, library_scope_known=False,
                                           status="metadata_unavailable")) for identity in identities]
    lineage_rows = [lineage.get(identity, dict(identity=list(identity), signal_group=None,
                                              status="metadata_unavailable")) for identity in identities]
    label_summary = summarize_cell_umi_rows(label_rows)
    lineage_summary = summarize_lineage_rows(lineage_rows)
    # These legacy summaries carry raw RG/QNAME identities. Membership below
    # uses explicit source-scoped fingerprints instead.
    label_summary.pop("segment_ids")
    lineage_summary.pop("segment_ids")
    signal_ids = sorted({_identifier("signal", [scope, row["signal_group"]]) for row in lineage_rows
                         if row["signal_group"] is not None})
    orientations = defaultdict(set)
    for o in observations:
        orientations[tuple(o["identity"][:2])].add(o["reverse_complement"])
    orientation_counts = dict(original_query=0, reverse_complement=0, mixed=0)
    for directions in orientations.values():
        key = "mixed" if len(directions) > 1 else "reverse_complement" if True in directions else "original_query"
        orientation_counts[key] += 1
    return dict(
        scope="union_of_full_interval_junction_linked_witnesses_within_input",
        segments=len(identities), fragments=len(fragments), segment_ids=segment_ids, fragment_ids=fragment_ids,
        evidence_set_id=_identifier("evidence", sorted(segment_ids)),
        signal_group_ids=signal_ids, cell_umi_support=label_summary, read_lineage=lineage_summary,
        independent_molecules=None, fragment_query_orientations=orientation_counts,
        missing_quality_segments=len({tuple(o["identity"]) for o in observations if o["missing_qualities"]}))


def _occurrence(path, candidate):
    row = {key: deepcopy(candidate[key]) for key in (
        "query_interval", "annotated_start", "start_context", "reference_comparisons", "crossed_junctions")}
    row.update({key: deepcopy(path[key]) for key in (
        "path_id", "frame_status", "unresolved_models")})
    row["reconstruction_scopes"] = deepcopy(path.get("reconstruction_scopes", []))
    row["start_evidence"] = deepcopy(candidate.get("start_evidence"))
    row["junctions"] = []
    crossings = {c["junction_index"]: c for c in candidate.get("junction_crossings", [])}
    for index in candidate["crossed_junctions"]:
        junction = path["junctions"][index]
        geometry = {key: deepcopy(junction[key]) for key in (
            "query_interval", "left", "right", "unplaced_bases", "forward_splice_geometry",
            "kinds", "annotated", "relation", "breakpoint_assignment") if key in junction}
        geometry.update(deepcopy(crossings.get(index, dict(junction_index=index, boundaries=None,
                                                          termination_only=None))))
        row["junctions"].append(geometry)
    row["witnesses"] = _unique(candidate["full_interval_support"]["witnesses"])
    row["orf_candidate_limit_reached"] = path["exploratory_orfs"]["candidate_limit_reached"]
    return row


def _flags(candidate, result):
    flags = {"initiation_unobserved", "translation_unobserved", "peptide_novelty_unassessed",
             "rna_strand_unresolved", "interval_base_quality_unassessed"}
    support, occurrences = candidate["rna_support"], candidate["occurrences"]
    if candidate["start_evidence_summary"]["status"] != "consistent":
        flags.add("start_tier_" + candidate["start_evidence_summary"]["status"])
    if not any(o["annotated_start"] for o in occurrences):
        flags.add("no_annotated_start_in_supplied_models")
    if any(o["frame_status"] != "translated" for o in occurrences):
        flags.add("annotated_path_frame_unresolved_or_ambiguous")
    if not candidate["ends_with_stop_codon"]:
        flags.add("partial_orf")
    if support["fragments"] < 2:
        flags.add("no_full_fragment_witness" if not support["fragments"] else "single_full_fragment_witness")
    if support["missing_quality_segments"]:
        flags.add("missing_base_qualities")
    labels, lineage = support["cell_umi_support"], support["read_lineage"]
    if labels["unknown_library_segments"]:
        flags.add("library_scope_unresolved")
    if labels["unresolved_segments"]:
        flags.add("cell_umi_labels_unresolved")
    if lineage["unresolved_segments"]:
        flags.add("signal_lineage_unresolved")
    if support["segments"] - lineage["unresolved_segments"] > lineage["resolved_signal_groups"]:
        flags.add("shared_signal_ancestry")
    directions = support["fragment_query_orientations"]
    if support["fragments"] and directions["reverse_complement"] == support["fragments"]:
        flags.add("reverse_complement_support_only")
    if directions["mixed"] or directions["original_query"] and directions["reverse_complement"]:
        flags.add("mixed_read_orientations")
    for occurrence in occurrences:
        if occurrence["orf_candidate_limit_reached"]:
            flags.add("orf_candidate_limit_reached")
        for junction in occurrence["junctions"]:
            if junction.get("unplaced_bases"):
                flags.add("unplaced_junction_sequence")
            if junction["termination_only"]:
                flags.add("termination_only_junction_crossing")
            if junction["boundaries"] is None:
                flags.add("junction_boundary_classification_unavailable")
            if junction["relation"] in ("splice_ambiguous_event_junction", "breakpoint_clip_partner_unplaced"):
                flags.add(junction["relation"])
    flags.update("reconstruction_limit:" + note for note in result["limitations"])
    defaults = dict(min_alternative_fragments=SV_MIN_ALTERNATIVE_FRAGMENTS,
                    min_alternative_fraction=SV_MIN_ALTERNATIVE_FRACTION,
                    min_local_variant_fraction=SV_MIN_LOCAL_VARIANT_FRACTION)
    if any(key in result["parameters"] and result["parameters"][key] != value for key, value in defaults.items()):
        flags.add("nondefault_branch_thresholds")
    return sorted(flags)


def export_sv_rna_orfs(result):
    """Export every exploratory ATG ORF from one SV reconstruction.

    Parameters
    ----------
    result : dict
        An ``isovar.sv_rna_candidates.v2`` result. Not modified.

    Returns
    -------
    dict
        ``isovar.sv_rna_orfs.v2`` sequences, occurrences, evidence and warnings.
        Exact nucleotide alternatives remain distinct even if their peptides
        agree. Candidate identity excludes sample/source; evidence identity
        includes them and excludes event. Hashes are references, not guarantees
        of anonymization or evidence independence across source aliases.
    """
    if result.get("schema") != "isovar.sv_rna_candidates.v2":
        raise ValueError("SV ORF export requires isovar.sv_rna_candidates.v2")
    groups = {}
    for path in result["paths"]:
        for candidate in path["exploratory_orfs"]["candidates"]:
            start, end = candidate["query_interval"]
            dna = path["sequence"][start:end]
            aa, stop = candidate["amino_acids"], candidate["ends_with_stop_codon"]
            if (not 0 <= start < end <= len(path["sequence"]) or not dna.startswith("ATG")
                    or len(dna) != 3 * (len(aa) + int(stop))
                    or standard_genetic_code.translate(dna, first_codon_is_start=True) != (aa, stop)):
                raise ValueError("ORF sequence, interval and translation disagree")
            key = _identifier("sv_orf", [result["event_id"], result["reference_name"], dna, aa, stop])
            if key not in groups:
                groups[key] = dict(candidate_id=key, nucleotide_sequence=dna, amino_acids=aa,
                                   nucleotide_sequence_id=_identifier("nt", dna), protein_sequence_id=_identifier("aa", aa),
                                   ends_with_stop_codon=stop, initiation_observed=False, translation_observed=False,
                                   peptide_novelty_assessed=False, occurrences=[])
            groups[key]["occurrences"].append(_occurrence(path, candidate))
    candidates = []
    for key in sorted(groups):
        candidate = groups[key]
        candidate["occurrences"] = _unique(candidate["occurrences"])
        witnesses = [w for occurrence in candidate["occurrences"] for w in occurrence["witnesses"]]
        candidate["rna_support"] = _support(result, witnesses)
        candidate["start_evidence_summary"] = summarize_orf_start_evidence(
            o["start_evidence"] for o in candidate["occurrences"])
        candidate["uncertainty_flags"] = _flags(candidate, result)
        candidates.append(candidate)
    export = {key: deepcopy(result[key]) for key in (
        "event_id", "reference_name", "sample_id", "source", "event_provenance", "donor", "acceptor",
        "reference_models", "parameters", "limitations")}
    export.update(schema="isovar.sv_rna_orfs.v2", reconstruction_schema=result["schema"], candidates=candidates,
                  interval_conventions=dict(orf_and_witness="zero_based_half_open",
                                            junction_query_interval="zero_based_flanking_base_offsets"),
                  evidence_identity_policy="sha256_of_domain_and_canonical_json; sample/source scoped; event excluded",
                  interpretation="Exploratory sequences; no inferred initiation, translation, presentation, "
                                 "independent molecules, abundance or independence across source aliases.")
    return export


def sv_rna_orf_output_paths(prefix):
    """Append format suffixes, preserving any dots in the supplied prefix."""
    return {kind: Path(str(Path(prefix).expanduser()) + suffix) for kind, suffix in (
        ("json", ".json"), ("tsv", ".tsv"), ("protein", ".protein.fasta"), ("nucleotide", ".nucleotide.fasta"))}


def write_sv_rna_orfs(export, prefix):
    """Write an exported result as JSON, TSV, protein FASTA and DNA FASTA.

    Unknown numeric values remain null in JSON and blank in TSV. DNA includes
    an observed terminal stop; protein excludes it. All alternatives are kept.
    Accepts v1/v2 exports and writes v2 using :func:`normalize_sv_rna_orf_export`;
    the input is not modified. Returns the output paths keyed by format.
    """
    export = normalize_sv_rna_orf_export(export)
    paths = sv_rna_orf_output_paths(prefix)
    for path in paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    paths["json"].write_text(json.dumps(export, indent=2) + "\n")
    columns = ["schema", "candidate_id", "event_id", "reference_name", "sample_id", "source", "amino_acids",
               "nucleotide_sequence", "ends_with_stop_codon", "fragments", "segments", "observed_cell_umi_labels",
               "complete_cell_umi_labels", "resolved_signal_groups", "independent_molecules",
               "original_query_fragments", "reverse_complement_fragments", "mixed_orientation_fragments",
               "missing_quality_segments", "evidence_set_id", "uncertainty_flags",
               "start_tier_status", "start_tier", "start_priority", "start_evidence_json"]
    with paths["tsv"].open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for candidate in export["candidates"]:
            support = candidate["rna_support"]
            row = {key: export[key] for key in ("schema", "event_id", "reference_name", "sample_id", "source")}
            row.update({key: candidate[key] for key in (
                "candidate_id", "amino_acids", "nucleotide_sequence", "ends_with_stop_codon")})
            row.update({key: support[key] for key in (
                "fragments", "segments", "independent_molecules", "missing_quality_segments", "evidence_set_id")})
            row.update(observed_cell_umi_labels=support["cell_umi_support"]["observed_labels"],
                       complete_cell_umi_labels=support["cell_umi_support"]["complete_label_count"],
                       resolved_signal_groups=support["read_lineage"]["resolved_signal_groups"],
                       original_query_fragments=support["fragment_query_orientations"]["original_query"],
                       reverse_complement_fragments=support["fragment_query_orientations"]["reverse_complement"],
                       mixed_orientation_fragments=support["fragment_query_orientations"]["mixed"],
                       uncertainty_flags=";".join(candidate["uncertainty_flags"]))
            summary = candidate.get("start_evidence_summary", summarize_orf_start_evidence([]))
            row.update(start_tier_status=summary["status"], start_tier=summary["tier"],
                       start_priority=summary["priority"],
                       start_evidence_json=_canonical([dict(path_id=o["path_id"],
                                                           start_evidence=o.get("start_evidence"))
                                                       for o in candidate["occurrences"]]))
            writer.writerow(row)
    for kind, field in (("protein", "amino_acids"), ("nucleotide", "nucleotide_sequence")):
        with paths[kind].open("w") as handle:
            for candidate in export["candidates"]:
                handle.write(">" + candidate["candidate_id"] + " event=" + quote(export["event_id"], safe="") + "\n")
                sequence = candidate[field]
                handle.write("\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + "\n")
    return paths
