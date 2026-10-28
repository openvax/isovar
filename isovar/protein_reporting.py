"""Shared report summaries of exported protein hypotheses and their paths."""

from collections.abc import Mapping, Sequence
from typing import Any

from .protein_window_comparison import ComparisonStatus


def comparison_status(comparisons: Sequence[Mapping[str, Any]]) -> ComparisonStatus:
    """Require compatible reference comparisons to agree; missing data stay unknown."""
    states = {item.get("status", "unresolved") for item in comparisons}
    if states == {"reference_identical"}:
        return "reference_identical"
    if states == {"amino_acid_divergent"}:
        return "amino_acid_divergent"
    return "unresolved"


def protein_comparison_status(protein: Mapping[str, Any], *, focal: bool = False) -> ComparisonStatus:
    """Classify one amino-acid hypothesis without pooling conflicting reference models."""
    key = "focal_allele_comparison" if focal else "annotation_comparison"
    comparisons = []
    for translation in protein["translations"]:
        record = translation.get("reference_comparison", {})
        comparison = record.get(key, {})
        # Older native receipts already recorded the focal counterfactual,
        # but do not establish annotation identity or changed-residue counts.
        if not comparison and focal and isinstance(record.get("focal_amino_acid_change"), bool):
            comparison = dict(status="amino_acid_divergent" if record["focal_amino_acid_change"] else
                              "reference_identical")
        comparisons.append(comparison)
    return comparison_status(comparisons)


def best_supported_translation(protein: Mapping[str, Any]) -> dict[str, Any]:
    """Choose one display path by focal change, exact coverage, then stable identity.

    This is presentation selection, not reconciliation of competing models.
    Reported coverage belongs to this nucleotide path, never their sum.
    """
    return max(protein["translations"], key=lambda translation: (
        translation["reference_comparison"]["focal_amino_acid_change"] is True,
        translation["coding_support"]["allele_observation_span_support"]["fragments"],
        translation["coding_support"]["minimum_fragment_depth"], translation["translation_id"]))


def protein_report_record(protein: Mapping[str, Any], call_id: str) -> dict[str, Any]:
    """Flatten report metrics while retaining distinct support and comparison units.

    The native export remains the authority for each path, interval and witness
    set. An exact count of changed residues requires all compatible comparisons
    to agree, even when every comparison establishes some divergence.
    """
    selected = best_supported_translation(protein)
    coverage = selected["coding_support"]
    comparisons = [t["reference_comparison"].get("annotation_comparison", {})
                   for t in protein["translations"]]
    changed_counts = {row.get("changed_residue_count") for row in comparisons}
    state = protein_comparison_status(protein)
    return dict(
        call_id=call_id, rank=str(protein["isovar_rank"]), protein_sequence=protein["amino_acids"],
        mutation_start_idx=str(protein["mutation_interval"][0]),
        mutation_end_idx=str(protein["mutation_interval"][1]),
        supporting_reads=str(protein["rna_support"]["reads"]),
        contributor_fragments=protein["rna_support"]["fragments"],
        comparison_status=state, focal_comparison_status=protein_comparison_status(protein, focal=True),
        changed_residue_count=next(iter(changed_counts))
            if len(changed_counts) == 1 and state != "unresolved" else None,
        displayed_translation_id=selected["translation_id"],
        coding_interval=coverage["coding_interval"],
        window_fragments=coverage["allele_observation_span_support"]["fragments"],
        single_read_fragments=coverage["single_read_span_support"]["fragments"],
        paired_template_fragments=coverage.get("paired_template_span_support", {}).get("fragments"),
        minimum_fragment_depth=coverage["minimum_fragment_depth"])
