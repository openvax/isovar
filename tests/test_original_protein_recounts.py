"""Pinned original-RNA cases distinguish contributors, exact spans and silent contexts."""

import gzip
import json
from pathlib import Path

import pytest

from isovar.protein_window_comparison import compare_protein_windows


@pytest.mark.parametrize("gene,contributors,spans,single,paired,state", [
    ("CDKN2A", 5248, 1027, 996, 31, "amino_acid_divergent"),
    ("FHL2", 654, 15, 0, 15, "amino_acid_divergent"),
    ("GOLGA5", 85, 9, 2, 7, "amino_acid_divergent"),
    ("TUBB2A", 564, 141, 60, 81, "reference_identical"),
    ("C8orf76", 2, 2, 2, 0, "amino_acid_divergent"),
])
def test_original_witnesses_and_comparison(gene, contributors, spans, single, paired, state):
    path = Path(__file__).resolve().parents[1] / "reports/hcc1395_colo829/focused_proteins/protein_evidence_recount.json.gz"
    with gzip.open(path, "rt") as stream:
        case, = [row for row in json.load(stream)["cases"] if row["gene"] == gene]
    metrics = case["metrics"]
    assert (metrics["contributor_fragments"], metrics["window_fragments"], metrics["single_read_fragments"],
            metrics["paired_template_fragments"]) == (contributors, spans, single, paired)
    translation = case["displayed_translation"]
    coverage = translation["coding_support"]
    sets = []
    for key, count in (("allele_observation_span_support", spans), ("single_read_span_support", single),
                        ("paired_template_span_support", paired)):
        support = coverage[key]
        witnesses = case["evidence_sets"][support["evidence_set_id"]]
        ids = set(witnesses["fragment_ids"])
        assert support["fragments"] == len(ids) == count
        assert witnesses["evidence_scope"] == [case["sample"], case["settings"]["library"]]
        sets.append(ids)
    assert sets[1] <= sets[0] and sets[2] <= sets[0]
    comparison = translation["reference_comparison"]
    observed = metrics["protein_sequence"]
    delta = comparison["reference_allele_length"] - comparison["observed_allele_length"]
    frameshift = delta % 3 != 0
    limit = len(observed) + (delta // 3 if not frameshift else 0)
    result = compare_protein_windows(comparison["reference_protein_window"][:limit], observed,
        [int(metrics["mutation_start_idx"]), int(metrics["mutation_end_idx"])],
        allele_length_delta=delta, frameshift=frameshift,
        reference_stop=comparison['reference_stop_in_context'] and len(comparison['reference_protein_window']) <= limit)
    assert result == comparison["annotation_comparison"]
    assert metrics["comparison_status"] == metrics["focal_comparison_status"] == state
    if gene == "FHL2":
        assert coverage["coding_interval"][1] - coverage["coding_interval"][0] == 105  # >101-bp reads
    if gene == "TUBB2A":
        assert metrics["changed_residue_count"] == 0
