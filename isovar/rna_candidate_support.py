"""Shared local RNA compatibility and catalog-relative fragment support.

Compatibility is an explicit edit/path rule, not a probability of expression.
Every scored fragment has a conserved unit budget across compatible hypotheses.
"""

from fractions import Fraction
from numbers import Integral

import edlib

from .dna import reverse_complement_dna
from .read_identity import fragment_ids, observation_groups, read_sort_key, source_read_ids
from .rna_evidence import canonical_json, content_identifier
from .transcript_compatibility import annotate_reads_with_transcript_compatibility

COMPARISON_POLICY = "anchored_edit_compatibility.v2"

def _candidate(translation, event_id):
    orf = translation.variant_orf
    sequence = orf.cdna_sequence
    start, end = orf.variant_cdna_interval_start, orf.variant_cdna_interval_end
    if translation.reference_context.strand == "-":
        sequence = reverse_complement_dna(sequence)
        start, end = len(sequence) - end, len(sequence) - start
    transcripts = tuple(sorted(translation.reference_context.transcripts, key=lambda t: t.id))
    identity = [event_id, sequence, start, end, sorted(t.id for t in transcripts)]
    return dict(candidate_id=content_identifier("cell_candidate", identity),
                sequence=sequence, variant_interval=[start, end],
                transcript_ids=[t.id for t in transcripts], pool_ids=[], hypothesis_ids=[]), transcripts


def _distance(a, b):
    if not a or not b:
        return max(len(a), len(b))
    return edlib.align(a, b, mode="NW", task="distance")["editDistance"]


def _flank_alignment(read, candidate):
    """Align outward from the allele until either flank ends.

    Both starts are fixed; only distal overhangs are free. Prefix alignment in
    both directions considers either sequence ending first, even at equal
    lengths. Retain all optimal (read bases, candidate bases) endpoint pairs
    so tied alignments cannot invent certain coverage.
    """
    best, endpoints = -1, set()
    directions = sorted(((read, candidate, False), (candidate, read, True)), key=lambda row: len(row[0]))
    for query, target, swapped in directions:
        result = edlib.align(query, target, mode="SHW", task="distance", k=best)
        cost = result["editDistance"]
        if cost < 0 or (best >= 0 and cost > best):
            continue
        if best < 0 or cost < best:
            best, endpoints = cost, set()
        for _, end in result["locations"]:
            pair = (len(query), end + 1)
            endpoints.add(pair[::-1] if swapped else pair)
    return dict(edit_distance=best, endpoints=[list(pair) for pair in sorted(endpoints)])


def _comparison(read, candidate, ref, transcripts, max_edits, min_overlap):
    sequence = candidate["sequence"]
    start, end = candidate["variant_interval"]
    left = _flank_alignment(read.prefix[::-1], sequence[:start][::-1])
    right = _flank_alignment(read.suffix, sequence[end:])
    before, after = [min(q for q, _ in flank["endpoints"]) for flank in (left, right)]
    candidate_before, candidate_after = [min(c for _, c in flank["endpoints"]) for flank in (left, right)]
    # Coverage uses the intersection of tied endpoint spans. Path and alphabet
    # checks use their union: no tie may hide an incompatible path/unknown base.
    query_before, query_after = [max(q for q, _ in flank["endpoints"]) for flank in (left, right)]
    target_before, target_after = [max(c for _, c in flank["endpoints"]) for flank in (left, right)]
    overlap = before + len(read.allele) + after
    compared = (read.prefix[len(read.prefix) - query_before:] + read.allele + read.suffix[:query_after]
                + sequence[start - target_before:end + target_after] + ref)
    if overlap < min_overlap or set(compared) - set("ACGT"):
        return None
    classified, = annotate_reads_with_transcript_compatibility(
        [read], transcripts, max_prefix_size=query_before, max_suffix_size=query_after)
    if not classified.compatible_transcript_ids:
        return None
    flanks = left["edit_distance"] + right["edit_distance"]
    alt_cost = flanks + _distance(read.allele, sequence[start:end])
    ref_cost = flanks + _distance(read.allele, ref)
    return dict(candidate_id=candidate["candidate_id"], edit_distance=alt_cost,
                reference_edit_distance=ref_cost, within_tolerance=alt_cost <= max_edits,
                candidate_interval=[start - candidate_before, end + candidate_after],
                spans_candidate=candidate_before == start and candidate_after == len(sequence) - end,
                compared_read_bases=overlap, flank_alignments=dict(prefix=left, suffix=right),
                reference_preferred=ref_cost <= alt_cost)


def _attribute(reads, candidates, transcripts, ref, max_edits, min_overlap):
    comparisons = []
    for read in reads:
        comparisons.append([comparison for candidate in candidates
                            if (comparison := _comparison(read, candidate, ref, transcripts[candidate["candidate_id"]],
                                                          max_edits, min_overlap)) is not None])
    # Alternative placements of one observation must all agree before it is
    # called unique. Keep the union visible if any placement disagrees.
    supports = [{c["candidate_id"] for c in rows if c["within_tolerance"] and not c["reference_preferred"]}
                for rows in comparisons]
    supported = set.union(*supports) if supports else set()
    if not candidates:
        status = "no_candidate"
    elif not any(comparisons):
        status = "insufficient_coverage_or_incompatible_path"
    elif not supported:
        status = ("reference_or_other_allele" if any(c["reference_preferred"] for rows in comparisons for c in rows)
                  else "outside_tolerance")
    elif len(supported) > 1 or any(s != supported for s in supports):
        status = "ambiguous"
    else:
        status = "unique_within_catalog"
    return dict(status=status, candidate_ids=sorted(supported), comparisons=comparisons,
                placement_ambiguous=any(s != supported for s in supports))


def _protein_assignment(attribution, candidates):
    ids = {hypothesis for candidate in candidates if candidate["candidate_id"] in attribution["candidate_ids"]
           for hypothesis in candidate["hypothesis_ids"]}
    status = attribution["status"]
    if ids:
        status = "unique_within_catalog" if len(ids) == 1 and not attribution["placement_ambiguous"] else "ambiguous"
    return dict(hypothesis_ids=sorted(ids), protein_status=status)


def validate_comparison_settings(max_edits, min_overlap):
    for name, value, minimum in (("max_edits", max_edits, 0), ("min_overlap", min_overlap, 1)):
        if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
            raise ValueError("%s must be an integer >= %d" % (name, minimum))


def candidate_catalog(proteins, protein_rows, event_id):
    """The nucleotide hypotheses behind proteins, in original SAM orientation."""
    candidates, transcripts = {}, {}
    for protein, row in zip(proteins, protein_rows):
        for translation in protein.translations:
            candidate, tx = _candidate(translation, event_id)
            key = candidate["candidate_id"]
            existing = candidates.setdefault(key, candidate)
            existing["hypothesis_ids"].append(row["hypothesis_id"])
            transcripts[key] = tx
    for candidate in candidates.values():
        candidate["hypothesis_ids"] = sorted(set(candidate["hypothesis_ids"]))
    return [candidates[key] for key in sorted(candidates)], transcripts


def _fragments(reads):
    groups = {}
    for read in sorted(reads, key=read_sort_key):
        keys = fragment_ids([read])
        if len(keys) != 1:
            raise ValueError("A read observation must describe one fragment")
        key, = keys
        # Keep legacy QNAME-only identities separate from collected RG/QNAMEs.
        identity = ("sam" if source_read_ids(read) else "legacy", key)
        groups.setdefault(identity, []).append(read)
    return [groups[key] for key in sorted(groups, key=canonical_json)]


def _fragment_assignment(reads, candidates, transcripts, ref, max_edits, min_overlap):
    # A retained view of an exact source alignment already covered by a merged
    # consensus is not an alternative placement. Keep conflicting sequences or
    # different placements; their uncertainty must remain visible.
    reads = sorted(set(reads), key=read_sort_key)
    reads = [read for read in reads if not any(
        set(getattr(read, "source_alignments", ())) < set(getattr(other, "source_alignments", ()))
        and source_read_ids(read) and read.allele == other.allele
        and other.prefix.endswith(read.prefix) and other.suffix.startswith(read.suffix)
        for other in reads)]
    observations = [_attribute(group, candidates, transcripts, ref, max_edits, min_overlap)
                    for group in observation_groups(reads)]
    sets = [set(o["candidate_ids"]) for o in observations]
    joint = set.intersection(*sets) if sets else set()
    if not candidates:
        status = "no_candidate"
    elif any(o["placement_ambiguous"] for o in observations):
        status = "alternative_placements_disagree"
    elif any(not s for s in sets):
        status = "observation_without_compatible_candidate"
    elif not joint:
        status = "conflicting_observations"
    else:
        status = "scored"
    ids = sorted(joint) if status == "scored" else []
    proteins = sorted({h for c in candidates if c["candidate_id"] in ids for h in c["hypothesis_ids"]})
    return dict(status=status, candidate_ids=ids, hypothesis_ids=proteins,
                candidate_weight=1 / len(ids) if ids else 0.0,
                protein_weight=1 / len(proteins) if proteins else 0.0,
                observations=observations)


def _quality_category(support):
    quality = support.get("base_quality")
    if quality is None:
        return None
    if quality["failed"]["reads"]:
        return "failed"
    if quality["unassessed"]["reads"] or not quality["passed"]["reads"]:
        return "unassessed"
    return "passed"


def summarize_fragments(fragments, candidate_ids, hypothesis_ids, quality_enabled=False):
    """Summarize the same fragment allocations at RNA and protein resolution.

    Fractions are accumulated exactly so input order cannot change scores or
    dense ranks. Projection to distinct protein IDs occurs before allocation.
    """
    fragments = list(fragments)

    def scores(keys, field, id_field):
        rows, totals = [], {}
        for key in sorted(set(keys)):
            selected = [f for f in fragments if key in f[field]]
            total = sum((Fraction(1, len(f[field])) for f in selected), Fraction())
            row = {id_field: key, "fractional_fragments": float(total),
                   "unique_fragments": sum(len(f[field]) == 1 for f in selected),
                   "ambiguous_fragments": sum(len(f[field]) > 1 for f in selected)}
            if quality_enabled:
                row["base_quality_fractional_fragments"] = {
                    status: float(sum((Fraction(1, len(f[field])) for f in selected
                                       if f["base_quality_status"] == status), Fraction()))
                    for status in ("passed", "failed", "unassessed")}
            rows.append(row)
            totals[key] = total
        ranks = {value: i for i, value in enumerate(sorted(set(totals.values()), reverse=True), 1)}
        for row in rows:
            row["rank"] = ranks[totals[row[id_field]]]
        return rows

    scored = sum(f["status"] == "scored" for f in fragments)
    return dict(input_fragments=len(fragments), scored_fragments=scored,
                unscored_fragments=len(fragments) - scored,
                candidate_scores=scores(candidate_ids, "candidate_ids", "candidate_id"),
                protein_scores=scores(hypothesis_ids, "hypothesis_ids", "hypothesis_id"))


def score_fragments(reads, candidates, transcripts, ref, evidence, labels=None, *,
                    max_edits=1, min_overlap=10, catalog_complete=None, eligible_candidates=None,
                    cell_for_fragment=None):
    """Score original partial reads without modifying assembly or raw counts.

    Parameters
    ----------
    reads : iterable of AlleleRead
        Original observations at the event (ref, alt and other).
    candidates, transcripts : list, dict
        Candidate catalog and corresponding transcript references.
    ref : str
        Trimmed focal reference allele, for reference-competitor comparisons.
    evidence : protein_hypotheses._EvidenceSets
        Scoped support/evidence-set collector, optionally with quality policy.
    labels : CellUmiEvidence, optional
        Original cell/UMI labels.
    max_edits, min_overlap : int
        Anchored compatibility settings, shared with cell attribution.
    catalog_complete : bool or None
        Whether all discovered hypotheses were retained; never asserts all
        possible expressed sequences have been discovered.
    eligible_candidates, cell_for_fragment : callable, optional
        Restrict comparison to eligible biological pools and resolve cell IDs.

    Returns
    -------
    dict
        Policy, catalog, fragment assignments and independently normalized RNA
        and protein scores. These are not calibrated abundance probabilities.
    """
    validate_comparison_settings(max_edits, min_overlap)
    candidates = sorted(candidates, key=lambda c: c["candidate_id"])
    fragments = []
    for group in _fragments(reads):
        selected = candidates if eligible_candidates is None else sorted(
            eligible_candidates(group), key=lambda c: c["candidate_id"])
        row = _fragment_assignment(group, selected, transcripts, ref, max_edits, min_overlap)
        row["rna_support"] = evidence.support(group, "original_fragment_compared_with_candidate_catalog", labels)
        identities = {key[:2] for read in group for key in source_read_ids(read)}
        row["fragment_id"] = (content_identifier("fragment", [evidence.scope, next(iter(identities))])
                              if len(identities) == 1 else None)
        if evidence.base_quality_policy is not None:
            row["base_quality_status"] = _quality_category(row["rna_support"])
        if cell_for_fragment is not None:
            row["cell_id"] = cell_for_fragment(group)
        fragments.append(row)
    summary = summarize_fragments(fragments, [c["candidate_id"] for c in candidates],
                                  [h for c in candidates for h in c["hypothesis_ids"]],
                                  evidence.base_quality_policy is not None)
    return dict(policy=dict(name="fractional_compatible_fragments.v1", unit="fragment",
                            comparison=COMPARISON_POLICY,
                            max_edits=max_edits, min_overlap=min_overlap,
                            allocation="equal_share_of_compatible_hypotheses",
                            length_weighted=False, calibrated_probability=False,
                            placement_disagreement="unscored", missing_comparison="unscored"),
                catalog_complete=catalog_complete, candidates=candidates,
                fragments=fragments, **summary)
