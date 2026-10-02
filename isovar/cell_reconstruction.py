"""Explicit reconstruction pools and conservative attribution of local RNA sequences.

Discovery retains the existing assembler's coverage floor. Attribution tests
original observations against all discovered candidates, including reference
allele competitors, without turning pooled support into per-cell confidence.
"""

from .base_quality import STATUSES
from .cell_evidence import CellUmiAlleles
from .cell_umi import UNTRUSTED_CELL_STATUSES
from .protein_hypotheses import _EvidenceSets, _event_id, _proteins
from .protein_sequence_creator import ProteinSequenceCreator
from .read_evidence import ReadEvidence
from .read_identity import observation_groups, read_sort_key, source_read_ids
from .rna_evidence import canonical_json, content_identifier
from .rna_candidate_support import (
    COMPARISON_POLICY, _attribute, _candidate, _comparison as _comparison, _protein_assignment,
    score_fragments, summarize_fragments, validate_comparison_settings,
)

SCHEMA = "isovar.cell_reconstruction.v1"
MODES = ("sample", "library", "cell", "group")


def _cell(row):
    if row["cell_barcode"] is None or row["status"] in UNTRUSTED_CELL_STATUSES:
        return None
    return dict(cell_id=content_identifier("cell", [row["scope"], row["cell_barcode"]]),
                scope=row["scope"], cell_barcode=row["cell_barcode"],
                library_scope_known=row["library_scope_known"])


def _read_cell(read, labels):
    cells = [_cell(labels.row(key)) for key in sorted(source_read_ids(read))]
    if not cells or any(cell is None for cell in cells):
        return None
    return cells[0] if len({cell["cell_id"] for cell in cells}) == 1 else None


def _pool(read, labels, mode, cell_groups):
    rows = [labels.row(key) for key in sorted(source_read_ids(read))]
    if not rows or any(row["scope"] is None for row in rows):
        return None
    if mode in ("cell", "group"):
        cell = _read_cell(read, labels)
        if cell is None:
            return None
        group = cell["cell_id"] if mode == "cell" else cell_groups.get(cell["cell_id"])
        if group is None:
            return None
        # Group names never erase the declared biological sample boundary.
        return dict(mode=mode, source=labels.source, sample_id=labels.sample_id,
                    header_sample=cell["scope"]["header_sample"], group=group)
    scopes = [row["scope"] for row in rows]
    if mode == "sample":
        scopes = [{key: scope[key] for key in ("source", "sample_id", "header_sample")} for scope in scopes]
    return dict(mode=mode, scope=scopes[0]) if len({canonical_json(s) for s in scopes}) == 1 else None


def _partition(reads, labels, mode, cell_groups):
    pools, excluded = {}, []
    for read in reads:
        description = _pool(read, labels, mode, cell_groups)
        if description is None:
            excluded.append(read)
            continue
        key = content_identifier("pool", description)
        pools.setdefault(key, (description, []))[1].append(read)
    return pools, excluded


def _cell_support(cells, key, ids_field, status_field, quality_enabled=False):
    result = {}
    if quality_enabled:
        result["base_quality"] = {status: {} for status in STATUSES}
    for category, status in (("unique", "unique_within_catalog"), ("ambiguous", "ambiguous")):
        cell_ids = sorted(cell_id for cell_id, cell in cells.items()
                          if any(o[status_field] == status and key in o[ids_field] for o in cell["observations"]))
        result[category + "_cell_ids"] = cell_ids
        result["num_" + category + "_cells"] = len(cell_ids)
        if quality_enabled:
            for quality_status in STATUSES:
                qualified = [cell_id for cell_id in cell_ids if any(
                    o[status_field] == status and key in o[ids_field]
                    and o["rna_support"]["base_quality"][quality_status]["reads"] > 0
                    for o in cells[cell_id]["observations"])]
                result["base_quality"][quality_status].update({
                    category + "_cell_ids": qualified, "num_" + category + "_cells": len(qualified)})
    return result


def _event(result, labels, evidence, mode, cell_groups, max_edits, min_overlap, transcript_id_whitelist,
           partial_read_support=False):
    read_evidence, variant = result.read_evidence, result.variant
    reads = sorted(set(read_evidence.ref_reads + read_evidence.alt_reads + read_evidence.other_reads), key=read_sort_key)
    event_id = _event_id(variant)
    pools, excluded = _partition(reads, labels, mode, cell_groups)
    settings = result.protein_sequence_settings
    if settings is None:
        raise ValueError("Cell reconstruction requires recorded protein_sequence_settings")
    # Attribution must see all alternatives, not only the ordinary top-protein cap.
    creator = ProteinSequenceCreator(**dict(settings, max_protein_sequences_per_variant=0))
    known_variants = next((p.known_variants for p in result.sorted_protein_sequences
                           if p.known_variants is not None), None)
    candidates, transcripts, pool_rows = {}, {}, []
    for pool_id, (description, pool_reads) in sorted(pools.items()):
        selected = set(pool_reads)
        subset = ReadEvidence(read_evidence.trimmed_base1_start, read_evidence.trimmed_ref, read_evidence.trimmed_alt,
                              *[[r for r in getattr(read_evidence, allele + "_reads") if r in selected]
                                for allele in ("ref", "alt", "other")])
        proteins = creator.sorted_protein_sequences_for_variant(variant, subset, transcript_id_whitelist)
        if known_variants is not None:
            proteins = [p.with_known_variants(known_variants) for p in proteins]
        exported = _proteins(proteins, event_id, evidence, labels)
        for protein, row in zip(proteins, exported):
            for translation in protein.translations:
                candidate, tx = _candidate(translation, event_id)
                key = candidate["candidate_id"]
                existing = candidates.setdefault(key, candidate)
                existing["pool_ids"].append(pool_id)
                existing["hypothesis_ids"].append(row["hypothesis_id"])
                transcripts[key] = tx
        pool_rows.append(dict(pool_id=pool_id, description=description,
                              rna_support=evidence.support(pool_reads, "reads_in_reconstruction_pool", labels),
                              protein_hypotheses=exported,
                              no_result_reason=None if proteins else "no_alt_reads" if not subset.alt_reads
                              else "no_protein_at_configured_coverage_and_reference_settings"))
    candidates = [candidates[key] for key in sorted(candidates)]
    for candidate in candidates:
        for field in ("pool_ids", "hypothesis_ids"):
            candidate[field] = sorted(set(candidate[field]))
    cells, unassigned = {}, []
    for observations in observation_groups(reads):
        cell_rows = [_read_cell(read, labels) for read in observations]
        cell = cell_rows[0]
        if cell is None or any(c != cell for c in cell_rows):
            unassigned.extend(observations)
            continue
        row = cells.setdefault(cell["cell_id"], dict(cell, observations=[]))
        eligible_pools = {pool_id for pool_id, (description, _) in pools.items()
                          if description.get("header_sample", description.get("scope", {}).get("header_sample"))
                          == cell["scope"]["header_sample"]}
        eligible_candidates = [c for c in candidates if eligible_pools.intersection(c["pool_ids"])]
        attribution = _attribute(observations, eligible_candidates, transcripts,
                                 read_evidence.trimmed_ref, max_edits, min_overlap)
        row["observations"].append(dict(attribution, **_protein_assignment(attribution, eligible_candidates),
                                        rna_support=evidence.support(
            observations, "original_observations_compared_with_candidate_catalog", labels)))
    # A synonymous RNA ambiguity may still identify one protein. A cell can
    # carry unique and ambiguous observations, or express multiple proteins.
    for candidate in candidates:
        candidate.update(_cell_support(cells, candidate["candidate_id"], "candidate_ids", "status",
                                       evidence.base_quality_policy is not None))
    protein_support = [dict(hypothesis_id=key, **_cell_support(cells, key, "hypothesis_ids", "protein_status",
                                                            evidence.base_quality_policy is not None))
                       for key in sorted({h for c in candidates for h in c["hypothesis_ids"]})]
    no_result = ("no_reads" if not reads else "no_eligible_pools" if not pools
                 else "no_reconstructed_candidates" if not candidates else None)
    event = dict(event_id=event_id, reconstruction_settings=creator.settings(), pools=pool_rows, candidates=candidates,
                no_result_reason=no_result, protein_attribution=protein_support,
                cells=[cells[key] for key in sorted(cells)],
                excluded_from_pool=evidence.support(excluded, "reads_without_requested_pool", labels),
                unresolved_cell_support=evidence.support(unassigned, "reads_without_resolved_cell", labels),
                input_support=evidence.support(reads, "all_input_observations", labels))

    if partial_read_support:
        def eligible(group):
            scopes = [_pool(read, labels, "sample", {}) for read in group]
            if any(scope is None for scope in scopes):
                return []
            samples = {scope["scope"]["header_sample"] for scope in scopes}
            if len(samples) != 1:
                return []
            sample, = samples
            eligible_pools = {key for key, (description, _) in pools.items()
                              if description.get("header_sample", description.get("scope", {}).get("header_sample"))
                              == sample}
            return [c for c in candidates if eligible_pools.intersection(c["pool_ids"])]

        def cell_id(group):
            rows = [_read_cell(read, labels) for read in group]
            return rows[0]["cell_id"] if rows and rows[0] is not None and all(r == rows[0] for r in rows) else None

        scored = score_fragments(reads, candidates, transcripts, read_evidence.trimmed_ref, evidence, labels,
                                 max_edits=max_edits, min_overlap=min_overlap, catalog_complete=True,
                                 eligible_candidates=eligible, cell_for_fragment=cell_id)
        event["partial_read_support"] = scored
        # The enclosing event already exports this exact candidate catalog.
        del scored["candidates"]
        scored["unresolved_cell_fragments"] = sum(f["cell_id"] is None for f in scored["fragments"])
        rna_scores = {r["candidate_id"]: r for r in scored["candidate_scores"]}
        protein_scores = {r["hypothesis_id"]: r for r in scored["protein_scores"]}
        for candidate in candidates:
            candidate["partial_read_support"] = rna_scores[candidate["candidate_id"]]
        for protein in protein_support:
            protein["partial_read_support"] = protein_scores[protein["hypothesis_id"]]
        for cell in event["cells"]:
            cell["partial_read_support"] = summarize_fragments(
                [f for f in scored["fragments"] if f["cell_id"] == cell["cell_id"]],
                rna_scores, protein_scores, evidence.base_quality_policy is not None)
    return event


def reconstruct_cell_groups(isovar_results, alignment_file, *, sample_id, source, mode="sample",
                            cell_groups=None, max_edits=1, min_overlap=10, read_collector=None,
                            transcript_id_whitelist=None, base_quality_policy=None, partial_read_support=False):
    """Reconstruct local protein candidates by pool and attribute original cell evidence.

    Parameters
    ----------
    isovar_results : iterable of IsovarResult
        Collected input evidence and recorded reconstruction settings. Pool
        discovery keeps all hypotheses regardless of the original protein cap.
    alignment_file : pysam.AlignmentFile
        Original alignments, used only to resolve cell/UMI labels.
    sample_id, source : str
        Explicit biological sample/timepoint and input namespaces.
    mode : {'sample', 'library', 'cell', 'group'}
        Discovery pool. Sample and group modes retain header SM boundaries;
        library and cell modes also retain the existing library/RG namespace.
    cell_groups : dict, optional
        Exported scoped cell ID to externally supplied group name. Required
        only for group mode; unmapped cells are excluded from discovery.
    max_edits : int
        Maximum anchored edit cost for compatibility (default 1). This is not
        a calibrated error probability. Every compatible alternative remains.
    min_overlap : int
        Minimum compared read bases including the allele (default 10).
    read_collector : ReadCollector, optional
        Collector used for the input, for metadata eligibility.
    transcript_id_whitelist : set of str, optional
        Restrict reconstruction to these reference transcripts.
    base_quality_policy : BaseQualityPolicy, optional
        Add reported focal-allele quality categories to support and cell counts.
        Does not filter reconstruction or attribution, or score the full ORF.
    partial_read_support : bool
        Add catalog-relative fractional fragment scores for RNA and protein
        hypotheses, including per-cell summaries. Default False.

    Returns
    -------
    dict
        Versioned JSON-ready export with pools, candidates, original observation
        assignments, scoped cell IDs and evidence sets. Unique means only within
        the discovered catalog; it does not establish full-ORF expression.
    """
    if mode not in MODES:
        raise ValueError("Unknown cell reconstruction mode: %s" % mode)
    validate_comparison_settings(max_edits, min_overlap)
    if mode == "group":
        if not isinstance(cell_groups, dict) or any(not isinstance(k, str) or not isinstance(v, str) or not v.strip()
                                                  for k, v in cell_groups.items()):
            raise ValueError("group mode requires a cell ID to nonempty group-name mapping")
    elif cell_groups is not None:
        raise ValueError("cell_groups applies only to group mode")
    evidence = _EvidenceSets([sample_id, source], base_quality_policy)
    resolver = CellUmiAlleles(alignment_file, sample_id=sample_id, source=source, read_collector=read_collector)
    events = []
    for result in isovar_results:
        reads = result.read_evidence
        identities = {key for read in reads.ref_reads + reads.alt_reads + reads.other_reads for key in source_read_ids(read)}
        labels = resolver.labels(result.variant, identities)
        events.append(_event(result, labels, evidence, mode, cell_groups, max_edits, min_overlap,
                             transcript_id_whitelist, partial_read_support))
    return dict(**({} if base_quality_policy is None else dict(base_quality_policy=base_quality_policy.description())),
                **({} if not partial_read_support else dict(partial_read_support=True)),
                schema=SCHEMA, mode=mode, sample_id=sample_id, source=source, max_edits=max_edits,
                min_overlap=min_overlap, cell_groups=cell_groups, score_policy=COMPARISON_POLICY,
                interpretation="Observed RNA compatibility within a discovered local candidate catalog; "
                               "not calibrated expression probabilities, independent molecules or cell prevalence.",
                events=events, evidence_sets={key: evidence.sets[key] for key in sorted(evidence.sets)})
