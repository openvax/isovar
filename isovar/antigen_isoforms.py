"""Variant-free, targeted RNA reconstruction with conservative attribution.

Uses the same observed placements, exact overlaps and original segment
identities as SV reconstruction. Reference models are competitors, never
sequence donors. Transcript abundance is estimated separately by Salmon.
"""

from collections import Counter, defaultdict
from contextlib import closing
from dataclasses import asdict
from hashlib import sha256
from pathlib import Path
import sqlite3
from tempfile import TemporaryDirectory

import pysam

from .antigen_reference import IsoformTarget
from .cell_umi import CellUmiEvidence
from .genetic_code import standard_genetic_code
from .read_collector import ReadCollector
from .read_identity import segment_identity
from .read_lineage import ReadLineage
from .read_metadata import record_evidence
from .sv_rna import (
    _GeometryCache, _Model, _ObservationIndex, _agrees, _build_segment, _extend_both,
    _offsets, _record_id, _runs,
)
from .sv_rna_dense import _DiskGroups, _spool


def _matches(sequence, positions, observation, *, contained=False):
    """All exact, coordinate-consistent overlaps with an original observation."""
    for _, offset in sorted(_offsets(_runs(positions), 0, observation.key, _runs(observation.positions))):
        if contained and (offset < 0 or offset + len(observation.sequence) > len(sequence)):
            continue
        lo, hi = max(0, offset), min(len(sequence), offset + len(observation.sequence))
        if lo < hi and _agrees(sequence, positions, observation, offset, lo, hi):
            yield lo, hi


def _contains(first, second):
    """Whether an entire (sequence, placements) is contained exactly in another."""
    sequence, positions = first
    other, mapping = second
    start = sequence.find(other)
    while start >= 0:
        if positions[start:start + len(other)] == mapping:
            return True
        start = sequence.find(other, start + 1)
    return False


def _coding(sequence, positions, models):
    readings = {}
    for model in models:
        reference = model.reference
        if reference.cds_start is None:
            continue
        # Transfer a frame through an observed, collinear complete codon.
        # Never choose a longest ORF or fill an uncovered start with reference.
        for q in range(len(sequence) - 2):
            offset = model.offsets.get(positions[q])
            if (offset is None or not reference.cds_start <= offset < reference.cds_end - 3
                    or (offset - reference.cds_start) % 3
                    or list(positions[q:q + 3]) != model.positions[offset:offset + 3]):
                continue
            complete_start = offset == reference.cds_start and sequence[q:q + 3] in standard_genetic_code.start_codons
            protein, stop = standard_genetic_code.translate(sequence[q:], first_codon_is_start=complete_start)
            if not protein:
                break
            key = (q, protein, complete_start, stop)
            row = readings.setdefault(key, dict(
                query_start=q, query_end=q + 3 * (len(protein) + stop), amino_acids=protein,
                complete_start=complete_start, observed_stop=stop,
                status="complete_CDS" if complete_start and stop else "partial_CDS", transcript_ids=[]))
            row["transcript_ids"].append(reference.transcript_id)
            break
    return sorted(readings.values(), key=lambda row: (row["query_start"], row["amino_acids"]))


def _candidate(sequence, positions, targets, models, annotated=()):
    positions = tuple(positions)
    genes = sorted({t.gene_id for t in targets if any(
        p is not None and p[0] == t.reference.contig and p[2] == t.reference.strand
        and t.reference.exons[0][0] <= p[1] < t.reference.exons[-1][1] for p in positions)})
    if annotated:
        genes = sorted({t.gene_id for t in targets if t.reference.transcript_id in annotated})
    key = repr((sequence, positions)).encode()
    reaches_ends = any(positions[0] == model.positions[0] and positions[-1] == model.positions[-1]
                       for model in models)
    return dict(id="isovar_isoform_" + sha256(key).hexdigest()[:24], sequence=sequence,
                placements=positions, gene_ids=genes, annotated_transcript_ids=sorted(annotated),
                origin="reference_model" if annotated else "observed_assembly",
                reaches_annotated_transcript_ends=reaches_ends,
                coding=_coding(sequence, positions, models))


def _relevant(read, targets):
    return not read.is_unmapped and any(
        read.reference_name == target.reference.contig
        and read.reference_start < target.reference.exons[-1][1]
        and read.reference_end > target.reference.exons[0][0] for target in targets)


def _target_orientation(observation, targets):
    return any(p is not None and p[0] == target.reference.contig and p[2] == target.reference.strand
               and target.reference.exons[0][0] <= p[1] < target.reference.exons[-1][1]
               for target in targets for p in observation.positions)


def _observations(template, collector, reasons, targets, cache):
    grouped = defaultdict(list)
    for read in template:
        grouped[segment_identity(read)].append(read)
    for identity, records in sorted(grouped.items()):
        if any(_relevant(read, targets) for read in records):
            # The placement interner is per template: full-input support does
            # not retain every background genomic base in an unbounded cache.
            yield identity, [o for o in _build_segment(identity, records, collector, reasons, {}, cache)
                             if _target_orientation(o, targets)]


def _discover(groups, targets, collector, max_discovery_reads, max_paths, min_overlap):
    chosen, reasons, cache = {}, Counter(), _GeometryCache()
    considered, selected = 0, 0
    for template in groups.templates():
        for _, observations in _observations(template, collector, reasons, targets, cache):
            considered += 1
            if selected >= max_discovery_reads:
                continue
            selected += 1
            for o in observations:
                chosen.setdefault((o.sequence, o.positions), o)
    index = _ObservationIndex(list(chosen.values()))
    paths, notes = {}, set()
    for seed in sorted(chosen.values(), key=lambda o: (-len(o.sequence), o.sequence, repr(o.positions))):
        key = seed.sequence, tuple(seed.positions)
        if any(_contains(path, key) for path in paths):
            continue
        if len(paths) >= max_paths:
            notes.add("path_limit")
            break
        for sequence, positions, _, _ in _extend_both(
                *key, index, min_overlap, (0, 0, 0), max_paths - len(paths), [], notes):
            paths[sequence, tuple(positions)] = None
            if len(paths) >= max_paths:
                notes.add("path_limit")
                break
    maximal = [p for p in paths if not any(len(q[0]) > len(p[0]) and _contains(q, p) for q in paths)]
    return maximal, dict(target_segments=considered, selected_segments=selected,
                         distinct_observations=len(chosen), max_discovery_reads=max_discovery_reads,
                         max_paths=max_paths, min_overlap=min_overlap,
                         truncated=considered > selected or bool(notes), limitations=sorted(notes),
                         observation_rejections=dict(reasons))


def _support(groups, candidates, targets, collector, header, sample_id, source):
    coverage = {c["id"]: [0] * len(c["sequence"]) for c in candidates}
    supporters, whole, assignments = defaultdict(set), defaultdict(set), Counter()
    whole_cds = defaultdict(set)
    witnesses, reasons, cache = [], Counter(), _GeometryCache()
    lineage = ReadLineage(groups, header)
    labels = CellUmiEvidence(groups, header, sample_id, source, template_loader=groups.template)
    target_fragments, unassigned = 0, 0
    for template in groups.templates():
        observations = list(_observations(template, collector, reasons, targets, cache))
        if not observations:
            continue
        target_fragments += 1
        fragment_matches, segment_rows = None, []
        outside = any(not _relevant(read, targets) for read in template)
        # NH advertises unreported competing placements; a unique retained
        # alignment cannot establish unique target attribution in that case.
        outside |= any(read.has_tag("NH") and read.get_tag("NH") > len(
            [r for r in template if segment_identity(r) == segment_identity(read) and not r.is_supplementary])
            for read in template)
        outside |= any(read.mapping_quality == 0 or (read.mapping_quality == 255 and not read.has_tag("NH"))
                       for read in template)
        for identity, observed in observations:
            compatible, supporting = set(), set()
            for candidate in candidates:
                cid = candidate["id"]
                covered, full = set(), False
                for o in observed:
                    for lo, hi in _matches(candidate["sequence"], candidate["placements"], o):
                        covered.update(range(lo, hi))
                        full |= lo == 0 and hi == len(candidate["sequence"])
                        for i, coding in enumerate(candidate["coding"]):
                            if coding["status"] == "complete_CDS" and lo <= coding["query_start"] and hi >= coding["query_end"]:
                                whole_cds[cid, i].add(identity)
                    if any(_matches(candidate["sequence"], candidate["placements"], o, contained=True)):
                        compatible.add(cid)
                if covered:
                    supporting.add(cid)
                    supporters[cid].add(identity)
                    for offset in covered:
                        coverage[cid][offset] += 1
                if full:
                    whole[cid].add(identity)
            fragment_matches = compatible if fragment_matches is None else fragment_matches & compatible
            segment_rows.append(dict(identity=list(identity), compatible_isoforms=sorted(compatible),
                                     overlapping_isoforms=sorted(supporting),
                                     lineage=lineage.segment(identity), cell_umi=labels.segment(identity)))
        ids = tuple(sorted(fragment_matches or ()))
        assignments[ids, outside] += 1
        unassigned += not ids
        witnesses.append(dict(fragment=list(observations[0][0][:2]), compatible_isoforms=list(ids),
                              competing_placement_outside_targets=outside, segments=segment_rows,
                              records=[dict(id=_record_id(read), **record_evidence(read)) for read in template]))
        # Witnesses hold the consulted metadata. Caches need not grow with the
        # total file; ancestry can still consult the disk-backed source index.
        lineage.cache.clear()
        labels.cache.clear()
    for candidate in candidates:
        cid = candidate["id"]
        candidate.update(coverage=coverage[cid], supporting_segments=len(supporters[cid]),
                         supporting_fragments=len({key[:2] for key in supporters[cid]}),
                         whole_path_segments=[list(key) for key in sorted(whole[cid])],
                         unique_fragments=assignments[((cid,), False)],
                         compatible_fragments=sum(n for (ids, _), n in assignments.items() if cid in ids))
        candidate["sequence_status"] = (
            "whole_path_observed" if whole[cid] else "fully_covered_hypothesis" if all(coverage[cid])
            else "partially_covered_reference" if any(coverage[cid]) else "not_observed")
        for i, coding in enumerate(candidate["coding"]):
            coding["whole_cds_segments"] = [list(key) for key in sorted(whole_cds[cid, i])]
            coding["all_bases_observed"] = all(coverage[cid][coding["query_start"]:coding["query_end"]])
        candidate["rna_support"] = labels.support(supporters[cid])
        candidate["signal_support"] = lineage.support(supporters[cid])
        labels.cache.clear()
        lineage.cache.clear()
    classes = [dict(isoform_ids=list(ids), outside_targets=outside, fragments=n)
               for (ids, outside), n in sorted(assignments.items())]
    return dict(scope="complete_declared_input", target_fragments=target_fragments,
                unassigned_target_fragments=unassigned, equivalence_classes=classes,
                observation_rejections=dict(reasons), witnesses=witnesses)


def _protein_groups(candidates, support):
    groups = defaultdict(lambda: dict(isoform_ids=[], fully_covered_isoform_ids=[], whole_cds_segments=set()))
    for row in candidates:
        products = {c["amino_acids"] for c in row["coding"] if c["status"] == "complete_CDS"}
        for protein in sorted(products):
            group = groups[protein]
            group["isoform_ids"].append(row["id"])
            matching = [c for c in row["coding"] if c["amino_acids"] == protein and c["status"] == "complete_CDS"]
            if any(c["all_bases_observed"] for c in matching):
                group["fully_covered_isoform_ids"].append(row["id"])
            group["whole_cds_segments"].update(tuple(key) for c in matching for key in c["whole_cds_segments"])
    output = []
    for protein, group in sorted(groups.items()):
        ids = set(group["isoform_ids"])
        group["whole_cds_segments"] = [list(key) for key in sorted(group["whole_cds_segments"])]
        group["compatible_fragments"] = sum(c["fragments"] for c in support["equivalence_classes"] if ids.intersection(c["isoform_ids"]))
        output.append(dict(amino_acids=protein, **group))
    return output


def reconstruct_antigen_isoforms(bam, *, targets, sample_id, source, reference_identity=None,
                                max_discovery_reads=10000, max_paths=256, min_overlap=30,
                                scratch_dir=None):
    """Reconstruct targeted RNA without a nominated mutation.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Complete declared input; BAM, SAM and CRAM use the existing disk spool.
    targets : iterable of IsoformTarget
        All annotated competitors of the explicitly selected genes. Include
        related loci when shared reads could otherwise be misattributed.
    sample_id, source : str
        Explicit namespaces retained alongside original RG/QNAME identities.
    reference_identity : dict, optional
        Alignment/annotation identity and file checksums supplied by the caller.
    max_discovery_reads, max_paths, min_overlap : int
        Discovery resource limits and minimum exact overlap. All eligible
        records are subsequently checked for support, including omitted reads.

    Returns
    -------
    dict
        Versioned sequence/evidence result. Compatible fragments are not TPM or
        independent molecules. Noisy reads may support structure without an
        exact sequence match; this initial workflow leaves them unassigned.
    """
    targets = tuple(targets)
    if not targets or not all(isinstance(t, IsoformTarget) for t in targets):
        raise ValueError("Provide at least one IsoformTarget")
    targets = tuple(sorted(targets, key=lambda t: t.reference.transcript_id))
    if not sample_id or not source:
        raise ValueError("Explicit sample_id and source are required")
    samples = {row["SM"] for row in bam.header.to_dict().get("RG", []) if row.get("SM")}
    if len(samples) > 1:
        raise ValueError("The BAM declares multiple samples; reconstruct each sample separately")
    if len({t.reference.reference_name for t in targets}) != 1:
        raise ValueError("Targets must share one alignment assembly")
    for name, value in (("max_discovery_reads", max_discovery_reads), ("max_paths", max_paths), ("min_overlap", min_overlap)):
        if type(value) is not int or value < 1:
            raise ValueError("%s must be a positive integer" % name)
    names = [t.reference.transcript_id for t in targets]
    if len(set(names)) != len(names):
        raise ValueError("Reference transcript IDs must be unique")
    lengths = dict(zip(bam.references, bam.lengths))
    for target in targets:
        r = target.reference
        if r.contig not in lengths or r.exons[-1][1] > lengths[r.contig]:
            raise ValueError("Target transcript lies outside the BAM reference: %s" % r.transcript_id)
    collector = ReadCollector(use_secondary_alignments=True, min_mapping_quality=0,
                              use_soft_clipped_bases=False, use_reads_without_base_qualities=True)
    models = [_Model(t.reference) for t in targets]
    with TemporaryDirectory(dir=scratch_dir, prefix="isovar-isoforms-") as directory, closing(
            sqlite3.connect(str(Path(directory) / "reads.sqlite"))) as db:
        scanned, excluded, unique, segments, path = _spool(bam, db, collector, directory)
        with pysam.AlignmentFile(path, "rb") as reader:
            groups = _DiskGroups(db, reader, bam.header, segments)
            paths, discovery = _discover(groups, targets, collector, max_discovery_reads, max_paths, min_overlap)
            candidates = {}
            for target, model in zip(targets, models):
                r = target.reference
                key = r.sequence, tuple(model.positions)
                if key in candidates:
                    candidates[key]["annotated_transcript_ids"].append(r.transcript_id)
                    candidates[key]["gene_ids"] = sorted(set(candidates[key]["gene_ids"]) | {target.gene_id})
                else:
                    candidates[key] = _candidate(*key, targets, models, [r.transcript_id])
            partial = []
            for sequence, positions in paths:
                key = sequence, positions
                if key in candidates:
                    continue
                row = _candidate(*key, targets, models)
                containers = [c["id"] for k, c in candidates.items() if _contains(k, key)]
                if containers or not row["reaches_annotated_transcript_ends"]:
                    row["compatible_models"] = sorted(containers)
                    partial.append(row)
                else:
                    candidates[key] = row
            candidates = sorted(candidates.values(), key=lambda c: c["id"])
            support = _support(groups, candidates + partial, targets, collector, bam.header.to_dict(), sample_id, source)
            full_ids = {c["id"] for c in candidates}
            support["fragments_without_full_transcript_assignment"] = sum(
                row["fragments"] for row in support["equivalence_classes"] if not full_ids.intersection(row["isoform_ids"]))
    for row in candidates + partial:
        # Compact placed runs leave insertions explicitly unplaced; query
        # sequence length remains available to reconstruct the placement map.
        row["mapping"] = [list(run) for run in _runs(row.pop("placements"))]
    return dict(schema="isovar.antigen_isoforms.v1", sample_id=sample_id, source=source,
                reference=reference_identity or {}, targets=[asdict(t) for t in targets],
                input=dict(records=scanned, unique_records=unique, eligible_segments=segments,
                           excluded_records=excluded, header=bam.header.to_dict()),
                discovery=discovery, isoforms=candidates, partial_sequences=partial, support=support,
                protein_groups=_protein_groups(candidates + partial, support), quantification=None,
                limitations=["exact_overlap_reconstruction", "RNA_does_not_establish_translation",
                             "fragment_counts_are_not_molecule_counts",
                             "unreported_alignments_not_recovered",
                             "overlap_assembly_does_not_establish_end_to_end_phase"])
