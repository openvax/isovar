"""Reconstruct RNA paths around one nominated structural variant.

Starting from an oriented DNA adjacency and annotated transcript models, this
module retrieves regional alignments, seeds candidate paths at unannotated RNA
junctions, extends them through exact, coordinate-anchored overlaps and
transfers reading frames from exact collinear CDS anchors.

Results are exploratory candidates, not the validated supplied-fusion schema.
RNA sequence support, reading frame and linkage to the nominated event are
independent evidence axes. A reference base never fills an unobserved RNA
base, an SA tag never substitutes for its record, and a branch that reads
cannot resolve is kept rather than decided by choosing one read.

Coordinates are 0-based. Each path base has one ``(contig, position, strand)``
placement in the path's 5'->3' orientation, or ``None`` when no single observed
placement exists (insertions, soft clips, junction homology).
"""

from collections import defaultdict
from dataclasses import asdict
from hashlib import sha256
from typing import TYPE_CHECKING

from .read_identity import sam_record_id
from .rna_observation import ObservationIndex, alignment_gaps
from .rna_observation import RnaObservation as RnaObservation
from .rna_path_assembly import extend_path_both_directions, supported_alternatives
from .rna_path_geometry import interval_contains, merge_intervals
from .sv_rna_collection import collect_sv_records
from .sv_rna_event import (
    OrientedAdjacency,
    SvObservationStore,
    contains_seed,
    sv_path_seeds,
)
from .sv_rna_result import build_sv_path_result

if TYPE_CHECKING:
    from .rna_library import RnaLibraryIndex

from .default_parameters import (
    FUSION_PEPTIDE_LENGTHS,
    SV_ANNOTATED_JUNCTION_TOLERANCE,
    SV_ASSEMBLE,
    SV_BREAKPOINT_WINDOW,
    SV_INCLUSION_MAPQ_255_IS_UNIQUE,
    SV_INCLUSION_MIN_BASE_QUALITY,
    SV_INCLUSION_MIN_MAPPING_QUALITY,
    SV_INCLUSION_MIN_SPLICE_ANCHOR_BASES,
    SV_MAX_BREAKPOINT_SHIFT,
    SV_MAX_EXTENSION_SEGMENTS,
    SV_MAX_ORF_CANDIDATES,
    SV_MAX_PATHS,
    SV_MAX_QUERIES,
    SV_MAX_RECORDS,
    SV_MIN_ALTERNATIVE_FRACTION,
    SV_MIN_ALTERNATIVE_FRAGMENTS,
    SV_MIN_ANCHOR_BASES,
    SV_MIN_LOCAL_VARIANT_FRACTION,
    SV_MIN_ORF_AMINO_ACIDS,
    SV_MIN_OVERLAP,
)
from .fusion import FusionBreakpoint, FusionReference
from .read_collector import ReadCollector
from .read_metadata import record_evidence
from .rna_path_geometry import reverse_position
from .sv_rna_orfs import exploratory_orfs
from .sv_rna_relations import EVENT_LINKED_RELATIONS, LINKAGE_STATUSES
from .transcript_path_frames import TranscriptPathModel as TranscriptPathModel
from .transcript_path_frames import sv_frame_evidence as sv_frame_evidence


def reconstruct_sv_rna(bam, *, event_id, reference_name, donor, acceptor, regions, references,
                       sample_id, source, event_provenance, read_collector=None,
                       min_anchor_bases=SV_MIN_ANCHOR_BASES, min_overlap=SV_MIN_OVERLAP,
                       min_alternative_fragments=SV_MIN_ALTERNATIVE_FRAGMENTS,
                       min_alternative_fraction=SV_MIN_ALTERNATIVE_FRACTION,
                       min_local_variant_fraction=SV_MIN_LOCAL_VARIANT_FRACTION, max_records=SV_MAX_RECORDS,
                       max_queries=SV_MAX_QUERIES, max_paths=SV_MAX_PATHS,
                       max_extension_segments=SV_MAX_EXTENSION_SEGMENTS, assemble=SV_ASSEMBLE,
                       include_regional_candidates=True, exon_anchored_clips=False,
                       breakpoint_window=SV_BREAKPOINT_WINDOW, max_breakpoint_shift=SV_MAX_BREAKPOINT_SHIFT,
                       annotated_junction_tolerance=SV_ANNOTATED_JUNCTION_TOLERANCE,
                       peptide_lengths=FUSION_PEPTIDE_LENGTHS,
                       min_orf_amino_acids=SV_MIN_ORF_AMINO_ACIDS, max_orf_candidates=SV_MAX_ORF_CANDIDATES,
                       inclusion_min_splice_anchor_bases=SV_INCLUSION_MIN_SPLICE_ANCHOR_BASES,
                       inclusion_min_base_quality=SV_INCLUSION_MIN_BASE_QUALITY,
                       inclusion_min_mapping_quality=SV_INCLUSION_MIN_MAPPING_QUALITY,
                       inclusion_mapq_255_is_unique=SV_INCLUSION_MAPQ_255_IS_UNIQUE,
                       dense_support=False, scratch_dir=None, max_evidence_groups=None,
                       library_index: "RnaLibraryIndex | None" = None, read_filter_identity: str | None = None):
    """Reconstruct RNA paths for one nominated, oriented DNA adjacency.

    Paths start at unannotated RNA junctions (and, when soft-clipped bases are
    used, unaligned sequence at a breakpoint) and extend through exact overlaps
    anchored by shared placements. With ``assemble=False`` only observations
    containing the seed junction contribute, so every base is observed in a
    read which also spans that junction. Frames are transferred from exact,
    collinear annotated CDS anchors into downstream observed sequence.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Indexed alignments from one sample/library source. Missing QUAL is kept.
    donor, acceptor : FusionBreakpoint
        Nominated oriented adjacency, retained-partner boundaries as in
        ``FusionTranscript``. It need not be an RNA junction.
    regions : iterable of (str, int, int)
        Additional genomic search intervals. Annotated exons and
        ``breakpoint_window`` bases around each breakpoint are always searched.
    references : iterable of FusionReference
        All candidate transcript models. Peptide novelty is relative to them.
        An explicitly empty iterable permits intergenic exploratory ORFs;
        annotated splice, start and frame assessment is then unavailable and
        reported in ``limitations``. No reference sequence is substituted.
    sample_id, source : str
        Sample and alignment-source identity.
    event_provenance : dict
        Caller-supplied DNA event evidence, retained verbatim.
    read_collector : ReadCollector, optional
        Read filters, soft-clip use and read-end handling.
    min_anchor_bases : int
        Exact collinear transcript match required to transfer a frame.
    min_overlap : int
        Overlap required to extend a path with another observation.
    min_alternative_fragments : int
        Fragments needed for an unattributed junction to seed a path, and for
        an extension branch to outweigh (and prune) a weaker alternative.
    min_alternative_fraction : float
        An alternative extension branch, or a sequence variant of an already
        seeded junction, is also pruned when its fragments are below this
        fraction of the best-supported alternative. Pruning is reported.
    min_local_variant_fraction : float
        A base or small-indel alternative (placed within 20 bases of the best,
        or a variant of a seeded junction's bases) also needs this fraction of
        the best's support, so systematic long-read errors do not fork paths.
    max_records, max_queries, max_paths : int
        Resource limits; reaching one is reported in ``limitations``.
    max_evidence_groups : int or None
        Optional dense-discovery limit on distinct exact structural/allelic
        groups. The independent ``max_records`` budget still applies; all
        eligible original records remain available to the support recount.
    max_extension_segments : int
        Segments built to extend one path end; those reaching furthest are
        used. Reaching it is reported as ``extension_segment_limit``.
    assemble : bool
        Extend through overlapping observations which do not span the seed.
        Also retain seed-spanning reconstructions not contained in a regional
        assembly, subject to the same path limit.
    include_regional_candidates : bool
        Also seed paths at novel regional junctions unrelated to the nominated
        event (default True). False still collects regional context, extends
        event paths through overlaps and assesses competing ordinary splices;
        skipped regional seeds remain explicit in the result.
    exon_anchored_clips : bool
        Also seed paths from terminal soft clips on versioned reference exons
        and extend them through exact sequence-only overlaps. Requires supplied
        reference models and ``use_soft_clipped_bases``. This exploratory mode
        can find mature junctions away from DNA breakpoints; a bridge alone
        does not establish event causation or a resolved isoform.
    breakpoint_window : int
        Bases searched to either side of each nominated breakpoint.
    max_breakpoint_shift : int
        Junction bases an aligner may place past a breakpoint (homology, not
        verified against a genome) while the join still counts as the
        breakpoint junction.
    annotated_junction_tolerance : int
        Unannotated joins within this many bases of an annotated junction are
        treated as alignment wobble and seed no ordinary-splicing or regional
        path (joins at the nominated adjacency are always used).
    peptide_lengths : iterable of int
    min_orf_amino_acids, max_orf_candidates : int
        Minimum length and per-path cap for separate exploratory ATG ORFs.
        Candidates are ordered by start offset; truncation is explicit.
    inclusion_min_splice_anchor_bases, inclusion_min_base_quality, inclusion_min_mapping_quality : int
        Gates for linking an intronic ORF start to an observed splice
        (priority 3); the MAPQ gate also selects the competing splices used
        for intron retention. See ``annotate_orf_inclusion``.
    inclusion_mapq_255_is_unique : bool
        Treat MAPQ 255 as a unique alignment, as STAR writes it (the default),
        in those gates. False treats it as unavailable, per the SAM
        specification; excluded 255 records are then reported as
        ``mapq_255_excluded_from_inclusion``.
    dense_support : bool
        False (default) searches a bounded set of regional reads. True uses
        temporary disk space to search large inputs: automatically group reads
        by full sequence and alignment, share repeated calculations within
        max_records, then count candidate support across all eligible original
        reads. Exon-anchored clip seeds, when enabled, are recounted separately
        from direct junctions. Grouping needs no separate tuning. Discovery and intron-competition
        limits remain explicit. An optional evidence-group cap can further
        bound the set of exact structural/allelic hypotheses searched.
    scratch_dir : str or pathlib.Path, optional
        Optional location for temporary files when dense_support is True.
        None uses the system temporary directory. Set an existing directory on
        another disk when more space is needed; this does not change results.
        A temporary original-record index is built here unless ``library_index``
        is supplied. Discovery work remains private to the current event view.
    library_index : RnaLibraryIndex, optional
        Open context from ``prepare_rna_library(...).open()`` for this same BAM
        and read eligibility policy. Reuse it across dense event views to avoid
        repeated whole-library preparation. Only complete original templates are
        acquired; exact CIGAR/sequence evidence is still checked for each event.
        Requires ``dense_support=True``. The caller owns its lifetime.
    read_filter_identity : str, optional
        Version/provenance of custom original-record eligibility, matching the
        reusable index. Includes custom predicates and overridden collector
        eligibility methods; update it when behavior or closed-over data changes.

    Returns
    -------
    dict
        JSON-serializable exploratory candidates with independent sequence,
        frame and event-linkage evidence, original records and limitations.
        This is not the validated supplied-fusion result consumed by Vaxrank.
    """
    options = locals().copy()
    if not isinstance(dense_support, bool):
        raise ValueError("dense_support must be True or False")
    if scratch_dir is not None and not dense_support:
        raise ValueError("scratch_dir requires dense_support")
    if not dense_support and (library_index is not None or read_filter_identity is not None):
        raise ValueError("library_index and read_filter_identity require dense_support")
    if max_evidence_groups is not None:
        if not dense_support:
            raise ValueError("max_evidence_groups requires dense_support")
        if type(max_evidence_groups) is not int or max_evidence_groups < 1:
            raise ValueError("max_evidence_groups must be a positive integer")
    if not all((event_id, reference_name, sample_id, source)) or not isinstance(event_provenance, dict) \
            or not event_provenance:
        raise ValueError("SV reconstruction requires event, reference, sample, source and event provenance")
    if not isinstance(donor, FusionBreakpoint) or not isinstance(acceptor, FusionBreakpoint):
        raise TypeError("Expected oriented FusionBreakpoint inputs")
    for value in (min_anchor_bases, min_overlap, min_alternative_fragments, max_records, max_queries, max_paths,
                  max_extension_segments, breakpoint_window, min_orf_amino_acids, max_orf_candidates,
                  inclusion_min_splice_anchor_bases, inclusion_min_base_quality, inclusion_min_mapping_quality):
        if type(value) is not int or value < 1:
            raise ValueError("SV evidence thresholds and search limits must be positive integers")
    if not isinstance(inclusion_mapq_255_is_unique, bool):
        raise ValueError("inclusion_mapq_255_is_unique must be True or False")
    inclusion_thresholds = dict(min_splice_anchor_bases=inclusion_min_splice_anchor_bases,
                                min_base_quality=inclusion_min_base_quality,
                                min_mapping_quality=inclusion_min_mapping_quality,
                                mapq_255_is_unique=inclusion_mapq_255_is_unique)
    for value in (max_breakpoint_shift, annotated_junction_tolerance):
        if type(value) is not int or value < 0:
            raise ValueError("Breakpoint shift and annotated-junction tolerance must be nonnegative integers")
    for value in (min_alternative_fraction, min_local_variant_fraction):
        if isinstance(value, bool) or not 0 <= value <= 1:
            raise ValueError("Alternative-support fractions must be between 0 and 1")
    if min_anchor_bases < 3:
        raise ValueError("A frame anchor needs at least one codon")
    lengths = tuple(sorted(set(peptide_lengths)))
    if not lengths or any(type(n) is not int or n < 1 for n in lengths):
        raise ValueError("Peptide lengths must be positive integers")
    if type(include_regional_candidates) is not bool:
        raise ValueError("include_regional_candidates must be a boolean")
    if type(exon_anchored_clips) is not bool:
        raise ValueError("exon_anchored_clips must be a boolean")
    references = tuple(references)
    if any(not isinstance(r, FusionReference) or r.reference_name != reference_name
           for r in references):
        raise ValueError("Matching, explicitly versioned reference models are required")
    if len({(r.annotation, r.transcript_id) for r in references}) != len(references):
        raise ValueError("Duplicate reference transcript identity")
    collector = ReadCollector() if read_collector is None else read_collector
    if exon_anchored_clips and (not references or not collector.use_soft_clipped_bases):
        raise ValueError("Exon-anchored clips require reference models and soft-clipped bases")
    if dense_support:
        from .sv_rna_dense import reconstruct_dense_sv_rna

        options.pop("bam")
        options.update(references=references, regions=tuple(regions), read_collector=collector,
                       peptide_lengths=lengths)
        return reconstruct_dense_sv_rna(bam, options)
    models = [TranscriptPathModel(r) for r in references]
    annotated = set().union(*(m.junctions() for m in models))
    reference_peptides = {m.protein[i:i + k] for m in models if m.protein
                          for k in lengths for i in range(len(m.protein) - k + 1)}
    event = OrientedAdjacency(donor, acceptor, max_breakpoint_shift, annotated, annotated_junction_tolerance)
    spans = defaultdict(list)
    exon_spans = defaultdict(list)
    for r in references:
        spans[r.contig, r.strand].append((r.exons[0][0], r.exons[-1][1]))
        exon_spans[r.contig, r.strand].extend(r.exons)
    if exon_anchored_clips:
        exon_spans = {key: merge_intervals(intervals) for key, intervals in exon_spans.items()}

    def anchored(position):
        """Within an annotated gene span, on its strand."""
        return any(a <= position[1] < b for a, b in spans.get((position[0], position[2]), ()))

    def exonic(position):
        return interval_contains(exon_spans.get((position[0], position[2]), ()), position[1])

    records, acquisition = collect_sv_records(bam, regions, references, max_records, max_queries,
                                              (donor, acceptor), breakpoint_window)
    header = bam.header.to_dict()
    store = SvObservationStore(records, collector, max_extension_segments, header, sample_id, source)
    store.sequence_overlaps = exon_anchored_clips
    inclusion_splices, mapq_255_excluded = defaultdict(set), False
    for identity, members in store.groups.items():
        for read in members:
            if read.mapping_quality == 255 and not inclusion_mapq_255_is_unique:
                mapq_255_excluded = True
                continue
            if read.mapping_quality != 255 and read.mapping_quality < inclusion_min_mapping_quality:
                continue
            for left, right, kind, _ in alignment_gaps(read, collector):
                if kind == "N":
                    inclusion_splices[left, right].add(identity[:2])
                    inclusion_splices[reverse_position(right), reverse_position(left)].add(identity[:2])
    inclusion_splices = [dict(left=list(left), right=list(right), fragments=len(fragments))
                         for (left, right), fragments in sorted(inclusion_splices.items())]
    now, lazy, competing = store.classify(event, annotated, anchored,
                                         exonic if exon_anchored_clips else None)
    initial = [o for identity in sorted(now) for o in store.build(identity)]
    first = sv_path_seeds(initial, event, annotated, anchored, competing, deferred=lazy,
                   exonic=exonic if exon_anchored_clips else None,
                   exon_clipped=store.exon_clipped)
    queue = [(relation, key, ids) for relation in ("splice_ambiguous_event_junction", "regional_novel_junction")
             for key, ids in sorted(lazy[relation].items(),
                                    key=lambda item: (-len({i[:2] for i in item[1]}), repr(item[0])))
             if len({i[:2] for i in ids}) >= min_alternative_fragments]
    thresholds = (min_alternative_fragments, min_alternative_fraction, min_local_variant_fraction)
    paths, pruned, notes, seed_rows, strongest, seen = {}, [], set(acquisition["limitations"]), [], {}, set()
    if not include_regional_candidates:
        seed_rows.extend(dict(relation=relation, left=list(left), right=list(right),
                              fragments=len({i[:2] for i in ids}), unexplored=True, paths=[],
                              reason="regional_candidate_search_disabled")
                         for relation, (left, right), ids in queue if relation == "regional_novel_junction")
        queue = [row for row in queue if row[0] != "regional_novel_junction"]
    if not references:
        notes.add("reference_models_unavailable")
    if mapq_255_excluded:
        notes.add("mapq_255_excluded_from_inclusion")
    reconstruction_scopes = defaultdict(set)
    end_reasons = defaultdict(set)

    def batches():
        """(seeds, gated): event seeds first; other joins are built only if reached."""
        yield first, False
        for k, (_, key, ids) in enumerate(queue):
            if len(paths) >= max_paths:
                notes.add("path_limit")
                seed_rows.extend(dict(relation=relation, left=list(left), right=list(right),
                                      fragments=len({i[:2] for i in ids}), unexplored=True, paths=[])
                                 for relation, (left, right), ids in queue[k:])
                return
            if len(ids) > max_extension_segments:
                notes.add("seed_segment_limit")
            chosen = sorted(ids)[:max_extension_segments]
            yield sv_path_seeds([o for identity in chosen for o in store.build(identity)], event, annotated, anchored,
                         competing, key=key), True

    for batch, gated in batches():
        for (core_sequence, core_positions), seed in batch:
            if (core_sequence, core_positions) in seen:
                continue
            seen.add((core_sequence, core_positions))
            row = dict(relation=seed["relation"], kinds=sorted(seed["kinds"]),
                       left=list(core_positions[0]) if core_positions[0] else None,
                       right=list(core_positions[-1]) if core_positions[-1] else None,
                       unplaced_bases=core_sequence[1:-1], fragments=len(seed["fragments"]),
                       reads=len({o.identity for o, _ in seed["observations"]}), paths=[])
            if seed["relation"] == "exon_anchored_clip_partner_unplaced":
                row["core_sequence"] = core_sequence
            seed_rows.append(row)
            # Joins at or across the event may rest on one fragment; ordinary
            # splices and regional joins need support.
            if gated and row["fragments"] < min_alternative_fragments:
                row["below_min_fragments"] = True
                continue
            # Seeds come strongest first. A base variant of a seeded junction is
            # pruned like a weak extension branch (only against a supported
            # best). Another placement of the nominated adjacency (a different
            # junction-base assignment or noisy near-miss) is the same event:
            # it seeds only if it is itself supported.
            variant = "adjacency" if seed["near_breakpoint"] else core_positions
            if variant in strongest:
                best, best_positions = strongest[variant]
                strong = supported_alternatives({("", None): best, ("", "this"): row["fragments"]}, *thresholds)
                if ("", "this") not in strong and (("", None) in strong or core_positions != best_positions):
                    row["minor_variant_of_seeded_junction"] = True
                    continue
            strongest.setdefault(variant, (row["fragments"], core_positions))
            if any(contains_seed(path, (core_sequence, core_positions)) for path in paths):
                row["represented_by_earlier_path"] = True
                continue
            if len(paths) >= max_paths:
                notes.add("path_limit")
                row["unexplored"] = True
                continue
            spanning = [o for o, _ in seed["observations"]]
            mirrors = {o.key: m for o in spanning for m in store.built[o.identity]
                       if m.key[:-1] == o.key[:-1] and m.key != o.key}
            pool = ObservationIndex(spanning + [mirrors[o.key] for o in spanning])
            # Regional assembly must not erase a seed-spanning reconstruction:
            # e.g. a 3' extender can occupy an upstream placement and the cycle
            # guard then prevents recovering the witnessed 5' sequence. Keep
            # both hypotheses when neither contains the other, without letting
            # one long junction read veto a well-supported regional isoform.
            for scope, index in ([("regional", store), ("seed_spanning", pool)] if assemble
                                 else [("seed_spanning", pool)]):
                scope_pruned = []
                for sequence, positions, voters, ends in extend_path_both_directions(
                        core_sequence, core_positions, index, min_overlap, thresholds,
                        max(1, max_paths - len(paths)), scope_pruned, notes):
                    path = sequence, positions
                    represented = [existing for existing in paths if contains_seed(existing, path)]
                    if represented:
                        for existing in represented:
                            paths[existing][0].append(row)
                            paths[existing][1].update(voters)
                            reconstruction_scopes[existing].add(scope)
                        continue
                    if path not in paths and len(paths) >= max_paths:
                        notes.add("path_limit")
                        break
                    rows, contributors = paths.setdefault(path, ([], set()))
                    rows.append(row)
                    contributors.update(voters)
                    reconstruction_scopes[path].add(scope)
                    end_reasons[path].update(ends.items())
                # Both scopes can prune exactly the same alternative. These
                # are descriptions of decisions, not additional evidence.
                for branch in scope_pruned:
                    if branch not in pruned:
                        pruned.append(branch)
    if store.limited:
        notes.add("extension_segment_limit")
    adjacency = [entry for (left, right), entries in store.joins.items() for entry in entries
                 if event.relation(left, right, entry[2] - entry[1] - 1)[0] == "breakpoint_junction"]
    clip_support = defaultdict(list)
    for o in store.observations.values():
        for i, j, assignment in event.clips(o.positions, len(o.sequence)):
            clip_support["donor" if assignment[0] == 0 else "acceptor"].append((o, i, j, "clip"))
    results, used, junction_cache, witness_groups = [], {}, {}, {}
    for (sequence, positions), (rows, voters) in paths.items():
        result, evidence = build_sv_path_result(
            sequence, positions, voters, store, event, annotated, adjacency, clip_support,
            lambda *path: sv_frame_evidence(*path, models, min_anchor_bases, reference_peptides, lengths),
            lambda *path: exploratory_orfs(*path, references, store.cell_umi.support, min_orf_amino_acids, max_orf_candidates,
                                           lineage=store.lineage.support, competing_splices=inclusion_splices,
                                           inclusion_thresholds=inclusion_thresholds,
                                           witness_groups=witness_groups), junction_cache=junction_cache)
        results.append(result)
        result["reconstruction_scopes"] = sorted(reconstruction_scopes[sequence, positions])
        # One path can be reached from several seeds or scopes, so an end can
        # have more than one reason; a natural end is not a truncation.
        result["end_reasons"] = {end: sorted(r for e, r in end_reasons[sequence, positions] if e == end)
                                 for end in ("5prime", "3prime")}
        for row in rows:
            if result["path_id"] not in row["paths"]:
                row["paths"].append(result["path_id"])
        used.update(evidence)
    rank = {status: i for i, status in enumerate(LINKAGE_STATUSES)}
    results.sort(key=lambda r: (rank.get(r["event_linkage"]["status"], len(rank)),
                                -max([j["direct_support"]["fragments"] for j in r["junctions"]], default=0),
                                r["sequence"], r["path_id"]))
    record_ids = {rid for o in used.values() for rid in o.records}
    cited_reads = {rid: r for identity in {o.identity for o in used.values()}
             for r in store.groups[identity] if (rid := sam_record_id(r)) in record_ids}
    best = next((s for s in LINKAGE_STATUSES if s in {r["event_linkage"]["status"] for r in results}), None)
    return dict(
        schema="isovar.sv_rna_candidates.v5",
        status=("event_linked_candidates" if best in EVENT_LINKED_RELATIONS
                else "splice_ambiguous_candidates" if best == "splice_ambiguous_event_junction"
                else "regional_candidates_only" if results else "no_candidate_paths"),
        event_id=event_id, reference_name=reference_name, sample_id=sample_id, source=source,
        alignment_metadata=dict(read_groups=header.get("RG", []), programs=header.get("PG", [])),
        read_lineage=dict(scope="eligible_records_within_input_read_group", reads=store.lineage.evidence()),
        cell_umi_evidence=store.cell_umi.evidence(),
        event_provenance=event_provenance, donor=asdict(donor), acceptor=asdict(acceptor),
        reference_models=[dict(transcript_id=r.transcript_id, annotation=r.annotation, contig=r.contig,
                               strand=r.strand, sequence_sha256=sha256(r.sequence.encode()).hexdigest())
                          for r in references],
        acquisition=acquisition, limitations=sorted(notes),
        parameters=dict(min_anchor_bases=min_anchor_bases, min_overlap=min_overlap,
                        min_alternative_fragments=min_alternative_fragments,
                        min_alternative_fraction=min_alternative_fraction,
                        min_local_variant_fraction=min_local_variant_fraction, max_records=max_records,
                        max_queries=max_queries, max_paths=max_paths,
                        max_extension_segments=max_extension_segments, assemble=assemble,
                        include_regional_candidates=include_regional_candidates,
                        exon_anchored_clips=exon_anchored_clips,
                        breakpoint_window=breakpoint_window, max_breakpoint_shift=max_breakpoint_shift,
                        annotated_junction_tolerance=annotated_junction_tolerance,
                        peptide_lengths=lengths, genetic_code=1,
                        min_orf_amino_acids=min_orf_amino_acids, max_orf_candidates=max_orf_candidates,
                        inclusion=dict(inclusion_thresholds),
                        read_collection=collector.collection_parameters()),
        seeds=seed_rows,
        competing_annotated_junctions=[dict(left=list(left), right=list(right), fragments=len(fragments))
                                       for (left, right), fragments in sorted(competing.items())],
        pruned_branches=pruned, paths=results,
        observations={key: dict(identity=list(o.identity), records=list(o.records), reverse_complement=o.reverse_complement,
                                original_query_interval=list(o.query_interval),
                                missing_qualities=o.missing_qualities, secondary=o.secondary)
                      for key, o in sorted(used.items())},
        observation_counts=dict(eligible_reads=len(store.groups), built_reads=len(store.built)),
        excluded_records=dict(store.excluded), read_path_notes=dict(store.reasons),
        record_evidence={rid: record_evidence(read) for rid, read in sorted(cited_reads.items())},
        original_records={rid: read.to_string() for rid, read in sorted(cited_reads.items())})


SV_RNA_INPUT_KEYS = frozenset((
    "event_id", "reference_name", "sample_id", "donor", "acceptor", "regions", "references", "event_provenance"))


def sv_rna_input_from_dict(data):
    """Decode the explicit JSON input used by ``isovar sv-rna``.

    Returns keyword arguments for ``reconstruct_sv_rna`` (other than the BAM,
    source and evidence parameters). Raises ValueError for a missing key or
    an unexpected field, at any level.
    """
    unknown = set(data) - SV_RNA_INPUT_KEYS
    if unknown:
        raise ValueError("Unexpected SV RNA input keys: %s" % ", ".join(sorted(unknown)))
    try:
        return _sv_rna_input_from_dict(data)
    except KeyError as error:
        raise ValueError("SV RNA input is missing %s" % error) from error
    except TypeError as error:
        raise ValueError("Invalid SV RNA input: %s" % error) from error


def _sv_rna_input_from_dict(data):
    return dict(
        event_id=data["event_id"], reference_name=data["reference_name"], sample_id=data["sample_id"],
        donor=FusionBreakpoint(**data["donor"]), acceptor=FusionBreakpoint(**data["acceptor"]),
        regions=[tuple(region) for region in data.get("regions", [])],
        references=[FusionReference(**r) for r in data["references"]],
        event_provenance=data["event_provenance"])

# Stable public entry points; shared implementations have category-neutral homes.
