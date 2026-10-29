"""Original RNA paths around literal alleles, including intronic nominations.

Discovery is independent of coding annotation. Original allele witnesses seed
shared RNA reconstruction; supplied transcript models compete afterward. A
coding prediction is conditional RNA evidence, not evidence of translation.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, replace
from hashlib import sha256
import logging
from tempfile import TemporaryDirectory
from typing import Any

import pysam
from pyensembl import Genome
from varcode import ContigResolver, Variant

from .allele_path_seeds import AllelePathSeed, AllelePathSeeds, allele_path_seeds
from .contig_resolution import require_resolved_contig, variant_reference_name, variant_source_contig
from .fusion import FusionReference
from .read_collector import ReadCollector
from .rna_library import RnaLibraryIndex, _original_library_context
from .rna_path_linkage import PathFeature
from .rna_reconstruction import _discover, _matches, _support
from .sv_rna import _Positions, _breaks, _collection_parameters, _runs
from .transcript_path_frames import Position, TranscriptPathModel, assess_transcript_path
from .variant_helpers import LiteralVariantLocus, require_literal_variant


Region = tuple[str, int, int]
ObservedPath = tuple[str, tuple[Position | None, ...]]
logger = logging.getLogger(__name__)


def variant_path_read_collector(settings: Mapping[str, Any] | None = None, *, files=None) -> ReadCollector:
    """Shared original-path defaults, with explicit cohort/API policy overrides."""
    configuration = dict(use_secondary_alignments=True, min_mapping_quality=0,
                         use_soft_clipped_bases=False, use_reads_without_base_qualities=True)
    configuration.update(settings or {})
    return ReadCollector.from_dict(configuration, files=files)


def _seed_offsets(sequence: str, positions: tuple, seeds: Iterable[AllelePathSeed]) -> list[tuple[int, AllelePathSeed]]:
    """Retain full original allele-bearing geometries, never sequence-only hits."""
    return [(low, seed) for seed in seeds for low, _ in
            _matches(sequence, positions, seed.observation, contained=True)]


def _permitted_departures(seed_offsets: Sequence[tuple[int, AllelePathSeed]]) -> dict[int, str]:
    """Only observed allele/splice/indel transitions justify departing frames."""
    departures = {}
    for offset, seed in seed_offsets:
        for left, right, kind in seed.observation.breaks:
            departures[offset + left + 1] = "original_RNA_" + kind
        start, end = seed.query_interval
        if start == end or any(p is None for p in seed.observation.positions[start:end]):
            departures[offset + start] = "observed_literal_indel"
    return departures


def _path_candidate(sequence: str, positions: tuple, seeds: AllelePathSeeds,
                    models: Sequence[TranscriptPathModel], minimum: int) -> dict:
    positions = _Positions(positions)
    matches = _seed_offsets(sequence, positions, seeds.seeds)
    identifier = "isovar_allele_path_" + sha256(repr((sequence, positions)).encode()).hexdigest()[:24]
    focal: dict[tuple[str, tuple[int, int], tuple[int, int]], list[AllelePathSeed]] = defaultdict(list)
    for offset, seed in matches:
        observed = offset + seed.query_interval[0], offset + seed.query_interval[1]
        canonical = offset + seed.canonical_query_interval[0], offset + seed.canonical_query_interval[1]
        focal[seed.state, observed, canonical].append(seed)
    features: list[PathFeature] = []
    focal_rows: list[dict[str, Any]] = []
    frames: list[dict[str, Any]] = []
    coding: list[dict[str, Any]] = []
    departures = {left + 1: "assembled_path_discontinuity" for left, _ in _breaks(positions)}
    departures.update(_permitted_departures(matches))
    for number, (left, right) in enumerate(_breaks(positions)):
        features.append(PathFeature("junction_%d" % number, left, right + 1))
    for index, ((state, observed, canonical), witnesses) in enumerate(sorted(focal.items())):
        feature_id = "allele_%d" % index
        features.append(PathFeature(feature_id, observed[0], observed[1], tuple(sorted({s.observation.key for s in witnesses}))))
        focal_rows.append(dict(id=feature_id, state=state, query_interval=observed,
            canonical_query_interval=canonical, segments=sorted({s.observation.identity for s in witnesses}),
            record_ids=sorted({s.record_id for s in witnesses}),
            observation_keys=sorted({s.observation.key for s in witnesses})))
        for model in models:
            assessed = assess_transcript_path(sequence, positions, model, min_anchor_bases=minimum,
                permitted_departures=departures, allow_retained_introns=True, focal_interval=observed)
            row: dict[str, Any] = dict(transcript_id=model.reference.transcript_id, focal_id=feature_id,
                anchors=[list(a) for a in assessed.anchors], reason=assessed.reason,
                unpermitted_departures=assessed.other_departures, readings=[])
            for reading in assessed.readings + assessed.collinear_readings:
                reading_id = "frame_%d" % len(coding)
                # Export reference protein separately. No reference residue is
                # substituted for an unobserved RNA codon.
                translated = dict(asdict(reading), id=reading_id,
                    transcript_id=model.reference.transcript_id, focal_id=feature_id,
                    reference_protein=model.protein, reference_cds_start=model.reference.cds_start,
                    query_start=reading.translation_start, query_end=reading.translation_end,
                    status="complete_CDS" if reading.complete_5prime and reading.ends_with_stop_codon else "partial_CDS",
                    protein_status="conditional_RNA_frame", focal_state=state)
                coding.append(translated)
                row["readings"].append(reading_id)
                for number, anchor in enumerate(reading.anchors):
                    features.append(PathFeature(reading_id + "_anchor_%d" % number,
                                                anchor.query_start, anchor.query_end))
            frames.append(row)
    return dict(id=identifier, sequence=sequence, placements=positions, origin="observed_assembly",
                focal_alleles=focal_rows, frame_assessments=frames, coding=coding, features=features,
                protein_interpretation="not_assessed" if not models else
                    "conditional_frames_available" if coding else "no_justified_CDS_frame")


def _reference_competitors(candidates: list[dict], models: Sequence[TranscriptPathModel]) -> list[dict]:
    """Count ordinary model-compatible RNA separately from allele-seeded paths.

    Model sequences are labelled competitors, never observed reconstructions.
    Identical observed/model geometries share an ID and support recount.
    """
    by_path = {(row["sequence"], tuple(row["placements"])): row for row in candidates}
    references = []
    for model in models:
        key = model.reference.sequence, model.positions
        if key not in by_path:
            row = dict(id="isovar_allele_path_" + sha256(repr(key).encode()).hexdigest()[:24],
                sequence=key[0], placements=key[1], origin="reference_model", coding=[])
            references.append(row)
            by_path[key] = row
        by_path[key].setdefault("reference_transcript_ids", []).append(model.reference.transcript_id)
    return references


def _resolved_context(bam: pysam.AlignmentFile, variant: Variant, reference_name: str,
                      references: Sequence[FusionReference], regions: Sequence[Region], flank: int
                      ) -> tuple[str, tuple[FusionReference, ...], dict, list[dict]]:
    resolver = ContigResolver(bam.references, reference_name=reference_name)
    diagnostics = []

    def resolve(name: str | int) -> str:
        resolution = resolver.resolve(name)
        diagnostics.append(asdict(resolution))
        return require_resolved_contig(resolution, bam)

    contig = resolve(variant_source_contig(variant))
    lengths = dict(zip(bam.references, bam.lengths))
    locus = LiteralVariantLocus.from_variant(variant)
    if locus.end > lengths[contig] or locus.start >= lengths[contig]:
        raise ValueError("Literal allele lies outside the BAM reference")
    resolved = []
    for reference in references:
        if reference.reference_name != reference_name:
            raise ValueError("Transcript and variant reference declarations differ; name resolution is not liftover")
        reference = replace(reference, contig=resolve(reference.contig))
        if reference.exons[-1][1] > lengths[reference.contig]:
            raise ValueError("Transcript lies outside the BAM reference: " + reference.transcript_id)
        resolved.append(reference)
    intervals: dict[str, set[tuple[int, int, str]]] = defaultdict(set)
    intervals[contig].update((max(0, locus.start - flank), min(lengths[contig], locus.end + flank), strand)
                            for strand in ("+", "-"))
    for name, start, end in regions:
        name = resolve(name)
        if type(start) is not int or type(end) is not int or not 0 <= start < end <= lengths[name]:
            raise ValueError("Search intervals must be nonempty half-open BAM-reference intervals")
        intervals[name].update((start, end, strand) for strand in ("+", "-"))
    for reference in resolved:
        intervals[reference.contig].update((reference.exons[0][0], reference.exons[-1][1], strand)
                                           for strand in ("+", "-"))
    return contig, tuple(resolved), dict(intervals), diagnostics


def reconstruct_variant_rna_paths(
        bam: pysam.AlignmentFile, *, variant: Variant, sample_id: str, source: str,
        references: Sequence[FusionReference] = (), regions: Sequence[Region] = (),
        provenance: Mapping[str, Any] | None = None, read_collector: ReadCollector | None = None,
        max_discovery_reads: int = 10000, max_seed_segments: int = 10000, max_paths: int = 256,
        min_overlap: int = 30, min_anchor_bases: int = 12, flank: int = 1000,
        library_index: RnaLibraryIndex | None = None, read_filter_identity: str | None = None,
        scratch_dir: str | None = None) -> dict[str, Any]:
    """Assemble literal-allele RNA before assessing any supplied coding model.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        One declared original sample/library. Regional inputs retain that scope.
    variant : varcode.Variant
        Literal allele with explicit reference context; no annotation download
        is needed. The original supplied contig is resolved against the BAM.
    sample_id, source : str
        Explicit identities. Different libraries/cultures require separate runs.
    references : sequence of FusionReference
        Versioned coding and noncoding competitors on the same assembly.
        Missing models never prevent RNA reconstruction or imply no possible ORF.
    regions : sequence of (contig, start, end)
        Additional half-open search intervals. The locus flank and full genomic
        spans of supplied models are also searched, in both orientations.
    provenance : mapping, optional
        Caller-supplied source/acquisition/reference checksums and scope.
    read_collector : ReadCollector, optional
        Shared eligibility, base-quality and technical-sequence policy.
    max_discovery_reads, max_seed_segments, max_paths : int
        Positive discovery budgets. ALT and REF use equal independent path caps.
        Allele counts and support use all eligible originals in the input index.
    min_overlap, min_anchor_bases, flank : int
        Exact extension overlap, coding-anchor length and locus-search radius.
    library_index : RnaLibraryIndex, optional
        Reuse an open index of this BAM across events with identical eligibility.
    read_filter_identity : str, optional
        Stable custom-filter identity required for persistent index reuse.
    scratch_dir : str, optional
        Parent directory for private temporary files.

    Returns
    -------
    dict
        Versioned RNA paths, conditional frames, original witnesses, separate
        support/linkage counts and explicit limits. This is not a mutant-protein
        funnel: an ALT RNA path or frame alone does not establish protein change,
        variant causation, mature RNA, translation or end-to-end molecular phase.
    """
    require_literal_variant(variant)
    reference_name = variant_reference_name(variant)
    if not reference_name or not all(isinstance(value, str) and value for value in (sample_id, source)):
        raise ValueError("Explicit variant reference, sample_id and source are required")
    samples = {row["SM"] for row in bam.header.to_dict().get("RG", []) if row.get("SM")}
    if len(samples) > 1:
        raise ValueError("Reconstruct each BAM sample separately")
    for name, value in (("max_discovery_reads", max_discovery_reads), ("max_seed_segments", max_seed_segments),
                        ("max_paths", max_paths), ("min_overlap", min_overlap), ("flank", flank)):
        if type(value) is not int or value < 1:
            raise ValueError(name + " must be a positive integer")
    if type(min_anchor_bases) is not int or min_anchor_bases < 3:
        raise ValueError("min_anchor_bases must be at least three")
    if not all(isinstance(reference, FusionReference) for reference in references):
        raise TypeError("references must contain FusionReference models")
    if len({r.transcript_id for r in references}) != len(references):
        raise ValueError("Transcript IDs must be unique")
    contig, resolved, search, diagnostics = _resolved_context(bam, variant, reference_name, references, regions, flank)
    models = tuple(TranscriptPathModel(reference) for reference in resolved)
    collector = read_collector or variant_path_read_collector()
    locus = LiteralVariantLocus.from_variant(variant)
    with TemporaryDirectory(dir=scratch_dir, prefix="isovar-allele-paths-") as directory, _original_library_context(
            bam, collector, directory, library_index=library_index, reference_name=reference_name,
            read_filter_identity=read_filter_identity) as groups:
        seeds = allele_path_seeds(groups.templates([(contig, max(0, locus.start - 1), locus.end + 1)]),
                                  locus, contig, read_collector=collector)
        logger.info("%s: allele witnesses: %d ALT / %d REF / %d OTHER segments", source,
                    len(seeds.segments("ALT")), len(seeds.segments("REF")), len(seeds.segments("OTHER")))
        paths: dict[ObservedPath, None] = {}
        discoveries: dict[str, Any] = {}
        for state in ("ALT", "REF"):
            identities = sorted(seeds.seed_segments(state))
            selected = set(identities[:max_seed_segments])
            allowed = {s.observation.key: s.observation for s in seeds.seeds
                       if s.state == state and s.observation.identity in selected}
            logger.info("%s: discovering %s paths from %d original seed segments", source, state, len(selected))
            found, discovery = _discover(groups, search, collector, max_discovery_reads, max_paths, min_overlap,
                                          seed_observations=tuple(allowed.values()))
            discovery.update(seed_segments=len(identities), selected_seed_segments=len(selected),
                             seed_limit_reached=len(identities) > len(selected))
            discovery["truncated"] |= discovery["seed_limit_reached"]
            discoveries[state] = discovery
            paths.update((path, None) for path in found)
        candidates = [_path_candidate(sequence, positions, seeds, models, min_anchor_bases)
                      for sequence, positions in sorted(paths, key=repr)]
        reference_candidates = _reference_competitors(candidates, models)
        logger.info("%s: recounting original support for %d retained RNA paths", source, len(candidates))
        support = _support(groups, candidates + reference_candidates, search, collector,
                           bam.header.to_dict(), sample_id, source, label="path")
        statistics, library_identity = groups.statistics, groups.prepared.identity
    for candidate in candidates + reference_candidates:
        candidate["mapping"] = [list(run) for run in _runs(candidate.pop("placements"))]
        if "features" in candidate:
            candidate["features"] = [asdict(feature) for feature in candidate["features"]]
    return dict(schema="isovar.variant_rna_paths.v1", sample_id=sample_id, source=source,
        reference_name=reference_name, variant=dict(contig=variant_source_contig(variant), start=variant.start,
            ref=variant.ref, alt=variant.alt), resolved_contig=contig, contig_diagnostics=diagnostics,
        references=[asdict(reference) for reference in resolved], provenance=dict(provenance or {}),
        input=dict(statistics=statistics, identity=library_identity, header=bam.header.to_dict()),
        collection=_collection_parameters(collector), search_regions={key: sorted(value) for key, value in search.items()},
        allele_counts={state: dict(segments=len(seeds.segments(state)), fragments=len(seeds.fragments(state)))
                       for state in ("ALT", "REF", "OTHER")},
        allele_observation_notes=dict(seeds.notes), original_allele_records=dict(seeds.original_records),
        original_allele_calls=[call._asdict() for call in seeds.calls],
        discovery=discoveries, paths=candidates, reference_paths=reference_candidates, support=support,
        limitations=["exact_overlap_reconstruction", "RNA_does_not_establish_translation",
                     "overlap_assembly_does_not_establish_end_to_end_phase", "unreported_alignments_not_recovered",
                     "protein_change_and_variant_causation_not_established_by_allele_coverage"])


def variant_rna_path_input_from_dict(data: Mapping[str, Any]) -> dict[str, Any]:
    """Decode one strict JSON input for the API, CLI and cohort runner.

    A reference-only PyEnsembl context supplies Varcode's explicit assembly;
    it supplies no annotation, transcript sequence or default human assembly.
    """
    allowed = {"variant", "reference_name", "sample_id", "references", "regions", "provenance"}
    unknown = set(data) - allowed
    if unknown:
        raise ValueError("Unexpected literal RNA input keys: " + ", ".join(sorted(unknown)))
    try:
        row = data["variant"]
        if set(row) != {"contig", "start", "ref", "alt"}:
            raise ValueError("A literal variant needs exactly contig, start, ref and alt")
        genome = Genome(reference_name=data["reference_name"], annotation_name="supplied_RNA_path_models")
        return dict(variant=Variant(**row, genome=genome, normalize_contig_names=False), sample_id=data["sample_id"],
            references=tuple(FusionReference(**reference) for reference in data.get("references", ())),
            regions=tuple(tuple(region) for region in data.get("regions", ())), provenance=data.get("provenance", {}))
    except (TypeError, KeyError) as error:
        raise ValueError("Invalid literal RNA input: " + str(error)) from error
