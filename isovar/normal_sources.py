"""Source-pinned matched-normal assessment using reusable VCF header resolvers."""

from collections.abc import Iterable, Iterator, Sequence
from dataclasses import asdict, dataclass, field
import logging
from pathlib import Path
from typing import Any, Literal, NotRequired, TypedDict

import pysam
from varcode import ContigResolver, Variant
from varcode.reference import infer_genome

from .contig_resolution import variant_source_contig, variant_reference_name
from .audit_bundle import software_identity
from .genomic_intervals import GenomicIntervalIndex
from .normal_evidence import (
    NormalEvidencePolicy, NormalRecordEvidence, NormalState, assess_normal_record, literal_allele_key,
    normal_source_state, reconcile_normal_states,
)
from .normal_regions import GenomicRegion, normal_background_regions, normal_query_interval, resolve_normal_regions
from .sv_rna_audit import file_sha256, read_json, write_json
from .rna_evidence import content_identifier


logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class NormalSource:
    """One caller/processing product and its explicitly reviewed sample assignment.

    ``acquired_regions`` uses actual source-contig names and zero-based half-open
    intervals. None declares the entire file; an empty tuple declares no acquired
    territory. A regional derivative must not masquerade as a complete callset.
    ``processing_group`` identifies shared input evidence, not an independent
    specimen. ``selection_reason`` is required for both selected and excluded files.
    """

    source_id: str
    reference_name: str | None
    evidence_kind: Literal["variant_only", "reference_confidence"]
    selected: bool
    selection_reason: str
    specimen_id: str | None = None
    library_id: str | None = None
    vcf_sample: str | None = None
    caller: str | None = None
    processing_group: str | None = None
    derived_from: tuple[str, ...] = ()
    path: str | None = None
    sha256: str | None = None
    index_path: str | None = None
    index_sha256: str | None = None
    acquired_regions: tuple[tuple[str, int, int], ...] | None = None
    policy: NormalEvidencePolicy = field(default_factory=NormalEvidencePolicy)
    provenance: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.source_id or not self.selection_reason:
            raise ValueError("Normal sources require identity and a selection/exclusion reason")
        if type(self.selected) is not bool:
            raise ValueError("Normal source selected must be a boolean")
        if not isinstance(self.policy, NormalEvidencePolicy):
            raise TypeError("Normal source policy must be a NormalEvidencePolicy")
        if self.evidence_kind not in ("variant_only", "reference_confidence"):
            raise ValueError("Unknown normal evidence kind")
        if self.selected and not all((self.reference_name, self.specimen_id, self.library_id, self.vcf_sample, self.caller,
                                       self.processing_group, self.path, self.sha256, self.index_path, self.index_sha256)):
            raise ValueError("Selected normal sources require reviewed specimen/library/sample/caller assignments and pinned file/index")
        for digest in (self.sha256, self.index_sha256):
            if digest is not None and (len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest)):
                raise ValueError("Normal source SHA-256 digests must contain 64 lowercase hexadecimal characters")
        for contig, start, end in self.acquired_regions or ():
            if not isinstance(contig, str) or not contig or type(start) is not int or type(end) is not int or not 0 <= start < end:
                raise ValueError("Invalid acquired normal-evidence region")


AlleleKey = tuple[str, int, str, str]


class NormalSourceAssessment(TypedDict):
    """One selected product's local evidence, including unavailable queries."""

    source_id: str
    state: NormalState
    interval: NotRequired[list[int]]
    contig_resolution: NotRequired[dict[str, Any]]
    records: NotRequired[list[NormalRecordEvidence]]
    reason: NotRequired[str]
    error: NotRequired[str]
    background_query_interval: NotRequired[list[int]]
    background_interval_fully_acquired: NotRequired[bool]


class NormalAlleleAssessment(TypedDict):
    """Focal allele plus reconciled state and the evidence behind it."""

    contig: str
    start: int
    ref: str
    alt: str
    normalized_allele: list[str | int]
    state: NormalState
    source_assessments: list[NormalSourceAssessment]
    origin_conflict: bool
    nomination: NotRequired[dict[str, Any]]


class MatchedNormalEvidence(TypedDict):
    """Versioned receipt shared by the Python API, CLI and cohort consumer."""

    schema: Literal["isovar.matched_normal.v1"]
    isovar_version: str
    software: dict[str, Any]
    sample_id: str
    reference_name: str
    background_flank: int
    background_regions: list[GenomicRegion]
    background_source_scopes: list[dict[str, Any]]
    source_inventory: list[dict[str, Any]]
    alleles: list[NormalAlleleAssessment]
    background_variants: list[dict[str, Any]]
    interpretation: str
    manifest: NotRequired[dict[str, Any]]


def _background_records(records: Iterable[pysam.VariantRecord], source: NormalSource,
                        territory: GenomicIntervalIndex
                        ) -> Iterator[tuple[AlleleKey, dict[str, Any]]]:
    """Extract only accepted, genotype-indexed literal alleles with provenance."""
    assert source.vcf_sample is not None  # selected sources were validated on construction
    for record in records:
        reference = record.ref
        if not reference:
            raise ValueError("Normal VCF record requires a nonempty REF allele")
        literal_interval = (record.contig, record.start, record.start + len(reference))
        if not territory.overlaps(*literal_interval):
            continue
        if not any(index is not None and index > 0 for index in record.samples[source.vcf_sample].get("GT") or ()):
            continue
        evidence = assess_normal_record(record, source.vcf_sample, None, source.policy)
        if not evidence["accepted"]:
            continue
        for index in evidence["called_alternate_indexes"]:
            alternate = evidence["alternates"][index - 1]
            if not set(reference + alternate) <= set("ACGT"):
                continue
            key = literal_allele_key(record.contig, record.pos, reference, alternate,
                                     reference_name=source.reference_name)
            yield key, dict(source_id=source.source_id, specimen_id=source.specimen_id,
                library_id=source.library_id, processing_group=source.processing_group,
                sample=source.vcf_sample,
                record={name: value for name, value in evidence.items() if name != "target_allele_called"},
                called_allele_index=index)


def _assess_source_allele(reader: pysam.VariantFile, source: NormalSource, sample_id: str,
                          resolver: ContigResolver, acquired: GenomicIntervalIndex | None,
                          variant: Variant, background_flank: int) -> NormalSourceAssessment:
    """Assess one focal query without changing independent background selection."""
    resolution = resolver.resolve(variant_source_contig(variant))
    interval = normal_query_interval(variant)
    row: NormalSourceAssessment = dict(source_id=source.source_id, state="unavailable", interval=list(interval),
                                        contig_resolution=asdict(resolution))
    contig = resolution.resolved_name
    if contig is None:
        row["reason"] = "unresolved_chromosome"
        return row
    contig_length = reader.header.contigs[contig].length
    if contig_length is not None and interval[1] > contig_length:
        row["reason"] = "outside_declared_contig_length"
        return row
    if acquired is not None and not acquired.covers(contig, *interval):
        row["reason"] = "outside_acquired_regions"
        return row
    try:
        records = list(reader.fetch(contig, *interval))
        assessed = [assess_normal_record(record, sample_id, variant, source.policy) for record in records]
        state = normal_source_state(assessed, interval,
            reference_confidence=source.evidence_kind == "reference_confidence")
        row["state"], row["records"] = state, assessed
        left, right = max(0, interval[0] - background_flank), interval[1] + background_flank
        if contig_length is not None:
            right = min(right, contig_length)
        row["background_query_interval"] = [left, right]
        row["background_interval_fully_acquired"] = acquired is None or acquired.covers(contig, left, right)
    except (OSError, ValueError) as error:
        row["state"], row["reason"], row["error"] = "unavailable", "normal_query_failed", str(error)
    return row


def _source_assessments(source: NormalSource, variants: Sequence[Variant], background_flank: int,
                        background_regions: Sequence[GenomicRegion]
                        ) -> tuple[list[NormalSourceAssessment], list[tuple[AlleleKey, dict[str, Any]]], dict[str, Any]]:
    """Open/index each selected source once; each query retains its actual scope."""
    if not source.selected or source.reference_name is None or source.path is None:
        raise ValueError("Source assessment requires a selected source with a declared assembly")
    if file_sha256(source.path) != source.sha256:
        raise ValueError("Normal source digest mismatch: " + source.source_id)
    if file_sha256(source.index_path) != source.index_sha256:
        raise ValueError("Normal source index digest mismatch: " + source.source_id)
    sample_id = source.vcf_sample
    if sample_id is None:
        raise ValueError("Selected source has no normal VCF sample column")
    acquired = GenomicIntervalIndex(source.acquired_regions) if source.acquired_regions is not None else None
    with pysam.VariantFile(source.path, index_filename=source.index_path) as reader:
        if sample_id not in reader.header.samples:
            raise ValueError("Declared normal VCF sample is absent: " + source.source_id)
        resolver = ContigResolver((str(name) for name in reader.header.contigs), reference_name=source.reference_name)
        rows = [_assess_source_allele(reader, source, sample_id, resolver, acquired, variant, background_flank)
                for variant in variants]
        scope = resolve_normal_regions(background_regions, reader.header, reference_name=source.reference_name)
        requested = GenomicIntervalIndex(scope.regions)
        territory = requested if acquired is None else requested.intersection(acquired)
        # A single sequential background pass is independent of focal query
        # success, ordering, and overlapping neighborhoods. Preserve complete
        # genotyped records selected by REF-span overlap; this is not a claim of
        # reference confidence across the record's full length.
        reader.reset()
        backgrounds = list(_background_records(reader, source, territory))
        background_scope = dict(source_id=source.source_id, resolved_region_count=len(scope.regions),
            resolved_bases=sum(end - start for _, start, end in scope.regions),
            acquired_requested_bases=sum(end - start for _, start, end in territory.intervals()),
            fully_acquired_region_count=sum(acquired is None or acquired.covers(*region) for region in scope.regions),
            unresolved_contigs=list(scope.unresolved_contigs), outside_contig_regions=list(scope.outside_contig_regions),
            clipped_regions=list(scope.clipped_regions))
    return rows, backgrounds, background_scope


def assess_matched_normal(variants: Sequence[Variant], sources: Sequence[NormalSource], *,
                          sample_id: str, reference_name: str, background_flank: int = 0,
                          additional_background_regions: Iterable[GenomicRegion] = ()) -> MatchedNormalEvidence:
    """Assess literal alleles against an explicit normal inventory and policy.

    Parameters
    ----------
    variants : sequence of varcode.Variant
        Unique focal alleles in the declared reference assembly.
    sources : sequence of NormalSource
        Exhaustive reviewed inventory, including excluded/unavailable products.
        Selected sources must be local, indexed VCF/BCF files with pinned hashes.
    sample_id : str
        Tumor/patient identity to which this matched-normal scope is assigned.
    reference_name : str
        Explicit common reference assembly; contig spelling is not liftover.
    background_flank : int, optional
        Additional bases on each side of the focal intervals for background calls.
    additional_background_regions : iterable of (str, int, int), optional
        Explicit contig, zero-based start and half-open end intervals in the same
        assembly. Include distant exons or ncORF blocks needed by reconstructed
        transcripts. Background records are selected by literal REF-span overlap
        with this territory or focal neighborhoods, and acquired territory.
        This selection does not establish an entire normal proteome.

    Returns
    -------
    MatchedNormalEvidence
        Versioned source inventory, per-allele states and original record evidence,
        genotype-indexed background alleles and their processing provenance.
        A confident reference call is local evidence, not proof of tumor specificity.
    """
    from . import __version__
    if not sample_id or not reference_name or type(background_flank) is not int or background_flank < 0:
        raise ValueError("Normal assessment requires identities and a nonnegative background flank")
    reference_name = infer_genome(reference_name)[0].reference_name
    if len({source.source_id for source in sources}) != len(sources):
        raise ValueError("Duplicate normal source identity")
    if any(source.reference_name is None or infer_genome(source.reference_name)[0].reference_name != reference_name
           for source in sources if source.selected):
        raise ValueError("Selected normal sources have incompatible reference assemblies")
    if any(variant_reference_name(variant) != reference_name for variant in variants):
        raise ValueError("Focal variants and normal evidence have incompatible reference assemblies")
    if any(set(variant.ref + variant.alt) - set("ACGT") for variant in variants):
        raise ValueError("Normal assessment requires literal A/C/G/T alleles")
    keys = [literal_allele_key(v.contig, v.start, v.ref, v.alt, reference_name=reference_name) for v in variants]
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate focal allele in normal assessment")
    regions = normal_background_regions(variants, background_flank=background_flank,
                                         additional_regions=additional_background_regions)
    source_rows = []
    background_scopes = []
    background: dict[AlleleKey, dict[tuple[str, str, int], dict[str, Any]]] = {}
    for source in sources:
        if not source.selected:
            # Exclusions apply to the whole run. Keep them once in the source
            # inventory, rather than copying them into every allele assessment.
            continue
        logger.info("Assessing normal source %s (%s): %d focal alleles", source.source_id, source.caller, len(variants))
        rows, observations, scope = _source_assessments(source, variants, background_flank, regions)
        logger.info("Completed normal source %s: %d local assessments, %d background observations",
                    source.source_id, len(rows), len(observations))
        source_rows.append(rows)
        background_scopes.append(scope)
        for key, evidence in observations:
            background.setdefault(key, {})[(source.source_id, content_identifier("normal_record", evidence["record"]),
                                            evidence["called_allele_index"])] = evidence
    alleles: list[NormalAlleleAssessment] = []
    for index, (variant, key) in enumerate(zip(variants, keys)):
        outcomes = [rows[index] for rows in source_rows]
        state = reconcile_normal_states(row["state"] for row in outcomes)
        alleles.append(dict(contig=variant.original_contig, start=variant.original_start,
                            ref=variant.original_ref, alt=variant.original_alt, normalized_allele=list(key),
                            state=state, source_assessments=outcomes,
                            origin_conflict=state in ("alternate_observed", "conflicting_calls")))
    return dict(schema="isovar.matched_normal.v1", isovar_version=__version__, software=software_identity(), sample_id=sample_id,
                reference_name=reference_name, background_flank=background_flank,
                background_regions=list(regions), background_source_scopes=background_scopes,
                source_inventory=[asdict(source) for source in sources], alleles=alleles,
                background_variants=[dict(contig=key[0], start=key[1], ref=key[2], alt=key[3],
                    origins=list(observations.values())) for key, observations in sorted(background.items())],
                interpretation="Local normal genotype/reference-confidence evidence; processing products are not independent replicates. "
                    "Per-allele assessments cover selected sources; exclusions and unavailable products remain in the source inventory.")


def run_normal_evidence_manifest(manifest_path: str | Path, output_path: str | Path) -> MatchedNormalEvidence:
    """Run the shared assessment from pinned inputs and write a new JSON receipt.

    Paths are relative to the manifest, not the working directory. The manifest
    identifies a hashed JSON/JSON.gz nomination list using cohort literal-allele
    fields. Original nominations are retained beside their new normal evidence;
    completed cohort receipts and input callsets are never rewritten.
    """
    manifest_path, output_path = Path(manifest_path).resolve(), Path(output_path).resolve()
    manifest = read_json(manifest_path)
    if manifest.get("schema") != "isovar.normal_manifest.v1":
        raise ValueError("Unsupported normal evidence manifest schema")
    root = manifest_path.parent
    nomination_file = (root / manifest["variants"]["path"]).resolve()
    if file_sha256(nomination_file) != manifest["variants"]["sha256"]:
        raise ValueError("Normal assessment nomination digest mismatch")
    sources = []
    inputs = {manifest_path, nomination_file}
    additional_regions = []
    if descriptor := manifest.get("additional_background_regions"):
        region_file = (root / descriptor["path"]).resolve()
        if file_sha256(region_file) != descriptor["sha256"]:
            raise ValueError("Normal background-region digest mismatch")
        inputs.add(region_file)
        additional_regions = read_json(region_file)
    for descriptor in manifest["sources"]:
        descriptor = dict(descriptor)
        descriptor["policy"] = NormalEvidencePolicy(**descriptor.get("policy", {}))
        descriptor["derived_from"] = tuple(descriptor.get("derived_from", ()))
        if descriptor.get("acquired_regions") is not None:
            descriptor["acquired_regions"] = tuple(tuple(region) for region in descriptor["acquired_regions"])
        for field_name in ("path", "index_path"):
            if descriptor.get(field_name):
                path = (root / descriptor[field_name]).resolve()
                inputs.add(path)
                descriptor[field_name] = str(path)
        sources.append(NormalSource(**descriptor))
    if output_path in inputs:
        raise ValueError("Normal assessment output must not overwrite an input")
    if output_path.exists():
        raise ValueError("Normal assessment output already exists; choose a new receipt path")
    nominations = read_json(nomination_file)
    reference = manifest["reference_name"]
    variants = [Variant(row["contig"], row["start"], row["ref"], row["alt"], genome=reference,
                        convert_ucsc_contig_names=True) for row in nominations]
    result = assess_matched_normal(variants, sources, sample_id=manifest["sample_id"], reference_name=reference,
                                   background_flank=manifest.get("background_flank", 0),
                                   additional_background_regions=additional_regions)
    for allele, nomination in zip(result["alleles"], nominations):
        allele["nomination"] = nomination
    result["manifest"] = dict(path=str(manifest_path), sha256=file_sha256(manifest_path),
                               snapshot=manifest.get("snapshot"), variants=manifest["variants"],
                               additional_background_regions=manifest.get("additional_background_regions"))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_json(output_path, result)
    return result
