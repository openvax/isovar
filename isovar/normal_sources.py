"""Source-pinned matched-normal assessment using reusable VCF header resolvers."""

from collections.abc import Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Literal, NotRequired, TypedDict

import pysam
from varcode import ContigResolver, Variant

from .contig_resolution import variant_source_contig, variant_reference_name
from .evidence_report import changed_reference_interval
from .genomic_intervals import GenomicIntervalIndex
from .normal_evidence import (
    NormalEvidencePolicy, NormalRecordEvidence, NormalState, assess_normal_record, literal_allele_key,
    normal_source_state, reconcile_normal_states,
)
from .sv_rna_audit import file_sha256, read_json, write_json
from .rna_evidence import content_identifier


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


def normal_query_interval(variant: Variant) -> tuple[int, int]:
    """Assess changed bases, or both flanks of an insertion, in zero-based coordinates."""
    start, end, _ = changed_reference_interval(variant.original_start, variant.original_ref, variant.original_alt)
    return start - 1, end


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
    sample_id: str
    reference_name: str
    background_flank: int
    source_inventory: list[dict[str, Any]]
    alleles: list[NormalAlleleAssessment]
    background_variants: list[dict[str, Any]]
    interpretation: str
    manifest: NotRequired[dict[str, Any]]


def _background_records(records: Sequence[pysam.VariantRecord], source: NormalSource, target: Variant,
                        policy: NormalEvidencePolicy) -> list[tuple[AlleleKey, dict[str, Any]]]:
    """Extract only accepted, genotype-indexed literal alleles with provenance."""
    assert source.vcf_sample is not None  # selected sources were validated on construction
    result = []
    for record in records:
        evidence = assess_normal_record(record, source.vcf_sample, target, policy)
        if not evidence["accepted"]:
            continue
        for index in evidence["called_alternate_indexes"]:
            alternate = record.alts[index - 1]
            if not set(record.ref + alternate) <= set("ACGT"):
                continue
            key = literal_allele_key(target.contig, record.pos, record.ref, alternate)
            result.append((key, dict(source_id=source.source_id, specimen_id=source.specimen_id,
                library_id=source.library_id, processing_group=source.processing_group,
                sample=source.vcf_sample,
                record={name: value for name, value in evidence.items() if name != "target_allele_called"},
                called_allele_index=index)))
    return result


def _source_assessments(source: NormalSource, variants: Sequence[Variant], background_flank: int
                        ) -> tuple[list[NormalSourceAssessment], list[tuple[AlleleKey, dict[str, Any]]]]:
    """Open/index each selected source once; each query retains its actual scope."""
    if not source.selected:
        return [dict(source_id=source.source_id, state="unavailable", reason=source.selection_reason)
                for _ in variants], []
    if file_sha256(source.path) != source.sha256:
        raise ValueError("Normal source digest mismatch: " + source.source_id)
    if file_sha256(source.index_path) != source.index_sha256:
        raise ValueError("Normal source index digest mismatch: " + source.source_id)
    sample_id = source.vcf_sample
    if sample_id is None:
        raise ValueError("Selected source has no normal VCF sample column")
    rows: list[NormalSourceAssessment] = []
    backgrounds: list[tuple[AlleleKey, dict[str, Any]]] = []
    acquired = GenomicIntervalIndex(source.acquired_regions) if source.acquired_regions is not None else None
    with pysam.VariantFile(source.path, index_filename=source.index_path) as reader:
        if source.vcf_sample not in reader.header.samples:
            raise ValueError("Declared normal VCF sample is absent: " + source.source_id)
        resolver = ContigResolver(reader.header.contigs, reference_name=source.reference_name)
        for variant in variants:
            resolution = resolver.resolve(variant_source_contig(variant))
            interval = normal_query_interval(variant)
            row: NormalSourceAssessment = dict(source_id=source.source_id, state="unavailable", interval=list(interval),
                                                contig_resolution=asdict(resolution))
            contig = resolution.resolved_name
            if contig is None:
                row["reason"] = "unresolved_chromosome"
                rows.append(row)
                continue
            if acquired is not None and not acquired.covers(contig, *interval):
                row["reason"] = "outside_acquired_regions"
                rows.append(row)
                continue
            try:
                records = list(reader.fetch(contig, *interval))
                assessed = [assess_normal_record(record, sample_id, variant, source.policy) for record in records]
                state = normal_source_state(assessed, interval,
                    reference_confidence=source.evidence_kind == "reference_confidence")
                row["state"], row["records"] = state, assessed
                left, right = max(0, interval[0] - background_flank), interval[1] + background_flank
                if background_flank:
                    background_records = list(reader.fetch(contig, left, right))
                else:
                    background_records = records
                # A derivative may contain additional physical records. Only the
                # explicitly acquired territory is eligible for background use.
                if acquired is not None:
                    background_records = [record for record in background_records if acquired.covers(
                        contig, record.start, record.start + len(record.ref))]
                backgrounds.extend(_background_records(background_records, source, variant, source.policy))
                row["background_query_interval"] = [left, right]
                row["background_interval_fully_acquired"] = acquired is None or acquired.covers(contig, left, right)
            except (OSError, ValueError) as error:
                row["state"], row["reason"], row["error"] = "unavailable", "normal_query_failed", str(error)
            rows.append(row)
    return rows, backgrounds


def assess_matched_normal(variants: Sequence[Variant], sources: Sequence[NormalSource], *,
                          sample_id: str, reference_name: str, background_flank: int = 0) -> MatchedNormalEvidence:
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
        Additional bases queried on each side for supplied background calls.
        This bounds attribution scope; it does not establish an entire normal proteome.

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
    if len({source.source_id for source in sources}) != len(sources):
        raise ValueError("Duplicate normal source identity")
    if any(source.reference_name != reference_name for source in sources if source.selected):
        raise ValueError("Selected normal sources have incompatible reference assemblies")
    if any(variant_reference_name(variant) != reference_name for variant in variants):
        raise ValueError("Focal variants and normal evidence have incompatible reference assemblies")
    keys = [literal_allele_key(v.contig, v.start, v.ref, v.alt) for v in variants]
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate focal allele in normal assessment")
    source_rows = []
    background: dict[AlleleKey, dict[tuple[str, str, int], dict[str, Any]]] = {}
    for source in sources:
        if not source.selected:
            # Exclusions apply to the whole run. Keep them once in the source
            # inventory, rather than copying them into every allele assessment.
            continue
        rows, observations = _source_assessments(source, variants, background_flank)
        source_rows.append(rows)
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
    return dict(schema="isovar.matched_normal.v1", isovar_version=__version__, sample_id=sample_id,
                reference_name=reference_name, background_flank=background_flank,
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
    variants = [Variant(row["contig"], row["start"], row["ref"], row["alt"], genome=reference) for row in nominations]
    result = assess_matched_normal(variants, sources, sample_id=manifest["sample_id"], reference_name=reference,
                                   background_flank=manifest.get("background_flank", 0))
    for allele, nomination in zip(result["alleles"], nominations):
        allele["nomination"] = nomination
    result["manifest"] = dict(path=str(manifest_path), sha256=file_sha256(manifest_path),
                               snapshot=manifest.get("snapshot"), variants=manifest["variants"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_json(output_path, result)
    return result
