"""Matched-normal genotypes and reference confidence, without absence inference.

This layer assesses already retrieved VCF records. Acquisition, specimen
assignment and caller-product selection are explicit source-manifest decisions.
"""

from collections.abc import Iterable
from dataclasses import dataclass
from math import isfinite
from typing import Literal, TypedDict

from pysam import VariantRecord
from varcode import Variant, normalize_contig_name

from .contig_resolution import variant_reference_name


NormalState = Literal["alternate_observed", "conflicting_calls", "confidently_reference",
                      "other_allele_observed", "uncalled", "variant_only_absence", "unavailable"]


@dataclass(frozen=True)
class NormalEvidencePolicy:
    """Declared genotype gates; thresholds are analysis policy, not biological facts.

    Reference blocks require MIN_DP, because block DP can be an average. The
    chosen quality field is explicit: callers may use GQ or calibrated GQX.
    Unfiltered records (FILTER='.') are not accepted as PASS by default.
    """

    minimum_genotype_quality: int = 20
    minimum_depth: int = 5
    genotype_quality_field: str = "GQ"
    allow_unfiltered: bool = False

    def __post_init__(self) -> None:
        for value in (self.minimum_genotype_quality, self.minimum_depth):
            if type(value) is not int or value < 0:
                raise ValueError("Normal confidence thresholds must be nonnegative integers")
        if type(self.allow_unfiltered) is not bool:
            raise ValueError("allow_unfiltered must be a boolean")
        if self.genotype_quality_field not in ("GQ", "GQX"):
            raise ValueError("Normal genotype quality field must be GQ or GQX")


class NormalRecordEvidence(TypedDict):
    """Original record provenance and the conclusion permitted by its genotype."""

    contig: str
    record_id: str | None
    record_quality: float | None
    start: int
    end: int
    reference: str
    alternates: list[str]
    genotype: list[int | None]
    phased: bool
    phase_set: int | str | None
    allele_depths: list[int | None] | None
    filters: list[str]
    sample_filter: list[str]
    quality_field: str
    genotype_quality: float | None
    depth_field: str
    depth: int | None
    rejection_reasons: list[str]
    called_alternate_indexes: list[int]
    target_allele_called: bool
    accepted: bool
    reference_interval: list[int] | None


def literal_allele_key(contig: str, start: int, reference: str, alternate: str, *,
                       reference_name: str | None = None) -> tuple[str, int, str, str]:
    """Normalize spelling and trim with Varcode; do not left-align or lift over."""
    value = Variant(normalize_contig_name(contig, reference_name=reference_name), start, reference, alternate)
    return value.contig, value.start, value.ref, value.alt


def assess_normal_record(record: VariantRecord, sample_id: str, target: Variant,
                         policy: NormalEvidencePolicy) -> NormalRecordEvidence:
    """Assess one normal sample's exact genotype and retain every rejection reason.

    Parameters
    ----------
    record : pysam.VariantRecord
        A record on the source contig already resolved to the target's context.
    sample_id : str
        Exact VCF column name from the reviewed source descriptor.
    target : varcode.Variant
        Focal literal allele. Reference build equivalence is a caller obligation.
    policy : NormalEvidencePolicy
        Required filter, quality and depth gates.

    Returns
    -------
    NormalRecordEvidence
        Genotype-indexed ALT evidence and an optional zero-based half-open
        confidently reference interval. Missing or partial genotypes never
        establish reference. All ALT entries are retained, but only called
        allele indexes can supply background variants.
    """
    if sample_id not in record.samples:
        raise ValueError(f"Normal VCF sample column is missing: {sample_id}")
    sample = record.samples[sample_id]
    genotype = list(sample.get("GT") or ())
    alternates = list(record.alts or ())
    filters = list(record.filter)
    raw_filter = sample.get("FT")
    sample_filter = ([raw_filter] if isinstance(raw_filter, str) else list(raw_filter or ()))
    raw_quality = sample.get(policy.genotype_quality_field)
    quality = float(raw_quality) if isinstance(raw_quality, (int, float)) and isfinite(raw_quality) else None
    is_block = record.stop > record.start + len(record.ref)
    depth_field = "MIN_DP" if is_block else "DP"
    raw_depth = sample.get(depth_field)
    depth = raw_depth if isinstance(raw_depth, int) else None
    reasons = []
    if not genotype or any(index is None for index in genotype):
        reasons.append("missing_or_partial_genotype")
    if any(index is not None and not 0 <= index <= len(alternates) for index in genotype):
        reasons.append("invalid_genotype_allele_index")
    if set(filters) - {"PASS"} or (not filters and not policy.allow_unfiltered):
        reasons.append("record_not_pass")
    if any(value not in (None, ".", "PASS") for value in sample_filter):
        reasons.append("sample_filter_failed")
    if quality is None or quality < policy.minimum_genotype_quality:
        reasons.append("genotype_quality_unavailable_or_low")
    if depth is None or depth < policy.minimum_depth:
        reasons.append("minimum_depth_unavailable_or_low")
    target_start = target.original_start - 1
    left = max(record.start, target_start)
    right = min(record.start + len(record.ref), target_start + len(target.original_ref))
    if right > left and record.ref[left - record.start:right - record.start] != target.original_ref[
            left - target_start:right - target_start]:
        reasons.append("reference_allele_mismatch")
    indexes = sorted({index for index in genotype if index is not None and 0 < index <= len(alternates)})
    assembly = variant_reference_name(target)
    key = literal_allele_key(target.contig, target.start, target.ref, target.alt, reference_name=assembly)
    target_called = any(
        set(record.ref + alternates[index - 1]) <= set("ACGT") and
        literal_allele_key(target.contig, record.pos, record.ref, alternates[index - 1], reference_name=assembly) == key
        for index in indexes)
    all_reference = bool(genotype) and all(index == 0 for index in genotype)
    base_confidence_record = all(alt in ("<NON_REF>", "<*>") or set(record.ref + alt) <= set("ACGT")
                                 for alt in alternates)
    # A 0/0 structural-variant genotype does not establish base-level reference
    # confidence over its END interval. Only non-variant block symbols qualify.
    if is_block and any(alt not in ("<NON_REF>", "<*>") for alt in alternates):
        base_confidence_record = False
    return dict(contig=record.contig, record_id=record.id,
        record_quality=record.qual if record.qual is not None and isfinite(record.qual) else None,
        start=record.start, end=record.stop,
        reference=record.ref, alternates=alternates, genotype=genotype, phased=sample.phased,
        phase_set=sample.get("PS"), allele_depths=list(sample["AD"]) if sample.get("AD") is not None else None,
        filters=filters, sample_filter=sample_filter, quality_field=policy.genotype_quality_field,
        genotype_quality=quality, depth_field=depth_field, depth=depth, rejection_reasons=reasons,
        called_alternate_indexes=indexes, target_allele_called=target_called, accepted=not reasons,
        reference_interval=[record.start, record.stop]
            if all_reference and not reasons and base_confidence_record else None)


def interval_is_covered(interval: tuple[int, int], pieces: Iterable[tuple[int, int]]) -> bool:
    """Require every base of one half-open interval, allowing adjacent blocks."""
    start, end = interval
    if not 0 <= start < end:
        raise ValueError("Normal query interval must be nonempty and nonnegative")
    covered_until = start
    for left, right in sorted(pieces):
        if left > covered_until:
            break
        covered_until = max(covered_until, right)
        if covered_until >= end:
            return True
    return False


def normal_source_state(records: list[NormalRecordEvidence], interval: tuple[int, int], *,
                        reference_confidence: bool) -> NormalState:
    """Summarize one caller product; retain disagreement instead of vote counting."""
    alternate = any(row["accepted"] and row["target_allele_called"] for row in records)
    reference = interval_is_covered(interval, ((row["reference_interval"][0], row["reference_interval"][1]) for row in records
                                               if row["reference_interval"] is not None))
    other_alternate = any(row["accepted"] and row["called_alternate_indexes"] for row in records)
    if other_alternate and reference:
        return "conflicting_calls"
    if alternate:
        return "alternate_observed"
    if reference:
        return "confidently_reference"
    if other_alternate:
        return "other_allele_observed"
    return "uncalled" if records or reference_confidence else "variant_only_absence"


def reconcile_normal_states(states: Iterable[NormalState]) -> NormalState:
    """Reconcile caller products without treating processing copies as replication.

    Positive ALT or reference calls can coexist with missing evidence elsewhere;
    both kinds of positive call together remain conflicting. Per-source outcomes
    must accompany this summary so unavailable territory is never erased.
    """
    values = set(states)
    if "conflicting_calls" in values or ("confidently_reference" in values and
                                         values.intersection(("alternate_observed", "other_allele_observed"))):
        return "conflicting_calls"
    for state in ("alternate_observed", "confidently_reference", "other_allele_observed",
                  "uncalled", "variant_only_absence"):
        if state in values:
            return state
    return "unavailable"
