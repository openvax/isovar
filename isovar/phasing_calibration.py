"""Conservative directional phasing-error bounds from truth-labeled RNA controls.

The controls must be independent of the sample being phased and must have a
known reference or alternate allele at every calibrated locus. This module
does not infer that truth from the observed reads.
"""

from collections import defaultdict
from dataclasses import dataclass
from math import exp, expm1, fsum, lgamma, log, log1p

from .read_identity import fragment_ids, source_read_ids
from .variant_helpers import trim_variant


@dataclass(frozen=True)
class PhasingControl:
    """One known-allele control's collected evidence at one locus and read group."""

    variant: object
    truth: str
    read_group: str
    source: str
    platform: str
    basecaller: str
    library: str
    biological_context: str
    reference: str
    context: str
    evidence: object


@dataclass
class PhasingCalibration:
    """Rates accepted by ``run_isovar`` and a provenance-rich audit report."""

    rates: dict
    report: list
    confidence: float

    def profile_rows(self):
        """Rows accepted by the existing ``--phasing-error-rates`` option."""
        return [dict(contig=v.contig, start=v.start, ref=v.ref, alt=v.alt,
                     ref_to_alt=rates[0], alt_to_ref=rates[1])
                for v, rates in sorted(self.rates.items(), key=lambda item: (
                    item[0].contig, item[0].start, item[0].ref, item[0].alt))]


def _binomial_cdf(errors, trials, rate):
    """P(X <= errors), computed in log space to retain low-error tails."""
    if rate <= 0:
        return 1.0
    if rate >= 1:
        return 0.0
    scale = lgamma(trials + 1)
    return min(1.0, fsum(exp(
        scale - lgamma(k + 1) - lgamma(trials - k + 1)
        + k * log(rate) + (trials - k) * log1p(-rate))
        for k in range(errors + 1)))


def one_sided_error_bound(errors, trials, confidence=0.95):
    """Exact Clopper-Pearson upper bound for a binomial wrong-allele rate."""
    if not isinstance(trials, int) or not isinstance(errors, int) or not 0 <= errors <= trials or trials == 0:
        raise ValueError("Need integer wrong-allele and callable-fragment counts, with 0 <= errors <= trials and trials > 0")
    if not 0 < confidence < 1:
        raise ValueError("confidence must be between 0 and 1")
    if errors == trials:
        return 1.0
    if errors == 0:
        return -expm1(log(1 - confidence) / trials)
    low, high = 0.0, 1.0
    alpha = 1 - confidence
    for _ in range(60):
        mid = (low + high) / 2
        if _binomial_cdf(errors, trials, mid) > alpha:
            low = mid
        else:
            high = mid
    return high


def _fragment_counts(control):
    evidence = control.evidence
    if (evidence.trimmed_base1_start, evidence.trimmed_ref, evidence.trimmed_alt) != trim_variant(control.variant):
        raise ValueError("Control evidence does not match its variant")
    buckets = []
    for allele in ("ref", "alt", "other"):
        reads = getattr(evidence, allele + "_reads")
        for read in reads:
            keys = source_read_ids(read)
            if not keys or any(key[0] != control.read_group for key in keys):
                raise ValueError("Control reads need source identities in the declared read group")
        buckets.append(fragment_ids(reads))
    ambiguous = (buckets[0] & buckets[1]) | (buckets[0] & buckets[2]) | (buckets[1] & buckets[2])
    ref, alt, other = (bucket - ambiguous for bucket in buckets)
    return len(ref), len(alt), len(other), len(ambiguous)


def calibrate_phasing_error_rates(controls, target_metadata, confidence=0.95):
    """Estimate directional rates from independently truth-labeled controls.

    ``target_metadata`` maps each variant to its included read groups. Each
    group maps to ``platform``, ``basecaller``, ``library``,
    ``biological_context``, ``reference`` and local sequence ``context``
    strings. Both true-ref and true-alt controls must match every group. A
    Bonferroni correction covers both directions, groups and control sources
    at each locus. The largest bound becomes the locus rate. Missing or fully
    erroneous controls fail instead of falling back to the scalar floor.
    """
    if not 0 < confidence < 1:
        raise ValueError("confidence must be between 0 and 1")
    if not target_metadata:
        raise ValueError("At least one target locus is required")
    totals = {}
    sources = defaultdict(list)
    seen = set()
    source_truth = {}
    for control in controls:
        if control.truth not in ("ref", "alt") or not control.source or not control.read_group:
            raise ValueError("Controls need known ref/alt truth, source and read group")
        if control.variant not in target_metadata or control.read_group not in target_metadata[control.variant]:
            raise ValueError("Control locus or read group is absent from the target")
        metadata = target_metadata[control.variant][control.read_group]
        metadata_keys = ("platform", "basecaller", "library", "biological_context",
                         "reference", "context")
        expected = {key: getattr(control, key) for key in metadata_keys}
        if any(not value or metadata.get(key) != value for key, value in expected.items()):
            raise ValueError("Control metadata differs from target")
        key = (control.variant, control.read_group, control.truth)
        unique = key + (control.source,)
        if unique in seen:
            raise ValueError("Duplicate control source at a locus and read group")
        seen.add(unique)
        source_key = (control.variant, control.read_group, control.source)
        if source_key in source_truth and source_truth[source_key] != control.truth:
            raise ValueError("One control source cannot have both ref and alt truth")
        source_truth[source_key] = control.truth
        totals[unique] = _fragment_counts(control)
        sources[key].append(control.source)

    rates, report = {}, []
    for variant, groups in target_metadata.items():
        if not groups:
            raise ValueError("Each target locus needs read groups")
        strata = [(variant, group, truth) for group in groups for truth in ("ref", "alt")]
        missing = [key for key in strata if not sources[key]]
        if missing:
            _, group, truth = missing[0]
            raise ValueError("Missing %s control for %s in read group %s" % (truth, variant, group))
        adjusted_confidence = 1 - (1 - confidence) / sum(len(sources[key]) for key in strata)
        group_rates = []
        for group, metadata in groups.items():
            if not group or any(not metadata.get(key) for key in (
                    "platform", "basecaller", "library", "biological_context",
                    "reference", "context")):
                raise ValueError("Target read groups need complete calibration metadata")
            directions = []
            for truth in ("ref", "alt"):
                key = (variant, group, truth)
                source_bounds = []
                for source in sorted(sources[key]):
                    ref, alt, other, ambiguous = totals[key + (source,)]
                    correct, wrong = (ref, alt) if truth == "ref" else (alt, ref)
                    callable_count = correct + wrong
                    if not callable_count:
                        raise ValueError("No callable control fragments for %s in read group %s" % (variant, group))
                    bound = one_sided_error_bound(wrong, callable_count, adjusted_confidence)
                    if bound >= 1:
                        raise ValueError("Control error bound is 1; phasing cannot be calibrated")
                    source_bounds.append(bound)
                    report.append(dict(contig=variant.contig, start=variant.start,
                                       ref=variant.ref, alt=variant.alt, read_group=group,
                                       truth=truth, source=source,
                                       platform=metadata["platform"], basecaller=metadata["basecaller"],
                                       library=metadata["library"],
                                       biological_context=metadata["biological_context"],
                                       reference=metadata["reference"], context=metadata["context"],
                                       correct=correct, wrong=wrong, other=other,
                                       ambiguous=ambiguous, callable=callable_count,
                                       upper_bound=bound, confidence=adjusted_confidence))
                directions.append(max(source_bounds))
            group_rates.append(directions)
        rates[variant] = tuple(max(group[i] for group in group_rates) for i in (0, 1))
    return PhasingCalibration(rates, report, confidence)
