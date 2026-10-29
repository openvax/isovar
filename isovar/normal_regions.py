"""Explicit background territory shared by normal acquisition and assessment."""

from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass
from typing import Any

from pysam import VariantHeader
from varcode import ContigResolver, Variant

from .contig_resolution import variant_source_contig
from .evidence_report import changed_reference_interval
from .genomic_intervals import GenomicIntervalIndex


GenomicRegion = tuple[str, int, int]


def normal_query_interval(variant: Variant) -> tuple[int, int]:
    """Assess changed bases, or both insertion flanks, in zero-based coordinates."""
    start, end, _ = changed_reference_interval(variant.original_start, variant.original_ref, variant.original_alt)
    return start - 1, end


def normal_background_regions(variants: Sequence[Variant], *, background_flank: int = 0,
                              additional_regions: Iterable[GenomicRegion] = ()) -> tuple[GenomicRegion, ...]:
    """Merge focal neighborhoods and explicit additional background territory.

    Parameters
    ----------
    variants : sequence of varcode.Variant
        Focal alleles whose changed intervals are included in the territory.
    background_flank : int, optional
        Nonnegative additional bases on each side of every focal interval.
    additional_regions : iterable of (str, int, int), optional
        Contig, zero-based start and half-open end in the variants' assembly.
        Include distant exons or ncORF blocks explicitly; flanking a genomic
        locus does not cover an entire spliced transcript.

    Returns
    -------
    tuple of (str, int, int)
        Sorted, merged intervals. Names remain as supplied until resolved against
        each source header; no acquisition or reference confidence is implied.
    """
    if type(background_flank) is not int or background_flank < 0:
        raise ValueError("Background flank must be a nonnegative integer")
    regions = list(additional_regions)
    for variant in variants:
        start, end = normal_query_interval(variant)
        regions.append((str(variant_source_contig(variant)), max(0, start - background_flank), end + background_flank))
    return GenomicIntervalIndex(regions).intervals()


@dataclass(frozen=True)
class ResolvedNormalRegions:
    """Source-header territory and explicit failures to acquire requested regions."""

    regions: tuple[GenomicRegion, ...]
    unresolved_contigs: tuple[dict[str, Any], ...]
    outside_contig_regions: tuple[GenomicRegion, ...]
    clipped_regions: tuple[GenomicRegion, ...]


def resolve_normal_regions(regions: Iterable[GenomicRegion], header: VariantHeader, *,
                           reference_name: str) -> ResolvedNormalRegions:
    """Resolve once per distinct name; bound intervals by declared contig lengths.

    Missing and ambiguous names remain unavailable. Exact source names retain
    priority, including viral, decoy and alternate contigs; no liftover occurs.
    """
    resolver = ContigResolver((str(name) for name in header.contigs), reference_name=reference_name)
    resolutions = {}
    resolved, outside, clipped = [], [], []
    for name, start, end in GenomicIntervalIndex(regions).intervals():
        if name not in resolutions:
            resolutions[name] = resolver.resolve(name)
        contig = resolutions[name].resolved_name
        if contig is None:
            continue
        length = header.contigs[contig].length
        if length is not None and start >= length:
            outside.append((name, start, end))
            continue
        if length is not None and end > length:
            clipped.append((name, start, end))
            end = length
        resolved.append((contig, start, end))
    return ResolvedNormalRegions(GenomicIntervalIndex(resolved).intervals(),
        tuple(asdict(value) for _, value in sorted(resolutions.items()) if value.resolved_name is None),
        tuple(outside), tuple(clipped))
