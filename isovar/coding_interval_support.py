"""Original RNA coverage of an explicitly positioned coding interval."""

from collections.abc import Callable, Hashable, Iterable, Sequence
from typing import TypedDict

from .allele_read import AlleleRead
from .dna import reverse_complement_dna
from .read_identity import fragment_ids
from .rna_evidence import RNASupport


class CodingIntervalSupport(TypedDict):
    scope: str
    coding_interval: list[int]
    fragment_depth: list[int]
    minimum_fragment_depth: int
    fully_tiled: bool
    allele_observation_span_support: RNASupport
    single_read_span_support: RNASupport
    paired_template_span_support: RNASupport


def coding_interval_support(
        reads: Iterable[AlleleRead], sequence: str, variant_start: int,
        interval: Sequence[int], strand: str,
        support: Callable[[list[AlleleRead]], RNASupport]) -> CodingIntervalSupport:
    """Measure exact allele-anchored coverage without counting partial contributors as spans.

    Parameters
    ----------
    reads : iterable of AlleleRead
        Original allele observations, possibly merged overlapping mates.
    sequence : str
        Candidate cDNA in transcript orientation.
    variant_start : int
        Zero-based start of the focal alternate allele in ``sequence``.
    interval : pair of int
        Half-open coding interval, including an observed stop when requested.
    strand : str
        Transcript strand relative to the observations' genomic orientation.
    support : callable
        Converts a read collection to its scoped RNA support record.

    Returns
    -------
    dict
        Per-base fragment depth, exact full-interval allele-observation support
        and its single-read and merged paired-template subsets. The latter
        counts observations made from two overlapping mates, not unsequenced
        inserts. Single-read spans exclude merged observations; they do not
        recount the original mates hidden within a merged consensus. Different placements
        and overlapping observations of a fragment count once at each base.
    """
    start, end = interval
    if strand not in ("+", "-") or not 0 <= start <= end <= len(sequence):
        raise ValueError("Invalid coding interval or strand")
    if not 0 <= variant_start <= len(sequence):
        raise ValueError("Invalid focal allele offset")
    covered: list[set[Hashable]] = [set() for _ in range(end - start)]
    spanning, single, paired = [], [], []
    for read in reads:
        observed = read.sequence
        allele_start = len(read.prefix)
        if strand == "-":
            observed = reverse_complement_dna(observed)
            allele_start = len(read.suffix)
        shift = variant_start - allele_start
        left, right = max(start, shift), min(end, shift + len(observed))
        if right <= left or sequence[left:right] != observed[left - shift:right - shift]:
            continue
        identities = fragment_ids([read])
        for position in range(left, right):
            covered[position - start].update(identities)
        if left == start and right == end:
            spanning.append(read)
            if read.source_read_count == 1:
                single.append(read)
            elif read.source_read_count == 2:
                paired.append(read)
    depths = [len(ids) for ids in covered]
    return dict(
        scope="exact_observed_coding_interval",
        coding_interval=[start, end], fragment_depth=depths,
        minimum_fragment_depth=min(depths) if depths else 0,
        fully_tiled=bool(depths) and all(depths),
        allele_observation_span_support=support(spanning),
        single_read_span_support=support(single),
        paired_template_span_support=support(paired))
