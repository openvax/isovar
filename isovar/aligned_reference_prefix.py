"""Transfer local coding phase through observed, transcript-compatible alignments."""

from bisect import bisect_right
from itertools import accumulate

from .dna import reverse_complement_dna
from .variant_helpers import base0_interval_for_variant


class SplicedCoordinates:
    """Strand-oriented transcript offsets of zero-based genomic coordinates.

    Parameters
    ----------
    transcript : object
        Supplies ``exon_intervals``: disjoint, one-based inclusive exons.
    reverse : bool
        Whether the transcript is on the negative strand.
    """

    def __init__(self, transcript, reverse):
        self.exons = sorted((a - 1, b) for a, b in transcript.exon_intervals)
        self.starts = [a for a, _ in self.exons]
        # Exonic bases before each exon; the last entry is the transcript length.
        self.before = list(accumulate((b - a for a, b in self.exons), initial=0))
        self.reverse = reverse

    def boundary(self, g):
        """Transcript offset of an interbase genomic coordinate.

        Parameters
        ----------
        g : int
            Zero-based interbase coordinate; exon edges, introns and
            insertion points are valid.

        Returns
        -------
        int
            Number of transcript bases upstream of the coordinate.
        """
        i = bisect_right(self.starts, g) - 1
        before = self.before[i] + min(g, self.exons[i][1]) - self.exons[i][0] if i >= 0 else 0
        return self.before[-1] - before if self.reverse else before

    def base(self, g):
        """Transcript offset of one genomic base.

        Parameters
        ----------
        g : int
            Zero-based genomic position.

        Returns
        -------
        int or None
            Offset of an exonic base; None for an intronic or flanking one.
        """
        i = bisect_right(self.starts, g) - 1
        if i < 0 or g >= self.exons[i][1]:
            return None
        before = self.before[i] + g - self.exons[i][0]
        return self.before[-1] - 1 - before if self.reverse else before


def aligned_reference_prefix(variant_sequence, reference_context):
    """Return RNA/reference prefixes, first-codon offset and edit distance.

    Only an observed aligned base can anchor the frame. Insertions consume RNA
    positions; deletions consume reference positions; annotated introns consume
    neither transcript positions nor phase (SAM v1.6 CIGAR semantics).

    Returns None for contradictory or missing alignment evidence. Callers must
    use the sequence-only compatibility path *only* when no reads have mappings.
    All offsets here are zero-based; genomic and query intervals are half-open.
    """
    reverse = reference_context.strand == "-"
    prefix = (reverse_complement_dna(variant_sequence.suffix) if reverse
              else variant_sequence.prefix)
    reference = reference_context.sequence_before_variant_locus
    # Read geometry belongs to the sequence, so every coding model reuses it.
    placements = variant_sequence.aligned_positions()
    last = len(variant_sequence) - 1
    mappings = {}
    for q in range(len(prefix)):
        observed = placements[last - q if reverse else q]
        if len(observed) > 1:
            return None
        if observed:
            mappings[q], = observed

    g_start, g_end = base0_interval_for_variant(reference_context.variant)
    results = set()
    for transcript in reference_context.transcripts:
        spliced = SplicedCoordinates(transcript, reverse)
        reference_start = spliced.boundary(g_end if reverse else g_start) - len(reference)
        pairs = {}
        for q, g in mappings.items():
            r = None if g is None else spliced.base(g)
            if r is not None and 0 <= r - reference_start < len(reference):
                pairs[q] = r - reference_start
        if not pairs:
            return None
        q_start = min(pairs)
        r_start = pairs[q_start]
        # No sequence-only bridge can establish an alignment/frame.
        if any(q not in mappings for q in range(q_start, len(prefix))):
            return None
        ordered = sorted(pairs.items())
        if any(r1 <= r0 for (_, r0), (_, r1) in zip(ordered, ordered[1:])):
            return None
        mismatches = sum(prefix[q] != reference[r] for q, r in ordered)
        mismatches += len(prefix) - q_start - len(pairs)
        mismatches += len(reference) - r_start - len(pairs)

        codon = reference_context.offset_to_first_complete_codon
        if reference_context.contains_start_codon and r_start <= codon:
            # A deleted/disrupted annotated start is not a justified new start.
            start_queries = [q for q, r in ordered if codon <= r < codon + 3]
            if (len(start_queries) != 3
                    or start_queries != list(range(start_queries[0], start_queries[0] + 3))
                    or prefix[start_queries[0]:start_queries[0] + 3] != reference[codon:codon + 3]):
                return None
            offset = start_queries[0] - q_start
        else:
            offset = (codon - r_start) % 3
        results.add((prefix[q_start:], reference[r_start:], offset, mismatches))
    # A context may group several transcripts; never select one silently.
    return next(iter(results)) if len(results) == 1 else None
