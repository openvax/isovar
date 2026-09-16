"""Transfer local coding phase through observed, transcript-compatible alignments."""

from .dna import reverse_complement_dna
from .variant_helpers import base0_interval_for_variant


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
    mappings = {}
    retained_start = len(variant_sequence) - len(prefix) if reverse else 0
    retained_end = len(variant_sequence) if reverse else len(prefix)
    for read in variant_sequence.reads:
        shift = len(variant_sequence.prefix) - len(read.prefix)
        blocks = read.reference_blocks
        for index, (q0, q1, g0, _) in enumerate(blocks):
            # Internal query gaps are insertions, not unobserved positions.
            start = blocks[index - 1][1] if index else q0
            for query in range(max(start, retained_start - shift), min(q1, retained_end - shift)):
                q = shift + query
                if reverse:
                    q = len(variant_sequence) - 1 - q
                if not 0 <= q < len(prefix):
                    continue
                g = g0 + query - q0 if query >= q0 else None
                if q in mappings and mappings[q] != g:
                    return None
                mappings[q] = g

    g_start, g_end = base0_interval_for_variant(reference_context.variant)
    results = set()
    for transcript in reference_context.transcripts:
        exons = sorted((a - 1, b) for a, b in transcript.exon_intervals)
        transcript_length = sum(b - a for a, b in exons)

        def spliced_boundary(g):
            # Interbase coordinate; also valid at exon edges and insertions.
            before = sum(max(0, min(g, b) - a) for a, b in exons)
            return transcript_length - before if reverse else before

        locus = spliced_boundary(g_end if reverse else g_start)
        reference_start = locus - len(reference)
        pairs = {}
        for q, g in mappings.items():
            if g is None or not any(a <= g < b for a, b in exons):
                continue
            r = spliced_boundary(g + 1 if reverse else g) - reference_start
            if 0 <= r < len(reference):
                pairs[q] = r
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
