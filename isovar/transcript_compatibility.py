"""Classify aligned RNA observations without choosing one transcript."""


def _base0_exons(transcript):
    """Return transcript exons as sorted 0-based half-open intervals."""
    return tuple(sorted((start - 1, end) for start, end in transcript.exon_intervals))


def _local_path(read, max_prefix_size, max_suffix_size):
    """Clip compact alignment blocks to the sequence entering assembly."""
    query_start = 0
    query_end = len(read.sequence)
    if max_prefix_size is not None:
        query_start = max(0, len(read.prefix) - max_prefix_size)
    if max_suffix_size is not None:
        query_end = min(
            query_end,
            len(read.prefix) + len(read.allele) + max_suffix_size)

    blocks = []
    for block_query_start, block_query_end, block_ref_start, block_ref_end in read.reference_blocks:
        start = max(query_start, block_query_start)
        end = min(query_end, block_query_end)
        if start < end:
            ref_start = block_ref_start + start - block_query_start
            blocks.append((ref_start, ref_start + end - start))
    block_gaps = {
        (left[1], right[0])
        for left, right in zip(blocks, blocks[1:])
    }
    junctions = frozenset(
        junction for junction in read.splice_junctions if junction in block_gaps)
    return tuple(blocks), junctions


def _compatible_with_path(blocks, observed_junctions, exons, transcript_junctions):
    """Whether every observed aligned block and splice agrees with a path."""
    if not exons:
        return False

    # Sequence beyond an annotated transcript end says nothing about its
    # internal splice path. Clip only at the outer boundaries: extension into
    # an intron between exons must still make the read incompatible.
    transcript_start, transcript_end = exons[0][0], exons[-1][1]
    blocks = tuple(
        (max(start, transcript_start), min(end, transcript_end))
        for start, end in blocks
        if start < transcript_end and end > transcript_start)
    observed_junctions = frozenset(
        junction
        for junction in observed_junctions
        if transcript_start <= junction[0] and junction[1] <= transcript_end)

    exon_indices = []
    exon_index = 0
    for block_start, block_end in blocks:
        while exon_index < len(exons) and exons[exon_index][1] <= block_start:
            exon_index += 1
        if exon_index == len(exons):
            return False
        exon_start, exon_end = exons[exon_index]
        if block_start < exon_start or block_end > exon_end:
            return False
        exon_indices.append(exon_index)

    if not observed_junctions.issubset(transcript_junctions):
        return False

    for i, (left, right) in enumerate(zip(blocks, blocks[1:])):
        if exon_indices[i] != exon_indices[i + 1] and (left[1], right[0]) not in observed_junctions:
            return False
    return True


def annotate_reads_with_transcript_compatibility(
        reads,
        transcripts,
        max_prefix_size=None,
        max_suffix_size=None):
    """Attach every compatible annotated transcript ID to each read.

    Reads without alignment-path metadata, such as manually created
    ``AlleleRead`` objects, are compatible with all supplied transcripts.
    BAM-derived reads are compatible only with transcripts sharing their path.
    """
    paths = []
    for transcript in transcripts:
        exons = _base0_exons(transcript)
        junctions = frozenset(
            (left[1], right[0])
            for left, right in zip(exons, exons[1:]))
        paths.append((transcript.id, exons, junctions))

    all_transcript_ids = frozenset(transcript_id for transcript_id, _, _ in paths)
    result = []
    for read in reads:
        blocks, observed_junctions = _local_path(
            read, max_prefix_size, max_suffix_size)
        if not blocks and not observed_junctions:
            compatible_ids = all_transcript_ids
        else:
            compatible_ids = frozenset(
                transcript_id
                for transcript_id, exons, junctions in paths
                if _compatible_with_path(
                    blocks, observed_junctions, exons, junctions))
        result.append(read.with_compatible_transcript_ids(compatible_ids))
    return result
