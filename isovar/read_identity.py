"""Sequenced segments are observations; SAM placements are hypotheses.

Collected reads carry immutable ``(segment_id, alignment_id)`` pairs.
Segment IDs are (read group, QNAME, first/last-segment flag bits); alignment
IDs are (reference ID, reference start, CIGAR, reverse strand). Legacy
caller-created reads without this metadata retain their original semantics.
"""


def source_alignments_from_pysam(read, name):
    group = read.get_tag("RG") if read.has_tag("RG") else ""
    segment = read.flag & 0xC0 if read.is_paired else 0
    return (((group, name, segment), (
        read.reference_id, read.reference_start, read.cigarstring,
        read.is_reverse)),)


def source_read_ids(read):
    return {key for key, _ in getattr(read, "source_alignments", ())}


def count_reads(reads):
    """Count each sequenced segment once, including merged mates."""
    identities = set()
    legacy_count = 0
    for read in reads:
        keys = source_read_ids(read)
        if keys:
            identities.update(keys)
        else:
            legacy_count += getattr(read, "source_read_count", 1)
    return len(identities) + legacy_count


def observation_groups(reads):
    """Deduplicate placements without changing merged-pair coverage semantics.

    A collapsed pair contributes one coverage observation, while two separate
    mates contribute two. A retained single-mate view cannot count again when
    its exact placement is already represented by a collapsed pair.
    """
    reads = list(reads)
    merged = {}
    for read in reads:
        keys = frozenset(source_read_ids(read))
        if len(keys) > 1:
            merged.update((key, keys) for key in keys)
    groups = {}
    for read in reads:
        keys = source_read_ids(read)
        if keys:
            key = next(iter(keys))
            key = merged.get(key, frozenset([key]))
        else:
            key = read
        groups.setdefault(key, []).append(read)
    return list(groups.values())


def fragment_ids(reads):
    """Read-group-aware fragment keys, with raw names for legacy objects."""
    result = set()
    for read in reads:
        keys = source_read_ids(read)
        if keys:
            result.update(key[:2] for key in keys)
        else:
            result.add(getattr(read, "name", read))
    return result


def alignment_constraints(reads):
    """Map segments to their selected placements, rejecting self-extension."""
    result = {}
    for read in reads:
        for key, placement in getattr(read, "source_alignments", ()):
            if key in result and result[key] != placement:
                raise ValueError("Alternative alignments of one segment cannot support one assembly")
            result[key] = placement
    return result


def compatible_alignments(constraints, reads):
    """Whether adding these observations selects at most one placement/segment."""
    combined = dict(constraints)
    for read in reads:
        for key, placement in getattr(read, "source_alignments", ()):
            if key in combined and combined[key] != placement:
                return False
            combined[key] = placement
    return True


def partition_read_alignments(reads):
    """Keep incompatible placements separate when grouping identical cDNA."""
    groups = []
    for read in sorted(reads, key=read_sort_key):
        placements = getattr(read, "source_alignments", ())
        for constraints, group in groups:
            if all(key not in constraints or constraints[key] == value
                   for key, value in placements):
                constraints.update(placements)
                group.add(read)
                break
        else:
            groups.append((dict(placements), {read}))
    return [group for _, group in groups]


def read_sort_key(read):
    """Stable ordering for alternative observations (also accepts legacy mocks)."""
    return (tuple(getattr(read, "source_alignments", ())),
            getattr(read, "name", read if isinstance(read, str) else ""),
            getattr(read, "prefix", ""), getattr(read, "allele", ""),
            getattr(read, "suffix", ""))
