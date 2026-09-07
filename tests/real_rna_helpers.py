"""Independent CIGAR observations for the public RNA regression fixtures.

Deliberately does not import isovar or perform repeat/indel normalization.
Coordinates in the fixture manifest are original, anchored, 1-based VCF
coordinates. Query offsets are in SAM orientation, including soft clips.
"""

from hashlib import sha256


def record_digest(read):
    return sha256(read.to_string().encode("ascii")).hexdigest()


def cigar_observation(read, variant):
    """Return the directly aligned allele and its query quality indices.

    A SNP needs its own aligned base. An indel needs both reference anchors;
    N (reference skip) must never become a deletion. Empty observations use
    anchor qualities, which are NOT a probability of the indel being correct.
    """
    pos, ref, alt = variant["pos"], variant["ref"], variant["alt"]
    if read.is_unmapped or read.query_sequence is None or read.reference_name != variant["chrom"]:
        return None
    reference_to_query = {}
    skips = []
    query_pos, reference_pos = 0, read.reference_start
    for operation, length in read.cigartuples:
        if operation in (0, 7, 8):  # M, =, X
            reference_to_query.update(
                (reference_pos + i, query_pos + i) for i in range(length))
        elif operation == 3:  # N, NOT D
            skips.append((reference_pos, reference_pos + length))
        if operation in (0, 1, 4, 7, 8):
            query_pos += length
        if operation in (0, 2, 3, 7, 8):
            reference_pos += length
    if len(ref) == len(alt) == 1:
        q = reference_to_query.get(pos - 1)
        if q is None:
            left = reference_to_query.get(pos - 2)
            right = reference_to_query.get(pos)
            if left is None or right is None or any(s <= pos - 1 < e for s, e in skips):
                return None
            start, end = left + 1, right
            indices = list(range(start, end)) if start < end else [left, right]
        else:
            # A trailing insertion is an additional change, not pure SNP
            # support. Keep it in the allele, with the base-only observation
            # separately available for a pileup-style quality audit.
            start, end = q, reference_to_query.get(pos, q + 1)
            indices = list(range(start, end))
    else:
        assert ref[0] == alt[0] and min(len(ref), len(alt)) == 1
        left = reference_to_query.get(pos - 1)
        right = reference_to_query.get(pos + len(ref) - 1)
        if left is None or right is None:
            return None
        if any(start < pos + len(ref) and end > pos - 1 for start, end in skips):
            return None
        start, end = left + 1, right
        indices = list(range(start, end)) if start < end else [left, right]
    qualities = read.query_qualities
    return {
        "record": record_digest(read),
        "name": read.query_name,
        "flag": read.flag,
        "allele": read.query_sequence[start:end],
        "base": read.query_sequence[q] if len(ref) == len(alt) == 1 and q is not None else None,
        "base_quality": (read.query_qualities[q] if len(ref) == len(alt) == 1
                         and q is not None and read.query_qualities is not None else None),
        "query_interval": [start, end],
        "quality_indices": indices,
        "qualities": None if qualities is None else [qualities[i] for i in indices],
        "mapq": read.mapping_quality,
        "nh": read.get_tag("NH") if read.has_tag("NH") else None,
    }
