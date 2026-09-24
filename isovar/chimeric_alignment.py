"""SAM supplementary paths are pieces of one segment, not mates or extra reads.

Keep this evidence separate from ``source_alignments``: ordinary cDNA assembly
still selects one linear placement per segment. Cross-record phasing additionally
requires reciprocal SA links, agreement on the remaining path, compatible query
spans, and the same primary/secondary hypothesis class.
"""

from collections import namedtuple
import re

from .logging import get_logger
from .read_identity import compatible_alignments

logger = get_logger(__name__)

_Alignment = namedtuple(
    "_Alignment", "contig start reverse cigar query_start query_end query_length reference_end")
_AlignmentPath = namedtuple("_AlignmentPath", "secondary supplementary alignment others query_interval")
_CIGAR = re.compile(r"(?:[1-9][0-9]*[MIDNSHP=X])+")


def _alignment(contig, start, reverse, cigar):
    """Normalize clipping and =/X, retaining the actual interior alignment."""
    if not isinstance(cigar, str) or not _CIGAR.fullmatch(cigar) or start < 0:
        raise ValueError("invalid alignment coordinates/CIGAR")
    operations = []
    for length, op in re.findall(r"([0-9]+)([MIDNSHP=X])", cigar):
        op = "S" if op == "H" else "M" if op in "=X" else op
        length = int(length)
        if operations and operations[-1][1] == op:
            operations[-1] = (operations[-1][0] + length, op)
        else:
            operations.append((length, op))
    if any(op == "S" for _, op in operations[1:-1]):
        raise ValueError("internal clipping")
    left = operations[0][0] if operations[0][1] == "S" else 0
    right = operations[-1][0] if operations[-1][1] == "S" else 0
    aligned_query = sum(n for n, op in operations if op in "MI")
    aligned_reference = sum(n for n, op in operations if op in "MDN")
    if not aligned_query or not aligned_reference:
        raise ValueError("empty aligned span")
    length = left + aligned_query + right
    start_query = right if reverse else left
    return _Alignment(contig, start, reverse, tuple(operations), start_query,
                      start_query + aligned_query, length, start + aligned_reference)


def sa_entries(tag):
    """Parse a SAM ``SA`` tag into (contig, 0-based start, strand, CIGAR, MAPQ, NM).

    Raises ValueError for a malformed entry. See the SAMtags definition.
    """
    if not isinstance(tag, str):
        raise ValueError("SA tag is not a string")
    entries = []
    for entry in tag.rstrip(";").split(";"):
        contig, position, strand, cigar, mapq, nm = entry.split(",")
        entries.append((contig, int(position) - 1, strand, cigar, int(mapq), int(nm)))
    return entries


def source_alignment_paths_from_pysam(read, source_alignments, query_interval=None):
    """Retain a validated SA declaration; malformed tags never discard reads.

    Each immutable entry is ``(segment_id, path)``. The path retains SAM's
    secondary/supplementary distinction and its own/other alignment descriptors.
    No SA declaration creates an observation: both variant-supporting records
    must actually be collected and pass the usual read filters.
    """
    if not read.has_tag("SA") or read.is_unmapped or read.cigarstring is None:
        return ()
    try:
        own = _alignment(read.reference_name, read.reference_start,
                         read.is_reverse, read.cigarstring)
        others = []
        for contig, start, strand, cigar, mapq, nm in sa_entries(read.get_tag("SA")):
            if (strand not in ("+", "-") or read.header.get_tid(contig) < 0
                    or not 0 <= mapq <= 255 or nm < 0):
                raise ValueError("invalid SA entry")
            others.append(_alignment(contig, start, strand == "-", cigar))
        members = sorted([own] + others, key=lambda a: (a.query_start, a.query_end))
        if (any(a.query_length != own.query_length for a in others)
                or any(a.query_start >= b.query_start or a.query_end >= b.query_end
                       for a, b in zip(members, members[1:]))):
            # Junction microhomology may overlap; a contained/equal query span
            # is another placement of those bases, not an additional piece.
            raise ValueError("inconsistent or contained query spans")
        if query_interval is not None:
            hard_clip = 0
            for op, n in read.cigartuples:
                if op != 5:
                    break
                hard_clip += n
            start, end = (hard_clip + i for i in query_interval)
            query_interval = ((own.query_length - end, own.query_length - start)
                              if own.reverse else (start, end))
        path = _AlignmentPath(bool(read.is_secondary), bool(read.is_supplementary),
                              own, tuple(sorted(others)), query_interval)
        return ((source_alignments[0][0], path),)
    except ValueError as error:
        logger.warning("Ignoring chimeric-path evidence for read %s: %s", read.query_name, error)
        return ()


def _matches_reported(actual, reported):
    """Match exact SA CIGARs or minimap2's documented span-only SA summary.

    The summary must match both reference and original-query coordinates. It
    never replaces the actual CIGAR used for allele extraction or assembly.
    See minimap2 format.c (SA writer) and the SAMtags SA definition.
    """
    if actual[:3] != reported[:3] or actual[4:] != reported[4:]:
        return False
    if actual.cigar == reported.cigar:
        return True
    query = actual.query_end - actual.query_start
    reference = actual.reference_end - actual.start
    left = actual.query_length - actual.query_end if actual.reverse else actual.query_start
    right = actual.query_start if actual.reverse else actual.query_length - actual.query_end
    summary = tuple((n, op) for n, op in (
        (left, "S"), (min(query, reference), "M"),
        (max(0, query - reference), "I"), (max(0, reference - query), "D"), (right, "S"),
    ) if n)
    return reported.cigar == summary


def _linked_member(path, other, placement, observed_placements):
    matches = [a for a in path.others if _matches_reported(other.alignment, a)]
    if len(matches) != 1:
        return None
    reported = matches[0]
    if len(observed_placements) > 1:
        # A lossy SA summary cannot choose between observed alternative CIGARs
        # with identical spans, even when one happens to equal the summary.
        try:
            candidates = {
                _alignment(reported.contig, start, reverse, cigar)
                for ref, start, cigar, reverse in observed_placements
                if ref == placement[0]
            }
        except ValueError:
            return None
        if sum(_matches_reported(a, reported) for a in candidates) != 1:
            return None
    return reported


def _unambiguous_locus(path):
    """A breakpoint overlap cannot assign the same variant bases to both pieces."""
    start, end = path.query_interval or (path.alignment.query_start, path.alignment.query_end)
    return not any(
        (start < a.query_end and a.query_start < end) if start != end
        else (a.query_start < start < a.query_end)
        for a in path.others)


def compatible_phasing_alignments(constraints, read, other_read, observed_placements):
    """Accept mates normally; relax same-segment conflicts only for linked pieces."""
    if compatible_alignments(constraints, (other_read,)):
        return True
    other_constraints = dict(other_read.source_alignments)
    paths = dict(getattr(read, "source_alignment_paths", ()))
    other_paths = dict(getattr(other_read, "source_alignment_paths", ()))
    for segment in constraints.keys() & other_constraints.keys():
        if constraints[segment] == other_constraints[segment]:
            continue
        first, second = paths.get(segment), other_paths.get(segment)
        if (first is None or second is None or first.secondary != second.secondary
                or not (first.supplementary or second.supplementary)
                or not (_unambiguous_locus(first) and _unambiguous_locus(second))):
            return False
        to_second = _linked_member(first, second, other_constraints[segment], observed_placements[segment])
        to_first = _linked_member(second, first, constraints[segment], observed_placements[segment])
        if (to_first is None or to_second is None
                or set(first.others) - {to_second} != set(second.others) - {to_first}):
            return False
    return True
