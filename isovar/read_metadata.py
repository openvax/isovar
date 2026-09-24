"""Native read evidence, kept separate from per-base and mapping qualities."""

from collections import Counter


# Select individual tags: get_tags() would also decode potentially large signal,
# kinetics and modification arrays. Lowercase/X tags have producer-specific
# meanings; retain their values without assuming a platform from tag presence.
EVIDENCE_TAGS = (
    "RG", "PG", "NM", "AS", "NH", "HI", "nM", "ms", "s1", "s2", "dv", "de", "rl", "tp",
    "qs", "dx", "pi", "sp",
    "mg", "rm", "rq", "np", "ec", "ff", "ic", "is", "im", "rc",
    "CB", "CR", "CY", "UB", "UR", "UY", "RX", "QX", "MI", "XM",
)


def unique_header_entries(entries):
    """Index @RG or @PG header lines by ID, omitting missing and repeated IDs.

    A duplicated ID cannot identify one entry, and SAM readers accept lines
    without one, so neither is guessed.
    """
    counts = Counter(entry.get("ID") for entry in entries)
    return {entry["ID"]: entry for entry in entries if entry.get("ID") and counts[entry["ID"]] == 1}


class ProgramHistory:
    """The @PG chains of one SAM header. Broken or cyclic chains are unresolved."""

    def __init__(self, header):
        entries = header.get("PG", [])
        self.programs = unique_header_entries(entries)
        self.unique = len(self.programs) == len(entries)
        parents = {program.get("PP") for program in self.programs.values()}
        self.leaves = sorted(set(self.programs) - parents)

    def chain(self, key):
        """Programs from ``key`` back to its root, newest first, or None."""
        chain, seen = [], set()
        while key is not None:
            if key in seen or key not in self.programs:
                return None
            seen.add(key)
            chain.append(self.programs[key])
            key = chain[-1].get("PP")
        return chain

    def leaves_cover_all(self):
        """Whether walking back from the leaves reaches every program."""
        covered, pending = set(), list(self.leaves)
        while pending:
            key = pending.pop()
            if key in self.programs and key not in covered:
                covered.add(key)
                pending.append(self.programs[key].get("PP"))
        return covered == self.programs.keys()

    @staticmethod
    def pointer(read, read_group):
        """The record's PG tag, else its read group's PG entry; None if neither."""
        return read.get_tag("PG") if read.has_tag("PG") else read_group.get("PG")


def record_evidence(read):
    """Return JSON-ready alignment, consensus, barcode and quality evidence.

    Parameters
    ----------
    read : pysam.AlignedSegment
        Unmodified SAM/BAM record. Interpret tags alongside its RG/PG header.

    Returns
    -------
    dict
        Missing tags are None, including unavailable MAPQ 255 in the normalized
        field. The raw MAPQ and flags are also retained. No value replaces QUAL;
        OQ presence only signals a possible source of original qualities.
        Call on demand; ordinary read collection does not extract these tags.
    """
    return dict(
        sam_flag=read.flag,
        mapping_quality_raw=read.mapping_quality,
        mapping_quality=None if read.mapping_quality == 255 else read.mapping_quality,
        base_qualities_available=read.query_qualities is not None,
        original_base_qualities_tag_present=read.has_tag("OQ"),
        query_length=read.query_length,
        aligned_query_bases=sum(n for op, n in (read.cigartuples or ()) if op in (0, 1, 7, 8)),
        tags={tag: read.get_tag(tag) if read.has_tag(tag) else None for tag in EVIDENCE_TAGS})
