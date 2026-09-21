"""Native read evidence, kept separate from per-base and mapping qualities."""


# Select individual tags: get_tags() would also decode potentially large signal,
# kinetics and modification arrays. Lowercase/X tags have producer-specific
# meanings; retain their values without assuming a platform from tag presence.
EVIDENCE_TAGS = (
    "RG", "PG", "NM", "AS", "NH", "HI", "nM", "ms", "s1", "s2", "dv", "de", "rl", "tp",
    "qs", "dx", "pi", "sp",
    "mg", "rm", "rq", "np", "ec", "ff", "ic", "is", "im", "rc",
    "CB", "CR", "CY", "UB", "UR", "UY", "RX", "QX", "MI", "XM",
)


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
