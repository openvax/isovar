"""Transcript-specific ATG origins, separate from observed initiation."""

from .genetic_code import standard_genetic_code


def summarize_orf_start_evidence(annotations):
    """Summarize start tiers without choosing the strongest overlapping model.

    Parameters
    ----------
    annotations : iterable of dict or None
        Results of :func:`annotate_orf_start`, for example from alternative
        occurrences of one exported ORF. ``None`` denotes unavailable metadata.

    Returns
    -------
    dict
        ``status``, ``tier`` and ``priority``. Lower priority values rank first
        as an annotation prior only. Empty or missing input is ``unavailable``;
        disagreement is ``ambiguous``. Both have null tier and priority. Equal
        tiers do not establish a unique transcript or biological RNA polarity.
    """
    annotations = list(annotations)
    if not annotations or any(a is None or not a.get("assessments") for a in annotations):
        return dict(status="unavailable", tier=None, priority=None)
    tiers = {(row["tier"], row["priority"]) for a in annotations for row in a["assessments"]}
    if len(tiers) != 1:
        return dict(status="ambiguous", tier=None, priority=None)
    tier, priority = tiers.pop()
    return dict(status="consistent", tier=tier, priority=priority)


def annotate_orf_start(sequence, positions, start, references):
    """Annotate an observed ATG against each supplied transcript it overlaps.

    Parameters
    ----------
    sequence : str
        Retained RNA-path sequence in the candidate's reading orientation.
    positions : sequence of tuple or None
        One ``(contig, zero-based genomic position, strand)`` per RNA base.
        Strand describes the path orientation, not proven biological polarity.
        ``None`` marks an unplaced base, which need not be a DNA insertion.
    start : int
        Zero-based query offset of an ATG. Invalid length, offset or codon
        raises ``ValueError``. The caller determines event linkage and stops.
    references : iterable of FusionReference
        Supplied spliced transcripts with annotation identity and optional CDS.
        An empty iterable makes placed starts outside the supplied models;
        it does not establish a genome-wide intergenic location.

    Returns
    -------
    dict
        Per-transcript ``assessments``, legacy-compatible reference ORF
        comparisons, and a conservative tier summary. Priorities 1, 2 and 4
        distinguish exact annotated CDS ATGs, 5-prime-UTR ATGs, and other
        sequence-only starts. Tier 3 is reserved for a future linked splice
        assessment and is never inferred here. Intronic coverage or a strong
        initiation context alone cannot promote a start. Initiation,
        translation and splice inclusion remain unassessed/unobserved.

        Reference readings distinguish an upstream ORF, overlapping ORF,
        unresolved stop, and an in-frame extension; a UTR ATG is not itself
        proof of a translated uORF. Alternative transcript assessments are
        retained, and conflicting tiers receive no single numeric priority.
    """
    if (len(sequence) != len(positions) or isinstance(start, bool) or not isinstance(start, int)
            or not 0 <= start <= len(sequence) - 3 or sequence[start:start + 3] != "ATG"):
        raise ValueError("ORF start requires an ATG offset and one placement per sequence base")
    codon_positions = positions[start:start + 3]
    placed = [p for p in codon_positions if p is not None]
    assessments, comparisons = [], []
    aa, stop = standard_genetic_code.translate(sequence[start:], first_codon_is_start=True)
    for reference in sorted(references, key=lambda r: (r.transcript_id, r.annotation, r.reference_name)):
        # Only overlapping models contribute; unrelated transcripts cannot
        # turn a supported UTR assignment into a spurious disagreement.
        if not any(p[0] == reference.contig and reference.exons[0][0] <= p[1] < reference.exons[-1][1]
                   for p in placed):
            continue
        row = dict(transcript_id=reference.transcript_id, annotation=reference.annotation,
                   reference_name=reference.reference_name, transcript_strand=reference.strand,
                   transcript_offset=None, annotated_cds_start=reference.cds_start,
                   tier="sequence_only", priority=4, region="junction_created",
                   reference_start_codon=None, reference_orf_type=None,
                   transcript_inclusion="unresolved", splice_inference="not_assessed")
        offsets = [reference.spliced_offset(p[1]) if p is not None and p[0] == reference.contig
                   and p[2] == reference.strand else None for p in codon_positions]
        offset = offsets[0]
        if len(placed) < 3:
            row["region"] = "partially_unplaced"
        elif all(p[0] == reference.contig and p[2] != reference.strand for p in placed):
            row["region"] = "antisense"
        elif offset is not None and offsets == [offset, offset + 1, offset + 2]:
            row.update(transcript_offset=offset, reference_start_codon=reference.sequence[offset:offset + 3],
                       transcript_inclusion="supplied_transcript_exon")
            if reference.cds_start is None:
                region = "noncoding_transcript"
            elif offset == reference.cds_start:
                region = "annotated_CDS_start"
            elif offset + 3 <= reference.cds_start:
                region = "five_prime_UTR"
            elif offset >= reference.cds_end:
                region = "three_prime_UTR"
            elif offset < reference.cds_start:
                region = "UTR_CDS_boundary"
            else:
                region = "internal_CDS"
            row["region"] = region
            if region == "five_prime_UTR":
                row.update(tier="five_prime_UTR_ATG", priority=2)
            if row["reference_start_codon"] == "ATG":
                ref_aa, ref_stop = standard_genetic_code.translate(
                    reference.sequence[offset:], first_codon_is_start=True)
                shared = 0
                for a, b in zip(aa, ref_aa):
                    if a != b:
                        break
                    shared += 1
                if region == "annotated_CDS_start":
                    row.update(tier="annotated_CDS_start", priority=1)
                if region == "five_prime_UTR":
                    if ref_stop and offset + 3 * (len(ref_aa) + 1) <= reference.cds_start:
                        row["reference_orf_type"] = "upstream_ORF"
                    elif (reference.cds_start - offset) % 3 == 0:
                        row["reference_orf_type"] = "in_frame_extension"
                    else:
                        row["reference_orf_type"] = "overlapping_ORF" if ref_stop else "stop_unresolved"
                comparisons.append(dict(
                    transcript_id=reference.transcript_id, annotation=reference.annotation,
                    transcript_offset=offset,
                    start_kind="annotated_start" if region == "annotated_CDS_start" else region,
                    annotated_cds_start=reference.cds_start, amino_acids=ref_aa, ends_with_stop_codon=ref_stop,
                    shared_prefix_amino_acids=shared, differs_from_reference_orf=(aa, stop) != (ref_aa, ref_stop)))
        else:
            for (_, left), (right, _) in zip(reference.exons, reference.exons[1:]):
                if all(p[0] == reference.contig and p[2] == reference.strand and left <= p[1] < right
                       for p in placed) and all(
                           b[1] - a[1] == (1 if reference.strand == "+" else -1)
                           for a, b in zip(placed, placed[1:])):
                    row.update(region="intronic", intron_interval=[left, right],
                               transcript_inclusion="not_established")
                    break
        assessments.append(row)
    if not assessments:
        if not placed:
            region = "unplaced"
        elif len(placed) < 3:
            region = "partially_unplaced"
        elif any(a[0] != b[0] or a[2] != b[2] or b[1] - a[1] != (1 if a[2] == "+" else -1)
                 for a, b in zip(placed, placed[1:])):
            region = "junction_created"
        else:
            region = "outside_supplied_transcripts"
        assessments.append(dict(transcript_id=None, annotation=None, reference_name=None,
                                region=region, tier="sequence_only", priority=4,
                                transcript_inclusion="unresolved", splice_inference="not_assessed"))
    result = dict(policy="isovar.orf_start.v1", query_offset=start,
                  genomic_positions=[list(p) if p is not None else None for p in codon_positions],
                  assessments=assessments, reference_comparisons=comparisons,
                  initiation_observed=False, translation_observed=False, rna_polarity="unresolved")
    result.update(summarize_orf_start_evidence([result]))
    return result
