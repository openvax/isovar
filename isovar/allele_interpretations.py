"""Which competing interpretations of one locus do the RNA reads support?

Variant callers and catalogues can describe the same locus differently: one
SNV, a multi-base substitution, two adjacent SNVs, or a nearby deletion.
`reconcile_allele_interpretations` compares each read's exact sequence across
the locus with every interpretation's haplotype. Interpretations that give the
same sequence cannot be told apart by any RNA and are reported together. The
result, ``isovar.allele_interpretations.v2``, is a candidate by evidence table;
it does not choose an interpretation or translate one.
"""

from collections import defaultdict

from varcode import Variant

from . import __version__
from .allele_read import AlleleRead
from .default_parameters import INTERPRETATION_FLANK, MIN_INTERPRETATION_FRAGMENTS
from .dna import complement_dna
from .read_collector import ReadCollector
from .read_identity import fragment_ids, segment_identity, source_read_ids
from .rna_evidence import EvidenceSets
from .variant_helpers import base0_interval_for_variant, require_literal_variant

SCHEMA = "isovar.allele_interpretations.v2"


def _exon_bases(variants, start, end):
    """Forward-strand reference bases from exons containing [start, end), keyed by 0-based position."""
    genome = variants[0].genome
    contig = variants[0].contig
    bases = {}
    try:
        transcripts = genome.transcripts_at_locus(contig, start + 1, max(start + 1, end))
    except (AttributeError, ValueError):
        return bases, {}
    exons = {}
    for transcript in transcripts:
        sequence = getattr(transcript, "sequence", None)
        if not sequence:
            continue
        for exon in transcript.exons:
            # Exon coordinates are 1-based inclusive; bases are keyed 0-based.
            if not (exon.start - 1 <= start and end <= exon.end):
                continue
            exons[(exon.start - 1, exon.end)] = True
            for position in range(exon.start - 1, exon.end):
                offset = transcript.spliced_offset(position + 1)
                base = sequence[offset]
                bases.setdefault(position, base if transcript.strand == "+" else complement_dna(base))
    return bases, exons


def _window(variants, flank):
    """The 0-based window, its reference sequence (None if unknown), a reason, and the padding."""
    intervals = [base0_interval_for_variant(v) for v in variants]
    start, end = min(a for a, _ in intervals), max(b for _, b in intervals)
    locus_start, locus_end = start, end
    known = {}
    for variant, (a, b) in zip(variants, intervals):
        for offset, base in enumerate(variant.ref):
            if known.setdefault(a + offset, base) != base:
                raise ValueError("Interpretations disagree on the reference base at %s:%d" % (
                    variant.contig, a + offset + 1))
    exon_bases, exons = _exon_bases(variants, start, end)
    for position, base in known.items():
        if position in exon_bases and exon_bases[position] != base:
            raise ValueError("Reference allele at %s:%d disagrees with the annotation" % (
                variants[0].contig, position + 1))
    if exons:
        # Pad within the exon containing the locus, never across a splice junction.
        exon_start = max(a for a, _ in exons)
        exon_end = min(b for _, b in exons)
        start, end = max(exon_start, start - flank), min(exon_end, end + flank)
    bases = {**exon_bases, **known}
    padding = [locus_start - start, end - locus_end]
    missing = [p for p in range(start, end) if p not in bases]
    if missing:
        return (start, end), None, "reference bases unknown at %d position(s); no annotated exon covers them" % len(
            missing), padding
    return (start, end), "".join(bases[p] for p in range(start, end)), None, padding


def _haplotype(reference, window_start, variants):
    """The window's sequence with ``variants`` applied."""
    edits = sorted((base0_interval_for_variant(v), v.alt) for v in variants)
    for ((_, end), _), ((start, _), _) in zip(edits, edits[1:]):
        if start < end:
            raise ValueError("An interpretation's variants overlap each other")
    sequence, position = [], window_start
    for (start, end), alt in edits:
        sequence.append(reference[position - window_start:start - window_start])
        sequence.append(alt)
        position = end
    sequence.append(reference[position - window_start:])
    return "".join(sequence)


def _difference(reference, window_start, allele):
    """The allele as one trimmed replacement of the reference window (1-based start), or None."""
    if allele == reference:
        return None
    prefix = 0
    while prefix < min(len(reference), len(allele)) and reference[prefix] == allele[prefix]:
        prefix += 1
    suffix = 0
    while (suffix < min(len(reference), len(allele)) - prefix
           and reference[-1 - suffix] == allele[-1 - suffix]):
        suffix += 1
    return dict(start=window_start + prefix + 1, ref=reference[prefix:len(reference) - suffix],
                alt=allele[prefix:len(allele) - suffix])


def _variant_row(variant):
    return dict(contig=variant.contig, start=variant.start, ref=variant.ref, alt=variant.alt,
                description=variant.short_description)


def reconcile_allele_interpretations(
        interpretations, alignment_file, *, sample_id, source, read_collector=None,
        flank=INTERPRETATION_FLANK, min_fragments=MIN_INTERPRETATION_FRAGMENTS, cell_umi_labels=False):
    """
    Compare competing interpretations of one locus with the RNA reads.

    Parameters
    ----------
    interpretations : dict
        Interpretation ID to a list of varcode.Variant, applied together as
        one haplotype. All variants must be literal, on one contig and one
        reference assembly.
    alignment_file : pysam.AlignmentFile
        RNA alignments.
    sample_id, source : str
        The evidence scope of the hashed read IDs, as in other exports.
    read_collector : ReadCollector, optional
        Read selection settings; the defaults otherwise.
    flank : int
        Exonic reference bases added on each side of the locus, within its
        exon, so that an indel an aligner placed elsewhere in a repeat still
        falls inside the compared window.
    min_fragments : int
        Fragments needed to call an interpretation supported, or to call it
        contradicted when no fragment carries it.
    cell_umi_labels : bool
        Also count the UMIs and cells behind each support from the reads' cell
        barcodes and UMIs; otherwise those fields are null (not assessed).

    Returns
    -------
    dict
        ``isovar.allele_interpretations.v2``: the compared window, every
        interpretation's haplotype and ``rna_status``, the reference allele's
        support, every observed allele with the interpretations it matches,
        and hashed ``evidence_sets``. Every support is the RNA support record
        (see `isovar.rna_evidence.rna_support`). Only reads spanning the
        whole window are informative; the rest are counted but set aside.

    Raises
    ------
    ValueError
        For no interpretations, symbolic or mixed-contig variants, overlapping
        variants within one interpretation, or disagreeing reference alleles.
    """
    if not interpretations:
        raise ValueError("No interpretations to compare")
    interpretations = {str(key): list(variants) for key, variants in interpretations.items()}
    variants = [v for vs in interpretations.values() for v in vs]
    if not variants or any(not vs for vs in interpretations.values()):
        raise ValueError("Every interpretation needs at least one variant")
    for variant in variants:
        require_literal_variant(variant)
    if len({(v.contig, v.reference_name) for v in variants}) > 1:
        raise ValueError("All interpretations must be on one contig of one reference")
    if read_collector is None:
        read_collector = ReadCollector()
    scope = [sample_id, source]
    (start, end), reference, reason, padding = _window(variants, flank)
    result = dict(
        schema=SCHEMA, isovar_version=__version__, sample_id=sample_id, source=source,
        evidence_scope=scope, interval_convention="zero_based_half_open",
        locus=dict(reference_name=variants[0].reference_name, contig=variants[0].contig,
                   window=[start, end], padding=padding, reference_sequence=reference),
        parameters=dict(flank=flank, min_fragments=min_fragments),
        interpretation="Exact read sequences across the window; no interpretation is chosen or translated. "
                       "Counts are reads, fragments, UMIs and cells, not molecules.")
    haplotypes = {} if reference is None else {
        key: _haplotype(reference, start, vs) for key, vs in interpretations.items()}
    if reference is None:
        result.update(status="unassessed", reason=reason, candidates=[
            dict(candidate_id=key, variants=[_variant_row(v) for v in vs], haplotype=None,
                 rna_status="unassessed", indistinguishable_from=[])
            for key, vs in interpretations.items()], observed_alleles=[], evidence_sets={})
        return result

    chromosome = read_collector._infer_chromosome_name(variants[0].contig, set(alignment_file.references))
    locus_reads = [] if chromosome is None else read_collector.get_locus_reads(
        alignment_file, chromosome, start, end)
    reads = [read for read in (AlleleRead.from_locus_read(r) for r in locus_reads) if read is not None]
    overlapping = set() if chromosome is None else {
        segment_identity(r) for r in alignment_file.fetch(chromosome, start, max(end, start + 1))
        if read_collector.alignment_filter_reason(r) is None}
    alleles_by_segment = defaultdict(set)
    for read in reads:
        for key in source_read_ids(read):
            alleles_by_segment[key].add(read.allele)
    conflicting = {key for key, alleles in alleles_by_segment.items() if len(alleles) > 1}
    informative = [read for read in reads if not source_read_ids(read) & conflicting]
    labels = None
    if cell_umi_labels:
        from .cell_evidence import CellUmiAlleles
        labels = CellUmiAlleles(alignment_file, sample_id=sample_id, source=source, read_collector=read_collector
                                ).labels_at(variants[0].contig, start, end,
                                            [k for read in informative for k in source_read_ids(read)])
    evidence = EvidenceSets(scope)

    def supports(reads):
        return evidence.support(reads, labels)

    by_allele = defaultdict(list)
    for read in informative:
        by_allele[read.allele].append(read)
    informative_fragments = len(fragment_ids(informative))

    def matching(allele):
        return (["reference"] if allele == reference else []) + sorted(
            key for key, haplotype in haplotypes.items() if haplotype == allele)

    candidates = []
    for key, vs in interpretations.items():
        haplotype = haplotypes[key]
        support = supports(by_allele.get(haplotype, []))
        status = ("supported" if support["fragments"] >= min_fragments else
                  "contradicted_by_informative_evidence"
                  if not support["fragments"] and informative_fragments >= min_fragments else "insufficient_RNA")
        candidates.append(dict(
            candidate_id=key, variants=[_variant_row(v) for v in vs], haplotype=haplotype,
            same_as_reference=haplotype == reference, rna_status=status, rna_support=support,
            indistinguishable_from=sorted(k for k, h in haplotypes.items() if h == haplotype and k != key)))
    result.update(
        status="assessed" if chromosome is not None else "unassessed",
        reason=None if chromosome is not None else "contig %s is not in the alignment file" % variants[0].contig,
        candidates=candidates,
        reference_allele=dict(haplotype=reference, rna_support=supports(by_allele.get(reference, []))),
        observed_alleles=[dict(haplotype=allele, matches=matching(allele),
                               difference_from_reference=_difference(reference, start, allele),
                               rna_support=supports(group))
                          for allele, group in sorted(by_allele.items(), key=lambda item: (-len(item[1]), item[0]))],
        informative=supports(informative),
        set_aside=dict(
            reads_not_spanning_window=len(overlapping - {k for r in reads for k in source_read_ids(r)}),
            conflicting_reads=len(conflicting)),
        evidence_sets={key: evidence.sets[key] for key in sorted(evidence.sets)})
    return result


INPUT_KEYS = frozenset(("reference_name", "interpretations"))
VARIANT_KEYS = frozenset(("contig", "start", "ref", "alt"))


def interpretations_from_dict(data, genome=None):
    """
    Interpretations from JSON: ``{"reference_name": ..., "interpretations":
    {ID: [{"contig", "start", "ref", "alt"}, ...]}}``, with ``start`` 1-based
    as in a VCF.

    Parameters
    ----------
    data : dict
    genome : str or pyensembl.Genome, optional
        Annotation for the variants; ``reference_name`` otherwise.

    Returns
    -------
    dict
        Interpretation ID to a list of varcode.Variant.

    Raises
    ------
    ValueError
        For missing or unknown keys.
    """
    if set(data) != INPUT_KEYS:
        raise ValueError("Interpretation input needs exactly the keys %s" % ", ".join(sorted(INPUT_KEYS)))
    genome = genome or data["reference_name"]
    result = {}
    for key, variants in data["interpretations"].items():
        for row in variants:
            if set(row) != VARIANT_KEYS:
                raise ValueError("Interpretation %s: each variant needs exactly %s" % (
                    key, ", ".join(sorted(VARIANT_KEYS))))
        result[key] = [Variant(str(row["contig"]), int(row["start"]), row["ref"], row["alt"], genome)
                       for row in variants]
    return result
