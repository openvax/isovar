"""Counts and reference comparisons for variant/RNA evidence reports.

The functions here operate on plain per-call records so source-specific audits
and the public VCF/BAM command use the same denominators. Exon overlap is an
annotation opportunity, not evidence for an expressed transcript.
"""

import csv
import gzip
from bisect import bisect_left
from collections import Counter, defaultdict
from pathlib import Path

from .protein_sequence_helpers import group_equivalent_translations, sort_protein_sequences
from .reference_context_helpers import reference_contexts_for_variant
from .transcript_compatibility import annotate_reads_with_transcript_compatibility
from .variant_sequence_helpers import initial_variant_sequences_from_reads, trim_variant_sequences


RNA_STATUSES = (
    "alternate_observed", "reference_only", "reference_and_other_only",
    "other_allele_only", "alignment_without_allele_call", "no_alignment_overlap",
    "unresolved_chromosome",
)
READ_BINS = ("0", "1", "2-4", "5-9", "10-49", "50-199", "200+")
AA_BINS = ("no local product", "0", "1", "2-4", "5-9", "10-19", "20+")


def direct_read_context_proteins(variant, alt_reads, creator):
    """Translate every high-depth ALT read without combinatorial overlap assembly."""
    contexts = reference_contexts_for_variant(
        variant, context_size=max(creator._cdna_sequence_length,
                                  creator.min_transcript_prefix_length))
    if not contexts:
        return []
    transcripts = {transcript.id: transcript
                   for context in contexts for transcript in context.transcripts}
    max_before, max_after = creator._variant_sequence_creator.flanking_sequence_lengths(
        len(alt_reads[0].allele))
    compatible = annotate_reads_with_transcript_compatibility(
        alt_reads, transcripts.values(), max_prefix_size=max_before,
        max_suffix_size=max_after)
    compatible = [read for read in compatible if read.compatible_transcript_ids]
    if not compatible:
        return []
    sequences = initial_variant_sequences_from_reads(
        compatible, max_nucleotides_before_variant=max_before,
        max_nucleotides_after_variant=max_after)
    sequences = trim_variant_sequences(sequences, creator.min_variant_sequence_coverage)
    if not sequences:
        return []
    translations = creator.all_pairs_translations(sequences, contexts)
    return sort_protein_sequences(
        group_equivalent_translations(translations),
        preference=creator.protein_sequence_preference,
        peptide_length=creator.protein_context_peptide_length,
        min_support_fraction=creator.min_protein_sequence_support_fraction)


def changed_reference_interval(pos, ref, alt):
    """Return a 1-based inclusive changed span and whether it is an insertion.

    Shared VCF padding and common suffix bases are removed. Pure insertions
    query the two flanking reference bases as candidate exon involvement.
    """
    if not ref or not alt or ref == alt:
        raise ValueError("Expected distinct literal REF and ALT alleles")
    left = 0
    while left < min(len(ref), len(alt)) and ref[left] == alt[left]:
        left += 1
    ref_end, alt_end = len(ref), len(alt)
    while ref_end > left and alt_end > left and ref[ref_end - 1] == alt[alt_end - 1]:
        ref_end -= 1
        alt_end -= 1
    first = pos + left
    if first == pos + ref_end:
        return max(1, first - 1), first, True
    return first, pos + ref_end - 1, False


def overlaps_exon(genome, contig, pos, ref, alt):
    """Whether any changed reference base (or insertion flank) meets an exon."""
    first, last, _ = changed_reference_interval(pos, ref, alt)
    return bool(genome.exon_ids_at_locus(contig.removeprefix("chr"), first, end=last))


def catalog_coding_index(records):
    """Index reference ncORF coding blocks for exact interval intersections."""
    blocks = defaultdict(list)
    for record in records:
        for first, last in record.coding_blocks:
            blocks[record.contig.removeprefix("chr")].append((first - 1, last, record.orf_id))
    result = {}
    for contig, intervals in blocks.items():
        intervals.sort()
        maximum = 0
        maxima = []
        for _, end, _ in intervals:
            maximum = max(maximum, end)
            maxima.append(maximum)
        result[contig] = ([item[0] for item in intervals], maxima, intervals)
    return result


def overlapping_catalog_orfs(index, contig, begin, end):
    """ORF IDs with a coding-block overlap in a zero-based half-open span."""
    contig = contig.removeprefix("chr")
    if contig not in index:
        return []
    starts, maxima, intervals = index[contig]
    i = bisect_left(starts, end) - 1
    ids = set()
    while i >= 0 and maxima[i] > begin:
        _, block_end, orf_id = intervals[i]
        if block_end > begin:
            ids.add(orf_id)
        i -= 1
    return sorted(ids)


def _yes(value):
    return value is True or value == "True"


def summarize(calls, annotations, proteins=(), sv=()):
    """Summarize source records with explicit, nested funnel denominators."""
    calls = list(calls)
    proteins = list(proteins)
    sv = list(sv)
    by_id = {row["call_id"]: row for row in calls}
    if len(by_id) != len(calls) or set(by_id) != set(annotations):
        raise ValueError("Calls and annotations must have the same unique call IDs")
    top = {}
    for row in proteins:
        key = row["call_id"]
        if key not in by_id:
            raise ValueError(f"Protein for unknown call {key}")
        if key not in top or int(row["rank"]) < int(top[key]["rank"]):
            top[key] = row
    counts = Counter()
    rna = Counter()
    classes = Counter()
    reads_by_class = {}
    exonic_reads_by_class = {}
    aa_counts = []
    for call in calls:
        ident = call["call_id"]
        annotation = annotations[ident]
        classes[call["class"]] += 1
        rna[call["rna_status"]] += 1
        if call["rna_status"] == "unresolved_chromosome":
            if ident in top or call.get("protein_count") not in (None, ""):
                raise ValueError(f"Unresolved chromosome cannot have assessed protein counts for {ident}")
            counts["unresolved_chromosomes"] += 1
            continue
        if "exon_overlap" not in annotation:
            raise ValueError("Exon overlap annotation is required for the funnel")
        exonic = _yes(annotation["exon_overlap"])
        alt = call["rna_status"] == "alternate_observed"
        product = int(call["protein_count"]) > 0
        if product and (not alt or ident not in top):
            raise ValueError(f"Protein count and RNA evidence disagree for {ident}")
        reads_by_class.setdefault(variant_type(call), []).append(int(call["alt_reads"]))
        if exonic:
            exonic_reads_by_class.setdefault(variant_type(call), []).append(int(call["alt_reads"]))
        counts["exonic"] += exonic
        counts["exonic_alt"] += exonic and alt
        counts["exonic_product"] += exonic and product
        counts["alt"] += alt
        counts["product"] += product
        counts["coding"] += _yes(annotation.get("reference_coding_effect"))
        counts["ncorf"] += _yes(annotation.get("orfdata_coding_block_overlap"))
        counts["frameshift"] += annotation.get("top_effect_type") == "FrameShift"
        counts["frameshift_alt"] += (
            annotation.get("top_effect_type") == "FrameShift" and alt)
        if alt:
            if product:
                selected = top[ident]
                count = int(selected["mutation_end_idx"]) - int(selected["mutation_start_idx"])
                if count < 0:
                    raise ValueError(f"Negative mutation interval for {ident}")
                aa_counts.append(count)
            else:
                aa_counts.append(None)
    counts["total"] = len(calls)
    counts["sv_total"] = len(sv)
    counts["sv_product"] = sum(row.get("protein_status", "").startswith("translated")
                               for row in sv)
    if not (counts["exonic_product"] <= counts["exonic_alt"] <= counts["exonic"] <= len(calls)):
        raise ValueError("Exonic funnel stages are not nested")
    return dict(counts=counts, rna=rna, classes=classes, top=top,
                reads_by_class=reads_by_class,
                exonic_reads_by_class=exonic_reads_by_class, aa_counts=aa_counts)


def read_bin(count):
    if count == 0:
        return "0"
    if count == 1:
        return "1"
    if count < 5:
        return "2-4"
    if count < 10:
        return "5-9"
    if count < 50:
        return "10-49"
    if count < 200:
        return "50-199"
    return "200+"


def variant_type(call):
    """Group literal small variants for read-support plots."""
    ref, alt = call["ref"], call["alt"]
    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(ref) > len(alt):
        return "deletion"
    if len(alt) > len(ref):
        return "insertion"
    return "replacement"


def aa_bin(count):
    if count is None:
        return "no local product"
    if count == 0:
        return "0"
    if count == 1:
        return "1"
    if count < 5:
        return "2-4"
    if count < 10:
        return "5-9"
    if count < 20:
        return "10-19"
    return "20+"


def mutant_ninemers(protein):
    """Exact 9-mers touching the mutant interval or a deletion junction."""
    sequence = protein["protein_sequence"]
    start, end = int(protein["mutation_start_idx"]), int(protein["mutation_end_idx"])
    if not 0 <= start <= end <= len(sequence):
        raise ValueError("Mutation interval lies outside the local protein sequence")
    if len(sequence) < 9:
        return set()
    return {sequence[i:i + 9] for i in range(len(sequence) - 8)
            if ((i < end and i + 9 > start) if start != end
                else (i < start < i + 9))}


def reference_ninemer_matches(proteins, fasta):
    """Find candidate mutant 9-mers present in a supplied reference FASTA.

    Only candidate 9-mers are retained, avoiding a whole-proteome k-mer index.
    The result is an exact reference screen, not a patient-specific self model.
    """
    targets = set().union(*(mutant_ninemers(p) for p in proteins))
    found = set()
    opener = gzip.open if str(fasta).endswith(".gz") else open
    with opener(fasta, "rt") as handle:
        sequence = []
        for line in handle:
            if line.startswith(">"):
                joined = "".join(sequence)
                found.update(joined[i:i + 9] for i in range(len(joined) - 8)
                             if joined[i:i + 9] in targets)
                sequence = []
            else:
                sequence.append(line.strip())
        joined = "".join(sequence)
        found.update(joined[i:i + 9] for i in range(len(joined) - 8)
                     if joined[i:i + 9] in targets)
    return found


def self_novelty(protein, reference_matches):
    """Return (absent, total) mutation-overlapping exact reference 9-mers."""
    kmers = mutant_ninemers(protein)
    return len(kmers - reference_matches), len(kmers)


def read_tsv(path):
    """Read a plain or gzip TSV into dictionaries."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(Path(path), "rt", newline="") as source:
        return list(csv.DictReader(source, delimiter="\t"))
