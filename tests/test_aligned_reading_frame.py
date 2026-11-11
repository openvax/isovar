"""Regression for #265: observed upstream indels must not be rephased away."""

from copy import copy
from types import SimpleNamespace
from unittest.mock import patch

import pytest
from varcode import Variant

from isovar import ReadCollector, ProteinSequenceCreator
from isovar.aligned_reference_prefix import SplicedCoordinates, aligned_reference_prefix
from isovar.allele_read import AlleleRead
from isovar.dna import reverse_complement_dna
from isovar.reference_context import ReferenceContext
from isovar.variant_orf import VariantORF
from isovar.variant_sequence import VariantSequence
from .mock_objects import MockAlignmentFile, make_pysam_read
from .data.osteosarc.expansion.references import translate


def inputs(prefix_length=30, strand="+", start_codon=False):
    prefix = ("CCC" + "ATG" + "A" * (prefix_length - 6) if start_codon
              else "A" * prefix_length)
    suffix = "CCCTTTGGG" * 3
    sequence = prefix + "G" + suffix
    variant = Variant("1", 130 if strand == "+" else 127,
                      "T" if strand == "+" else "A", "G" if strand == "+" else "C", "GRCh38")
    context = ReferenceContext(
        strand=strand, sequence_before_variant_locus=("CCCATG" + "A" * 24 if start_codon else "A" * 30),
        sequence_at_variant_locus="T", sequence_after_variant_locus=suffix,
        offset_to_first_complete_codon=3 if start_codon else 0,
        contains_start_codon=start_codon, overlaps_start_codon=False,
        contains_five_prime_utr=start_codon, amino_acids_before_variant="K" * 10,
        variant=variant, transcripts=[SimpleNamespace(id="TX", name="TX", strand=strand,
                                                      exon_intervals=((100, 157),))])
    return variant, context, sequence


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("prefix_length,cigars,protein", [
    (30, ("58M", "58M"), "KKKKKKKKKKALWALWALW"),
    (31, ("15M1I43M", "43M1I15M"), "KKKKKKKKKKSPLGPLGPL"),
    (29, ("15M1D42M", "42M1D15M"), "KKKKKKKKKKPFGPFGPFG"),
])
def test_upstream_cigar_indel_transfers_frame(strand, prefix_length, cigars, protein):
    variant, context, sequence = inputs(prefix_length, strand)
    genomic = sequence if strand == "+" else reverse_complement_dna(sequence)
    records = [make_pysam_read(genomic, cigars[strand == "-"], name="read%d" % i,
                               reference_start=99) for i in range(3)]
    evidence = ReadCollector().read_evidence_for_variant(variant, MockAlignmentFile(["1"], records))
    with patch("isovar.protein_sequence_creator.reference_contexts_for_variant", return_value=[context]):
        translations = ProteinSequenceCreator(protein_sequence_preference="support").translate_variant_reads(
            variant, evidence.alt_reads)
    assert translations
    assert translate(sequence)[0] == protein
    assert {t.amino_acids for t in translations} == {protein}
    assert {t.variant_orf.num_mismatches_before_variant for t in translations} == {abs(prefix_length - 30)}


def test_conflicting_alignment_is_not_a_sequence_only_fallback():
    _, context, sequence = inputs()
    reads = [AlleleRead(sequence[:30], "G", sequence[31:], name="r%d" % i,
                        reference_blocks=((0, 58, 99 + i, 157 + i),)) for i in range(2)]
    variant_sequence = VariantSequence(sequence[:30], "G", sequence[31:], reads)
    assert VariantORF.from_variant_sequence_and_reference_context(variant_sequence, context) is None


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("disrupt_start", [False, True])
def test_annotated_start_is_transferred_or_rejected(strand, disrupt_start):
    variant, context, sequence = inputs(strand=strand, start_codon=True)
    position = 4 if disrupt_start else 2
    sequence = sequence[:position] + "C" + sequence[position:]
    cigar = "%dM1I%dM" % ((58 - position, position) if strand == "-" else (position, 58 - position))
    genomic = reverse_complement_dna(sequence) if strand == "-" else sequence
    record = make_pysam_read(genomic, cigar, name="start", reference_start=99)
    evidence = ReadCollector().read_evidence_for_variant(variant, MockAlignmentFile(["1"], [record]))
    read, = evidence.alt_reads
    candidate = VariantSequence(read.prefix, read.allele, read.suffix, [read])
    orf = VariantORF.from_variant_sequence_and_reference_context(candidate, context)
    if disrupt_start:
        assert orf is None
    else:
        assert orf.offset_to_first_complete_codon == 4
        assert orf.cdna_sequence[4:7] == "ATG"


def test_spliced_assembly_transfers_frame_across_two_partial_reads():
    variant, context, sequence = inputs(31)
    # Split the first exon after 15 nt, insert a base at the next exon start.
    variant = Variant("1", 230, "T", "G", "GRCh38")
    context.variant = variant
    context.transcripts = (SimpleNamespace(exon_intervals=((100, 114), (215, 257))),)
    records = [make_pysam_read(sequence[:44], "15M100N1I28M", name="left", reference_start=99),
               make_pysam_read(sequence[21:], "38M", name="right", reference_start=219)]
    evidence = ReadCollector().read_evidence_for_variant(variant, MockAlignmentFile(["1"], records))
    assert len(evidence.alt_reads) == 2
    orf = VariantORF.from_variant_sequence_and_reference_context(
        VariantSequence(sequence[:31], "G", sequence[32:], evidence.alt_reads), context)
    assert orf.offset_to_first_complete_codon == 0
    assert orf.num_mismatches_before_variant == 1
    assert translate(orf.in_frame_cdna_sequence)[0] == "KKKKKKKKKKSPLGPLGPL"


@pytest.mark.parametrize("read_id,protein", [
    ("1634caa5-09e3-4e97-bfa3-e117841b46cb_0", "VKENKDDLNHVDLNVCTSFRARV"),
    ("b0884250-7f95-40e9-bfd6-fde0ff6a6806_0", "VKENKDDLNHVDLNVCTSFWALV"),
])
def test_sid_cd109_observed_deletion_does_not_reset_frame(read_id, protein):
    # Original ONT records, digest-pinned with their acquisition provenance.
    # This tests translation of *observed* local RNA, not whether the deletion
    # is biological rather than a sequencing/alignment error.
    from examples.osteosarc_context_figures import load_context, original_reads

    data = load_context()
    model = data["models"]["ENST00000287097"]
    record = next(r for r in original_reads(data["cd109"][0]) if r.query_name == read_id)
    coordinates = {g: q for q, g in record.get_aligned_pairs(matches_only=True)}
    a, v, b = [coordinates[g] for g in (73815070, 73818404, 73818475)]
    assert v - a == 69  # 70 reference bases: one directly observed deletion.
    sequence = record.query_sequence[a:b + 1]
    prefix, alt, suffix = sequence[:v - a], sequence[v - a], sequence[v - a + 1:]
    blocks = AlleleRead._reference_blocks(record.get_reference_positions(full_length=True)[a:b + 1])
    read = AlleleRead(prefix, alt, suffix, read_id, reference_blocks=blocks)
    variant = Variant("6", 73818405, "G", "T", "GRCh38")
    positions = [p for start, end in model["exons"] for p in range(start - 1, end)]
    t0, tv, t1 = [positions.index(g) for g in (73815070, 73818404, 73818475)]
    assert (t0 - model["cds_start"]) % 3 == 0
    context = ReferenceContext(
        strand="+", sequence_before_variant_locus=model["cdna"][t0:tv],
        sequence_at_variant_locus="G", sequence_after_variant_locus=model["cdna"][tv + 1:t1 + 1],
        offset_to_first_complete_codon=0, contains_start_codon=False, overlaps_start_codon=False,
        contains_five_prime_utr=False, amino_acids_before_variant="",
        variant=variant, transcripts=[SimpleNamespace(exon_intervals=model["exons"])])
    orf = VariantORF.from_variant_sequence_and_reference_context(
        VariantSequence(prefix, alt, suffix, [read]), context)
    assert orf is not None
    assert orf.offset_to_first_complete_codon == 0
    assert orf.cdna_sequence == sequence
    assert translate(orf.in_frame_cdna_sequence)[0] == protein


def test_sequence_geometry_keeps_insertions_gaps_and_disagreements():
    # Query 2 is inserted; the genomic jump is a splice or deletion.
    long = AlleleRead("ACGT", "G", "TT", name="long", reference_blocks=((0, 2, 10, 12), (3, 7, 112, 116)))
    # Its first base is unaligned and its last lies beyond the sequence.
    short = AlleleRead("T", "G", "TTA", name="short", reference_blocks=((1, 5, 113, 117),))
    sequence = VariantSequence("ACGT", "G", "TT", [long, short])
    assert sequence.aligned_positions() == ((10,), (11,), (None,), (112,), (113,), (114,), (115,))
    assert sequence.aligned_positions() is sequence.aligned_positions()
    shifted = AlleleRead("ACGT", "G", "TT", name="shifted", reference_blocks=((0, 7, 109, 116),))
    disputed = VariantSequence("ACGT", "G", "TT", [long, shifted]).aligned_positions()
    assert (disputed[0], disputed[2], disputed[3]) == ((10, 109), (111, None), (112,))
    bare = VariantSequence("ACGT", "G", "TT", [AlleleRead("ACGT", "G", "TT", name="bare")])
    assert bare.aligned_positions() == ((),) * 7
    # A sequence pickled before the cache slot existed still computes.
    del bare._aligned_positions_cache
    assert bare.aligned_positions() == ((),) * 7


@pytest.mark.parametrize("reverse", [False, True])
def test_spliced_coordinates_follow_transcript_orientation(reverse):
    # One-base, adjacent and unsorted exons; 1-based inclusive as in PyEnsembl.
    exons = ((40, 44), (11, 11), (20, 24), (25, 30))
    spliced = SplicedCoordinates(SimpleNamespace(exon_intervals=exons), reverse)
    bases = sorted(g for a, b in exons for g in range(a - 1, b))
    for g in range(5, 50):
        upstream = [b for b in bases if (b >= g if reverse else b < g)]
        assert spliced.boundary(g) == len(upstream)
        assert spliced.base(g) == ((bases[::-1] if reverse else bases).index(g) if g in bases else None)


class WalkedBlocks(tuple):
    """Alignment blocks that count each traversal of their geometry."""

    walks = 0

    def __iter__(self):
        WalkedBlocks.walks += 1
        return super().__iter__()


@pytest.mark.parametrize("strand", ["+", "-"])
def test_coding_models_share_one_walk_of_each_read(strand):
    # #690: alternative ORFs of one RNA sequence must not rescan its reads.
    _, context, sequence = inputs(strand=strand)
    genomic = sequence if strand == "+" else reverse_complement_dna(sequence)
    cut = 30 if strand == "+" else 27
    reads = [AlleleRead(genomic[:cut], genomic[cut], genomic[cut + 1:], name="r%d" % i) for i in range(4)]
    for read in reads:
        read.reference_blocks = WalkedBlocks(((0, 58, 99, 157),))
    variant_sequence = VariantSequence(genomic[:cut], genomic[cut], genomic[cut + 1:], reads)
    WalkedBlocks.walks = 0
    for size, codon in [(30, 0), (24, 1), (18, 2), (12, 0), (6, 1)]:
        model = copy(context)
        model.sequence_before_variant_locus = "A" * size
        model.offset_to_first_complete_codon = codon
        assert aligned_reference_prefix(variant_sequence, model) == ("A" * size, "A" * size, codon, 0)
    assert WalkedBlocks.walks == len(reads)
