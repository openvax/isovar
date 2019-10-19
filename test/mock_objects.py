from __future__ import print_function, division, absolute_import

import pysam

from isovar.translation import Translation
from isovar.protein_sequence import ProteinSequence
from isovar.variant_orf import VariantORF

class MockAlignmentFile(object):
    """
    Used instead of real AlignmentFile objects for test.
    """
    def __init__(self, references, reads):
        self.references = tuple(references)
        self.reads = reads

    def fetch(self, *args, **kwargs):
        return self.reads

    @property
    def filename(self):
        return "MOCK-READS.bam"

def make_pysam_read(
        seq,
        cigar,
        mdtag=None,
        name="dummy",
        mapq=10,
        baseq=30,
        reference_start=0,
        reference_id=0):
    read = pysam.AlignedSegment()
    read.seq = seq
    read.cigarstring = cigar
    if mdtag:
        read.set_tag("MD", mdtag)
    read.qname = name
    read.mapq = mapq
    read.reference_start = reference_start
    read.reference_id = reference_id
    qualities_string = pysam.qualities_to_qualitystring([baseq] * len(seq))
    read.qual = qualities_string.encode("ascii")
    return read


class MockAlleleRead(object):
    def __init__(self, name="mock-read"):
        self.name = name

class MockVariantSequence(object):
    def __init__(self, n_reads):
        self.n_reads = n_reads

    @property
    def reads(self):
        return [MockAlleleRead("mock-read-%d" % i) for i in range(self.n_reads)]

class MockReferenceContext(object):
    pass

def make_dummy_translation(
        amino_acids="MKHW",  # ATG=M|AAA=K|CAC=H|TGG=W
        cdna_sequence="CCCATGAAACACTGGTAG",
        offset_to_first_complete_codon=3,
        variant_cdna_interval_start=8,  # assuming variant was AAC>AAA
        variant_cdna_interval_end=9,
        variant_aa_interval_start=1,
        variant_aa_interval_end=2,
        num_mismatches=1,
        n_variant_reads=1):
    """
    Create mock Translation object with minimal information needed to
    get used successfully by ProteinSequence.
    """
    varseq_in_orf = VariantORF(
        cdna_sequence=cdna_sequence,
        offset_to_first_complete_codon=offset_to_first_complete_codon,
        variant_cdna_interval_start=variant_cdna_interval_start,
        variant_cdna_interval_end=variant_cdna_interval_end,
        reference_cdna_sequence_before_variant=cdna_sequence[:variant_cdna_interval_start],
        reference_cdna_sequence_after_variant=cdna_sequence[variant_cdna_interval_end:],
        num_mismatches_before_variant=num_mismatches,
        num_mismatches_after_variant=0)

    return Translation(
        variant_orf=varseq_in_orf,
        amino_acids=amino_acids,
        contains_mutation=True,
        mutation_start_idx=variant_aa_interval_start,
        mutation_end_idx=variant_aa_interval_end,
        frameshift=False,
        ends_with_stop_codon=False,
        untrimmed_variant_sequence=MockVariantSequence(n_reads=n_variant_reads),
        reference_context=MockReferenceContext())


def make_dummy_protein_sequence(
        n_supporting_variant_reads,
        n_supporting_variant_sequences,
        n_supporting_reference_transcripts,
        n_total_variant_sequences=None,
        n_total_variant_reads=None,
        n_total_reference_transcripts=None,
        amino_acids="MKHW",  # ATG=M|AAA=K|CAC=H|TGG=W
        cdna_sequence="CCCATGAAACACTGGTAG",
        variant_cdna_interval_start=8,  # assuming variant was AAC>AAA
        variant_cdna_interval_end=9,
        variant_aa_interval_start=1,
        variant_aa_interval_end=2,
        num_mismatches=1):
    """
    Creates ProteinSequence object with None filled in for most fields
    """
    if n_total_variant_reads is None:
        n_total_variant_reads = n_supporting_variant_reads

    if n_total_variant_sequences is None:
        n_total_variant_sequences = n_supporting_variant_sequences

    if n_total_reference_transcripts is None:
        n_total_reference_transcripts = n_total_reference_transcripts

    assert n_supporting_variant_sequences <= n_supporting_variant_reads
    assert n_supporting_variant_sequences <= n_total_variant_sequences
    assert n_supporting_reference_transcripts <= n_total_reference_transcripts

    n_translations = n_total_reference_transcripts * n_total_variant_sequences

    translation = make_dummy_translation(
        amino_acids=amino_acids,
        cdna_sequence=cdna_sequence,
        offset_to_first_complete_codon=3,
        variant_cdna_interval_start=variant_cdna_interval_start,  # assuming variant was AAC>AAA
        variant_cdna_interval_end=variant_cdna_interval_end,
        variant_aa_interval_start=variant_aa_interval_start,
        variant_aa_interval_end=variant_aa_interval_end,
        num_mismatches=num_mismatches,
        n_variant_reads=n_total_variant_reads)

    return ProteinSequence.from_translations([translation] * n_translations)