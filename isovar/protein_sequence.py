# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function, division, absolute_import

from collections import namedtuple, OrderedDict

import pandas as pd
from skbio import DNA

from .logging import create_logger
from .reference_context import (
    ReferenceContext,
    reference_contexts_for_variant
)
from .variant_sequence import variant_sequences_generator, VariantSequence

logger = create_logger(__name__)

# When multiple distinct RNA sequence contexts and/or reference transcripts
# give us the same translation then we group them into a single object
# which summarizes the supporting read count and reference transcript ORFs
# for a unique protein sequence.
#
# This object combines the information of a ReferenceContext
# (which establishes the reading frame for particular reference sequence) with
# a VariantSequence (which contains a cDNA sequence and its supporting reads).

ProteinSequence = namedtuple(
    "ProteinSequence",
    [
        # translated sequence of a variant sequence in the ORF established
        # by a reference context
        "amino_acids",
        # half-open interval coordinates for variant amino acids
        # in the translated sequence
        "variant_aa_interval_start",
        "variant_aa_end_offset",
        "reference_context"
        "variant_sequence",
    ])

MIN_READS_SUPPORTING_RNA_SEQUENCE = 3
MIN_TRANSCRIPT_PREFIX_LENGTH = 15
MAX_TRANSCRIPT_MISMATCHES = 2
PROTEIN_FRAGMENT_LEGNTH = 25
MAX_SEQUENCES_PER_VARIANT = 5


def translate(
        reference_context,
        cdna_sequence_before_variant,
        cdna_sequence_of_variant,
        cdna_sequence_after_variant,
        transcript_id,
        transcript_reading_frame_at_start_of_context_sequence,
        transcript_full_cdna,
        transcript_full_protein,
        variant_base1_offset_in_transcript):
    """
    Try to  translate the detected cDNA sequence context around a variant using
    the reading frame of reference transcripts (which must share sequences
    around the variant).

    Parameters
    ----------
    reference_context : ReferenceContext
        Has the following fields:
         - strand
         - transcript_ids
         - transcript_names
         - gene
         - reference_cdna_sequence_around_variant
         - reading_frame_at_start_of_context
         - number_cdna_bases_before_variant
         - number_cdna_bases_after_variant
         - reference_protein_sequence_around_variant

    To determine the reading frame we assume that the reference transcript
    sequence up to the variant largely matches the sequence we detected
    from RNA. If this is violated then the returned
    TranslationFromReferenceContext will be incomplete.
    """
    if transcript_reading_frame_at_start_of_context_sequence == 1:
        # if we're 1 nucleotide into the codon then we need to shift
        # over two more to restore the ORF
        orf_offset = 2
    elif transcript_reading_frame_at_start_of_context_sequence == 2:
        orf_offset = 1
    else:
        orf_offset = 0

    logger.info("ORF offset into sequence %s_%s_%s from transcript %s: %d" % (
        cdna_sequence_before_variant,
        cdna_sequence_of_variant,
        cdna_sequence_after_variant,
        transcript_id,
        orf_offset))

    # translate the variant cDNA sequence we detected from spanning reads
    # using the ORF offset from the current reference transcript
    combined_variant_cdna_sequence = DNA(
        cdna_sequence_before_variant + cdna_sequence_of_variant + cdna_sequence_after_variant)

    variant_protein_fragment_sequence = str(
        combined_variant_cdna_sequence[orf_offset:].translate())

    logger.info("Combined variant cDNA sequence %s translates to protein %s" % (
        combined_variant_cdna_sequence,
        variant_protein_fragment_sequence))

    # translate the reference sequence using the given ORF offset,
    # we can probably sanity check this by making sure it matches a
    # substring of the transcript.protein_sequence field
    combined_transcript_sequence = DNA(
        transcript_full_cdna[
            query_sequence_start_idx:
            query_sequence_start_idx + len(combined_variant_cdna_sequence)])

    transcript_protein_fragment_sequence = str(
        combined_transcript_sequence[orf_offset].translate())

    fragment_aa_start_offset_in_protein = (
        query_sequence_start_idx - start_codon_idx) // 3
    fragment_aa_end_offset_in_protein = (
        query_sequence_end_idx - start_codon_idx) // 3 + 1

    n_prefix_codons = n_prefix_nucleotides // 3

    if variant_is_deletion and variant_is_frameshift:
        # after a deletion that shifts the reading frame,
        # there is a partial codon left, which may be different from
        # the original codon at that location. Additionally,
        # the reading frame for all subsequent codons will be different
        # from the reference
        variant_codons_start_offset_in_fragment = None

    # TODO: change variant_codons_start_offset_in_fragment in the
    # TranslationFromReferenceORF to
    # - variant_amino_acids_start_offset_in_fragment
    # - variant_amino_acids_end_offset_in_fragment
    # which are computed from _codons_ offsets by checking which amino
    # acids differ from the reference sequence(s).
    # If a mutated sequence ends up having no differences then
    # we should skip it as a synonymous mutation.

    # TODO #2:
    # Instead of doing this once per sequence/transcript pair, we
    # should precompute the reference transcript sequences and the
    # ORFs they imply and pass in a list of ReferenceTranscriptContext
    # objects with the following field:
    #   sequence_before_variant_locus : cDNA
    #   reading_frame_at_start_of_sequence : int
    #   transcript_ids : str set
    #   transcript_names : str set
    #   gene : str

    variant_aa_start_offset_in_fragment = 0
    last_variant_codon = 0

    """
    # the number of non-mutated codons in the prefix (before the variant)
    # has to trim the ORF offset and then count up by multiples of 3


    aa_prefix = variant_protein_fragment_sequence[:n_prefix_codons]
    aa_variant = variant_protein_fragment_sequence[
        n_prefix_codons:n_prefix_codons + n_variant_codons]
    aa_suffix = variant_protein_fragment_sequence[
        n_prefix_codons + n_variant_codons:]

    assert aa_prefix + aa_variant + aa_suffix == variant_protein_fragment_sequence
                    aa_prefix=aa_prefix,
            aa_variant=aa_variant,
            aa_suffix=aa_suffix)
    """



def translate_variant_on_transcript(
        dna_sequence_prefix,
        dna_sequence_variant,
        dna_sequence_suffix,
        variant_is_insertion,
        variant_is_deletion,
        variant_is_frameshift,
        base1_variant_start_location,
        base1_variant_end_location,
        transcript,
        max_transcript_mismatches=MAX_TRANSCRIPT_MISMATCHES,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH):
    if transcript.strand == "+":
        cdna_prefix = dna_sequence_prefix
        cdna_suffix = dna_sequence_suffix
        cdna_variant = dna_sequence_variant
    else:
        # if the transcript is on the reverse strand then we have to
        # take the sequence PREFIX|VARIANT|SUFFIX
        # and take the complement of XIFFUS|TNAIRAV|XIFERP
        cdna_prefix = str(DNA(dna_sequence_suffix).reverse_complement())
        cdna_suffix = str(DNA(dna_sequence_prefix).reverse_complement())
        cdna_variant = str(
            DNA(dna_sequence_variant).reverse_complement())

    n_mismatch_before_variant = sum(
        xi != yi
        for (xi, yi) in zip(
            transcript_sequence_before_variant, cdna_prefix))

    if n_mismatch_before_variant > max_transcript_mismatches:
        logger.info(
            "Skipping transcript %s, too many mismatching bases (%d)",
            transcript,
            n_mismatch_before_variant)
        return None

    transcript_reading_frame = (query_sequence_start_idx - start_codon_idx) % 3

    return TranslationFromReferenceORF(
        cdna_prefix=cdna_prefix,
        cdna_variant=cdna_variant,
        cdna_suffix=cdna_suffix,
        transcript_id=transcript.id,
        transcript_name=transcript.name,
        number_transcript_sequence_mismatches=n_mismatch,
        reading_frame_at_start_of_cdna_sequence=reading_frame,
        transcript_sequence_before_variant=transcript_sequence_before_variant,
        variant_protein_sequence=variant_protein_fragment_sequence,
        reference_protein_sequence=transcript_protein_fragment_sequence,
        fragment_aa_start_offset_in_protein=fragment_aa_start_offset_in_protein,
        fragment_aa_end_offset_in_protein=fragment_aa_end_offset_in_protein,
        variant_aa_start_offset_in_fragment=variant_aa_start_offset_in_fragment,
        variant_aa_end_offset_in_fragment=variant_aa_end_offset_in_fragment)

def translate_variant_sequence(
        variant_sequence,
        reference_context,
        max_transcript_mismatches):
    """
    Attempt to translate a single VariantSequence using the reading frame
    from a single ReferenceContext.

    Parameters
    ----------
    variant_sequence : VariantSequence

    reference_context : ReferenceContext

    max_transcript_mismatches : int
        Don't use the reading frame from a context where the cDNA variant
        sequences disagrees at more than this number of positions before the
        variant nucleotides.

    Returns either a ProteinSequence object or None.
    """
    pass


def translate_multiple_protein_sequences_from_variant_sequences(
        variant_sequences,
        reference_contexts,
        max_transcript_mismatches,
        max_protein_sequences_per_variant):
    """
    Given all detected VariantSequence objects for a particular variant
    and all the ReferenceContext objects for that locus, translate
    multiple protein sequences, up to the number specified by the argument
    max_protein_sequences_per_variant.

    Parameters
    ----------
    variant_sequences : list of VariantSequence objects
        Variant sequences overlapping a single original variant

    reference_contexts : list of ReferenceContext objects
        Reference sequence contexts from the same variant as the variant_sequences

    max_transcript_mismatches : int

    max_protein_sequences_per_variant
    """
    protein_sequences = []
    for reference_context in reference_contexts:
        for variant_sequence in variant_sequences:
            # assuming that the reference contexts and variant sequences
            # were sorted in order of priority, we will also get back
            # translations in order of quality, so it's OK to only keep
            # the top k
            n_protein_sequences = len(protein_sequences)
            if n_protein_sequences >= max_protein_sequences_per_variant:
                return protein_sequences
            protein_sequence = translate_variant_sequence(
                variant_sequence=variant_sequence,
                reference_context=reference_context,
                max_transcript_mismatches=max_transcript_mismatches)
            if protein_sequence is not None:
                protein_sequences.append(protein_sequence)
    return protein_sequences

def translate_variants(
        variants,
        samfile,
        transcript_id_whitelist=None,
        protein_fragment_length=PROTEIN_FRAGMENT_LEGNTH,
        min_reads_supporting_rna_sequence=MIN_READS_SUPPORTING_RNA_SEQUENCE,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
        max_transcript_mismatches=MAX_TRANSCRIPT_MISMATCHES,
        max_protein_sequences_per_variant=MAX_SEQUENCES_PER_VARIANT):
    """
    Translates each coding variant in a collection to one or more protein
    fragment sequences (if the variant is not filtered and its spanning RNA
    sequences can be given a reading frame).

    Parameters
    ----------
    variants : varcode.VariantCollection

    samfile : pysam.AlignmentFile

    transcript_id_whitelist : set, optional
        If given, expected to be a set of transcript IDs which we should use
        for determining the reading frame around a variant. If omitted, then
        try to use all overlapping reference transcripts.

    protein_fragment_length : int
        Try to translate protein sequences of this length, though sometimes
        we'll have to return something shorter (depending on the RNAseq data,
        and presence of stop codons).

    min_reads_supporting_rna_sequence : int
        Drop variant sequences supported by fewer than this number of reads.

    min_transcript_prefix_length : int
        Minimum number of bases we need to try matching between the reference
        context and variant sequence.

    max_transcript_mismatches : int
        Don't try to determine the reading frame for a transcript if more
        than this number of bases differ.

    max_protein_sequences_per_variant : int
        Ignore any translations once we have this many for a variant.

    Returns a dictionary mapping each variant to a list of ProteinSequence
    records.
    """

    # adding 2nt to total RNA sequence length  in case we need to clip 1 or 2
    # bases of the sequence to match a reference ORF but still want to end up
    # with the desired number of amino acids
    rna_sequence_length = protein_fragment_length * 3 + 2

    variant_to_protein_sequences_dict = OrderedDict()

    for variant, variant_sequences in variant_sequences_generator(
            variants=variants,
            samfile=samfile,
            sequence_length=rna_sequence_length,
            min_reads=min_reads_supporting_rna_sequence):
        # include every variant in the dictionary, even if it doesn't yield any
        # translated proteins, just in case we later want to count how many
        # translations succeeded vs. failed
        variant_to_protein_sequences_dict[variant] = []

        if len(variant_sequences) == 0:
            logger.info(
                "Skipping variant %s, no cDNA sequences detected" % (
                    variant,))
            continue

        # try translating the variant sequences from the same set of
        # ReferenceContext objects, which requires using the longest
        # context_size to be compatible with all of the sequences. Some
        # sequences maybe have fewer nucleotides than this before the variant
        # and will thus have to be trimmed.
        context_size = max(
            len(variant_sequence.prefix)
            for variant_sequence in variant_sequences)

        if context_size < min_transcript_prefix_length:
            logger.info(
                "Skipping variant %s, none of the cDNA sequences have sufficient context" % (
                    variant,))
            continue

        reference_contexts = reference_contexts_for_variant(
            variant,
            context_size=context_size,
            transcript_id_whitelist=transcript_id_whitelist)
        variant_to_protein_sequences_dict[variant] = \
            translate_multiple_protein_sequences_from_variant_sequences(
                variant_sequences=variant_sequences,
                reference_contexts=reference_contexts,
                max_transcript_mismatches=max_transcript_mismatches,
                max_protein_sequences_per_variant=max_protein_sequences_per_variant)
    return variant_to_protein_sequences_dict

def translate_variants_dataframe(
        variants,
        samfile,
        transcript_id_whitelist=None,
        protein_fragment_length=PROTEIN_FRAGMENT_LEGNTH,
        min_reads_supporting_rna_sequence=MIN_READS_SUPPORTING_RNA_SEQUENCE,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
        max_transcript_mismatches=MAX_TRANSCRIPT_MISMATCHES,
        max_sequences_per_variant=MAX_SEQUENCES_PER_VARIANT):
    """
    Given a collection of variants and a SAM/BAM file of overlapping reads,
    returns a DataFrame of translated protein fragments with the following
    columns:
        chr : str
            Chromosome of variant

        base1_pos : int
            First reference nucleotide affected by variant (or position before
            insertion)

        ref : str
            Reference nucleotides

        alt : str
            Variant nucleotides

        cdna_sequence : str
            cDNA sequence context detected from RNAseq BAM

        cdna_mutation_start_offset : int
            Interbase start offset for variant nucleotides in the cDNA sequence

        cdna_mutation_end_offset : int
            Interbase end offset for variant nucleotides in the cDNA sequence

        supporting_read_count : int
            How many reads fully spanned the cDNA sequence

        variant_protein_sequence : str
            Translated variant protein fragment sequence

        variant_protein_sequence_length : int
            Number of amino acids in each sequence

        reference_transcript_ids : str list

        reference_transcript_names : str list

        reference_protein_sequences : str list

        cdna_sequence_to_support_reads : dict
            Maps distinct (prefix, variant, suffix) triplets to
            names of reads supporting these cDNA sequences

        total_supporting_read_count : int

        """
    # construct a dictionary incrementally which we'll turn into a
    # DataFrame
    columns = [
        # fields related to variant
        ("chr", []),
        ("pos", []),
        ("ref", []),
        ("alt", []),
    ]

    # flattened representation of fields in ProteinSequence, and its member
    # objects of type ReferenceContext and VariantSequence
    for field in ProteinSequence._fields:
        if field != "reference_context" and field != "variant_sequence":
            columns.append((field, []))

    for field in ReferenceContext._fields + VariantSequence._fields:
        columns.append((field, []))

    column_dict = OrderedDict(columns)

    for (variant, protein_sequences) in translate_variants(
            variants,
            samfile,
            transcript_id_whitelist=transcript_id_whitelist,
            protein_fragment_length=protein_fragment_length,
            min_reads_supporting_rna_sequence=min_reads_supporting_rna_sequence,
            min_transcript_prefix_length=min_transcript_prefix_length,
            max_transcript_mismatches=max_transcript_mismatches,
            max_sequences_per_variant=max_sequences_per_variant).items():
        for protein_sequence in protein_sequences:
            # variant info
            column_dict["chr"].append(variant.contig)
            column_dict["pos"].append(variant.original_start)
            column_dict["ref"].append(variant.original_ref)
            column_dict["alt"].append(variant.original_alt)
            for field in ProteinSequence._fields:
                if field == "variant_sequence" or field == "reference_context":
                    continue
                column_dict[field].append(
                    getattr(protein_sequence, field))
            # include fields of protein_sequence.reference_context directly
            # with the primary columns
            for field in ReferenceContext._fields:
                column_dict[field].append(
                    getattr(protein_sequence.reference_context, field))
            # include fields of protein_sequence.VariantSequence directly
            # with the primary columns
            for field in VariantSequence._fields:
                column_dict[field].append(
                    getattr(protein_sequence.variant_sequence, field))

    df = pd.DataFrame(column_dict)
    return df
