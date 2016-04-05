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
import math
import logging

import pandas as pd
from skbio import DNA

from .reference_context import (
    ReferenceContext,
    reference_contexts_for_variant
)
from .variant_sequence import variant_sequences_generator, VariantSequence

# When multiple distinct RNA sequence contexts and/or reference transcripts
# give us the same translation then we group them into a single object
# which summarizes the supporting read count and reference transcript ORFs
# for a unique protein sequence.

################################
#
# VariantSequenceInReadingFrame
# ------------------------------
#
# Combines a VariantSequence with the reading frame implied by a
# ReferenceContext, reverse complementing if necessary and finding the
# offset to the first complete codon in the cDNA sequence.
#
#################################


VariantSequenceInReadingFrame = namedtuple(
    "VariantSequenceInReadingFrame",
    (
        # since the reference context and variant sequence may have
        # different numbers of nucleotides before the variant, the cDNA prefix
        # gets truncated to the shortest length. To avoid having to recompute
        # that sequence again, let's just cache the full cDNA sequence we used
        # for translation here, along with an interval indicating which
        # nucleotides are from the variant of interest
        "cdna_sequence",
        "offset_to_first_complete_codon",
        "variant_cdna_interval_start",
        "variant_cdna_interval_end",
        "reference_cdna_sequence_before_variant",
        "number_mismatches",
    ))


##########################
#
# ProteinSequence
# ---------------
#
# Translation of a VariantSequenceInReadingFrame to amino acids.
#
##########################


ProteinSequence = namedtuple(
    "ProteinSequence", VariantSequenceInReadingFrame._fields + (
        # translated sequence of a variant sequence in the ORF established
        # by a reference context
        "amino_acids",
        # half-open interval coordinates for variant amino acids
        # in the translated sequence
        "variant_aa_interval_start",
        "variant_aa_interval_end",
        # did the amino acid sequence end due to a stop codon or did we
        # just run out of sequence context around the variant?
        "ends_with_stop_codon",
        # was the variant a frameshift relative to the reference sequence?
        "frameshift",
    ))


MIN_READS_SUPPORTING_RNA_SEQUENCE = 3
MIN_TRANSCRIPT_PREFIX_LENGTH = 30
MAX_TRANSCRIPT_MISMATCHES = 2
PROTEIN_FRAGMENT_LEGNTH = 20
MAX_SEQUENCES_PER_VARIANT = 1

def trim_sequences(variant_sequence, reference_context):
    """
    A VariantSequence and ReferenceContext may contain a different number of
    nucleotides before the variant locus. Furthermore, the VariantSequence is
    always expressed in terms of the positive strand against which it aligned,
    but reference transcripts may have sequences from the negative strand of the
    genome. Take the reverse complement of the VariantSequence if the
    ReferenceContext is from negative strand transcripts and trim either
    sequence to ensure that the prefixes are of the same length.

    Parameters
    ----------
    variant_sequence : VariantSequence

    reference_context : ReferenceContext

    Returns a tuple with the following fields:
        1) cDNA prefix of variant sequence, trimmed to be same length as the
           reference prefix. If the reference context was on the negative
           strand then this is the trimmed sequence *after* the variant from
           the genomic DNA sequence.
        2) cDNA sequence of the variant nucleotides, in reverse complement if
           the reference context is on the negative strand.
        3) cDNA sequence of the nucleotides after the variant nucleotides. If
           the reference context is on the negative strand then this sequence
           is the reverse complement of the original prefix sequence.
        4) Reference sequence before the variant locus, trimmed to be the
           same length as the variant prefix.
        5) Number of nucleotides trimmed from the reference sequence, used
           later for adjustint offset to first complete codon.
    """
    cdna_prefix = DNA(variant_sequence.prefix)
    cdna_alt = DNA(variant_sequence.alt)
    cdna_suffix = DNA(variant_sequence.suffix)

    # if the transcript is on the reverse strand then we have to
    # take the sequence PREFIX|VARIANT|SUFFIX
    # and take the complement of XIFFUS|TNAIRAV|XIFERP
    if reference_context.strand == "-":
        cdna_prefix = cdna_suffix.reverse_complement()
        cdna_alt = cdna_alt.reverse_complement()
        cdna_suffix = cdna_prefix.reverse_complement()

    reference_sequence_before_variant = reference_context.sequence_before_variant_locus

    # trim the reference sequence and the RNA-derived sequence to the same length
    if len(reference_sequence_before_variant) > len(cdna_prefix):
        n_trimmed_from_reference = len(reference_sequence_before_variant) - len(cdna_prefix)
        n_trimmed_from_variant = 0
    elif len(reference_sequence_before_variant) < len(cdna_prefix):
        n_trimmed_from_variant = len(cdna_prefix) - len(reference_sequence_before_variant)
        n_trimmed_from_reference = 0
    else:
        n_trimmed_from_variant = 0
        n_trimmed_from_reference = 0

    reference_sequence_before_variant = reference_sequence_before_variant[
        n_trimmed_from_reference:]

    cdna_prefix = cdna_prefix[n_trimmed_from_variant:]

    return (
        cdna_prefix,
        cdna_alt,
        cdna_suffix,
        reference_sequence_before_variant,
        n_trimmed_from_reference
    )


def count_mismatches(reference_prefix, cdna_prefix):
    """
    Computes the number of mismatching nucleotides between two cDNA sequences.

    Parameters
    ----------
    reference_prefix : str or skbio.DNA
        cDNA sequence of a reference transcript before a variant locus

    cdna_prefix : str or skbio.DNA
        cDNA sequence detected from RNAseq before a variant locus
    """
    if len(reference_prefix) != len(cdna_prefix):
        raise ValueError(
            "Expected reference prefix '%s' to be same length as %s" % (
                reference_prefix, cdna_prefix))
    # converting to str before zip because skbio.DNA is weird and sometimes
    # returns elements which don't compare equally despite being the same
    # nucleotide
    return sum(xi != yi for (xi, yi) in zip(
        str(reference_prefix), str(cdna_prefix)))


def compute_offset_to_first_complete_codon(
        offset_to_first_complete_reference_codon,
        n_trimmed_from_reference):
    """
    Once we've aligned the variant sequence to the ReferenceContext, we need
    to transfer reading frame from the reference transcripts to the variant
    sequences.

    Parameters
    ----------
    offset_to_first_complete_codon_in_reference_context : ReferenceContext

    n_trimmed_from_reference : int

    Returns an offset into the variant sequence that starts from a complete
    codon.
    """
    # if the reference sequence is longer then add in the number of
    # of codons (full or partial) that we trimmed
    n_reference_codons_trimmed = int(math.ceil(n_trimmed_from_reference / 3.0))
    n_reference_nucleotides_trimmed = n_reference_codons_trimmed * 3
    return offset_to_first_complete_reference_codon + n_reference_nucleotides_trimmed

def align_variant_sequence_to_reference_context(
        variant_sequence,
        reference_context,
        max_transcript_mismatches=MAX_TRANSCRIPT_MISMATCHES):
    """
    Parameters
    ----------
    variant_sequence : VariantSequence

    reference_context : ReferenceContext

    Returns a VariantSequenceInReadingFrame object or, if the max number of
    mismatches is exceeded, then None.
    """
    cdna_prefix, cdna_alt, cdna_suffix, reference_prefix, n_trimmed_from_reference = \
        trim_sequences(variant_sequence, reference_context)

    n_mismatch_before_variant = count_mismatches(reference_prefix, cdna_prefix)

    if n_mismatch_before_variant > max_transcript_mismatches:
        logging.info(
            "Skipping reference context %s for %s, too many mismatching bases (%d)",
            reference_context,
            variant_sequence,
            n_mismatch_before_variant)
        return None

    # ReferenceContext carries with an offset to the first complete codon
    # in the reference sequence. This may need to be adjusted if the reference
    # sequence is longer than the variant sequence (and thus needs to be trimmed)
    offset_to_first_complete_codon = compute_offset_to_first_complete_codon(
        offset_to_first_complete_reference_codon=reference_context.offset_to_first_complete_codon,
        n_trimmed_from_reference=n_trimmed_from_reference)

    cdna_sequence = DNA.concat([cdna_prefix, cdna_alt, cdna_suffix])
    variant_interval_start = len(cdna_prefix) + 1
    variant_interval_end = variant_interval_start + len(cdna_alt)

    return VariantSequenceInReadingFrame(
        cdna_sequence=cdna_sequence,
        offset_to_first_complete_codon=offset_to_first_complete_codon,
        variant_cdna_interval_start=variant_interval_start,
        variant_cdna_interval_end=variant_interval_end,
        reference_cdna_sequence_before_variant=reference_prefix,
        number_mismatches=n_mismatch_before_variant)

def find_mutant_amino_acid_interval(
        cdna_sequence,
        cdna_first_codon_offset,
        cdna_variant_start_offset,
        cdna_variant_end_offset,
        n_ref,
        n_amino_acids):
    """
    Parameters
    ----------
    cdna_sequence : skbio.DNA or str
        cDNA sequence found in RNAseq data

    cdna_first_codon_offset : int
        Offset into cDNA sequence to first complete codon, lets us skip
        past UTR region and incomplete codons.

    cdna_variant_start_offset : int
        Interbase start offset into cDNA sequence for selecting mutant
        nucleotides.

    cdna_variant_end_offset : int
        Interbase end offset into cDNA sequence for selecting mutant
        nucleotides.

    n_ref : int
        Number of reference nucleotides

    n_amino_acids : int
        Number of translated amino acids

    Returns tuple with three fields:
        1) Start offset for interval of mutant amino acids in translated sequence
        2) End offset for interval of mutant amino acids in translated sequence
        3) Boolean flag indicating whether the variant was a frameshift.
    """
    cdna_alt_nucleotides = cdna_sequence[
        cdna_variant_start_offset:cdna_variant_end_offset]

    n_alt = len(cdna_alt_nucleotides)

    # sequence of nucleotides before the variant starting from the first codon
    cdna_coding_prefix = cdna_sequence[cdna_first_codon_offset:cdna_variant_start_offset]

    # rounding down since a change in the middle of a codon should count
    # toward the variant codons
    n_coding_nucleotides_before_variant = len(cdna_coding_prefix)
    n_complete_prefix_codons = cdna_variant_start_offset // 3

    frame_of_variant_nucleotides = n_coding_nucleotides_before_variant % 3
    frameshift = abs(n_ref - n_alt) % 3 != 0
    indel = n_ref != n_alt

    variant_aa_interval_start = n_complete_prefix_codons + 1
    if frameshift:
        # if mutation is a frame shift then every amino acid from the
        # first affected codon to the stop is considered mutant
        #
        # TODO: what if the first k amino acids are synonymous with the
        # reference sequence?
        variant_aa_interval_end = n_amino_acids
    else:
        n_alt_codons = int(math.ceil(n_alt / 3.0))
        if indel:
            # We need to adjust the number of affected codons by whether the
            # variant is aligned with codon boundaries, since in-frame indels
            # may still be split across multiple codons.
            #
            # Example of in-frame deletion of 3 nucleotides which leaves
            # 0 variant codons in the sequence (interval = 1:1)
            #   ref = CCC|AAA|GGG|TTT
            #   alt = CCC|GGG|TTT
            #
            # Example of in-frame deletion of 3 nucleotides which leaves
            # 1 variant codon in the sequence (interval = 1:2)
            #   ref = CCC|AAA|GGG|TTT
            #   alt = CCC|AGG|TTT
            #
            # Example of in-frame insertion of 3 nucleotides which
            # yields two variant codons:
            #   ref = CCC|AAA|GGG|TTT
            #   alt = CTT|TCC|AAA|GGG|TTT
            extra_affected_codon = int(frame_of_variant_nucleotides != 0)
            variant_aa_interval_end = (
                variant_aa_interval_start + n_alt_codons + extra_affected_codon)
        else:
            # if the variant is a simple substitution then it only affects
            # as many codons as are in the alternate sequence
            variant_aa_interval_end = variant_aa_interval_start + n_alt_codons
    return variant_aa_interval_start, variant_aa_interval_end, frameshift

def translate_cdna(cdna_sequence_from_first_codon):
    """
    Given a cDNA sequence which is aligned to a reading frame, returns
    the translated protein sequence and a boolean flag indicating whether
    the translated sequence ended on a stop codon (or just ran out of codons).

    Parameters
    ----------
    cdna_sequence_from_first_codon : skbio.DNA
        cDNA sequence which is expected to start and end on complete codons.
    """
    # once we drop some of the prefix nucleotides, we should be in a reading frame
    # which allows us to translate this protein
    variant_amino_acids = cdna_sequence_from_first_codon.translate()

    # by default, the translate method emits "*" for
    # the stop codons TAA, TAG and TGG
    ends_with_stop_codon = "*" in variant_amino_acids

    if ends_with_stop_codon:
        # truncate the amino acid sequence at the stop codon
        stop_codon_index = variant_amino_acids.index("*")
        variant_amino_acids = variant_amino_acids[:stop_codon_index]

    return variant_amino_acids, ends_with_stop_codon


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

    variant_sequence_in_reading_frame = align_variant_sequence_to_reference_context(
        variant_sequence=variant_sequence,
        reference_context=reference_context,
        max_transcript_mismatches=max_transcript_mismatches)

    if variant_sequence_in_reading_frame is None:
        logging.debug("Failed to align %s to reference context %s" % (
            variant_sequence,
            reference_context))
        return None

    cdna_sequence = variant_sequence_in_reading_frame.cdna_sequence
    cdna_codon_offset = variant_sequence_in_reading_frame.offset_to_first_complete_codon
    cdna_sequence_from_first_codon = cdna_sequence[cdna_codon_offset:]

    variant_amino_acids, ends_with_stop_codon = translate_cdna(cdna_sequence_from_first_codon)

    # get the offsets into the cDNA sequence which pick out the variant nucleotides
    cdna_variant_start_offset = variant_sequence_in_reading_frame.variant_cdna_interval_start
    cdna_variant_end_offset = variant_sequence_in_reading_frame.variant_cdna_interval_end

    variant_aa_interval_start, variant_aa_interval_end, frameshift = \
        find_mutant_amino_acid_interval(
            cdna_sequence=cdna_sequence,
            cdna_first_codon_offset=cdna_codon_offset,
            cdna_variant_start_offset=cdna_variant_start_offset,
            cdna_variant_end_offset=cdna_variant_end_offset,
            n_ref=len(reference_context.sequence_at_variant_locus),
            n_amino_acids=len(variant_amino_acids))

    reference_sequence_before_variant = (
        variant_sequence_in_reading_frame.reference_cdna_sequence_before_variant)

    return ProteinSequence(
        cdna_sequence=cdna_sequence,
        offset_to_first_complete_codon=cdna_codon_offset,
        variant_cdna_interval_start=cdna_variant_start_offset,
        variant_cdna_interval_end=cdna_variant_end_offset,
        reference_cdna_sequence_before_variant=reference_sequence_before_variant,
        number_mismatches=variant_sequence_in_reading_frame.number_mismatches,
        amino_acids=variant_amino_acids,
        frameshift=frameshift,
        ends_with_stop_codon=ends_with_stop_codon,
        variant_aa_interval_start=variant_aa_interval_start,
        variant_aa_interval_end=variant_aa_interval_end)


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
            logging.info(
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
            logging.info(
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
