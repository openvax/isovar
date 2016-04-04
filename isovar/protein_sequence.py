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
    cdna_prefix = DNA(variant_sequence.prefix)
    cdna_alt = DNA(variant_sequence.alt)
    cdna_suffix = DNA(variant_sequence.suffix)

    # if the transcript is on the reverse strand then we have to
    # take the sequence PREFIX|VARIANT|SUFFIX
    # and take the complement of XIFFUS|TNAIRAV|XIFERP
    if reference_context.strand == "-":
        cdna_prefix = cdna_prefix.reverse_complement()
        cdna_alt = cdna_alt.reverse_complement()
        cdna_suffix = cdna_suffix.reverse_complement()

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
    assert len(reference_sequence_before_variant) == len(cdna_prefix), \
        "Something has gone wrong! %s should be same length as %s" % (
            reference_sequence_before_variant, cdna_prefix)

    n_mismatch_before_variant = sum(
        xi != yi
        for (xi, yi) in zip(
            reference_sequence_before_variant,
            cdna_prefix))

    if n_mismatch_before_variant > max_transcript_mismatches:
        logger.info(
            "Skipping reference context %s for %s, too many mismatching bases (%d)",
            reference_context,
            variant_sequence,
            n_mismatch_before_variant)
        return None

    # ReferenceContext carries with an offset to the first complete codon
    # in the reference sequence. This may need to be adjusted if the reference
    # sequence is longer than the variant sequence (and thus needs to be trimmed)
    offset_to_first_complete_codon = reference_context.offset_to_first_complete_codon

    # if the reference sequence is longer then add in the number of
    # of codons (full or partial) that we trimmed
    n_reference_codons_trimmed = int(math.ceil(n_trimmed_from_reference / 3.0))
    offset_to_first_complete_codon += n_reference_codons_trimmed * 3

    cdna_prefix_from_first_codon = cdna_prefix[offset_to_first_complete_codon:]
    combined_variant_cdna_sequence = DNA.concat([
        cdna_prefix_from_first_codon, cdna_alt, cdna_suffix])

    variant_amino_acids = combined_variant_cdna_sequence.translate()

    # rounding down since a change in the middle of a codon should count
    # toward the variant codons
    n_complete_prefix_codons = len(cdna_prefix_from_first_codon) // 3

    if len(cdna_alt) % 3 != 0:
        pass

    n_variant_codons = int(math.ceil(len(cdna_alt) / 3.0))

    return ProteinSequence(
        amino_acids=variant_amino_acids,
        variant_aa_interval_start=n_complete_prefix_codons + 1,
        variant_aa_interval_end=variant_aa_interval_end,
        reference_context=reference_context,
        variant_sequence=variant_sequence)


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
