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

import logging
from collections import namedtuple, defaultdict, OrderedDict

from skbio import DNA
import numpy as np
import pandas as pd

from .variant_reads import gather_variant_reads
from .sequence_counts import sequence_counts

# information related to the translation of a RNA sequence context
# in a reading frame determined by a particular reference transcript
TranslationFromReferenceORF = namedtuple(
    "TranslationFromReferenceORF",
    [
        "cdna_prefix",
        "cdna_variant",
        "cdna_suffix",
        "transcript_id",
        "transcript_name",
        "number_transcript_sequence_mismatches",
        "fraction_transcript_sequence_mismatches",
        "reading_frame_at_start_of_cdna_sequence",
        "transcript_sequence_before_variant",
        "variant_protein_sequence",
        "reference_protein_sequence",
        "protein_fragment_start_offset",
        "protein_variant_offset",
    ])

# if multiple distinct RNA sequence contexts and/or reference transcripts
# give us the same translation then we group them into a single object
# which summarizes the supporting read count and reference transcript ORFs
# for a unique protein sequence
ProteinFragment = namedtuple(
    "ProteinFragment",
    [
        # number of reads supporting any RNA sequence which translates to
        # this protein sequence
        "number_supporting_reads",
        "variant_protein_sequence",
        "reference_transcript_ids",
        "reference_transcript_names",
        "reference_protein_sequences",
        "cdna_sequences_to_read_names",
    ])

MIN_READS_SUPPORTING_RNA_SEQUENCE = 3
MIN_TRANSCRIPT_PREFIX_LENGTH = 15
MAX_TRANSCRIPT_MISMATCHES = 2
PROTEIN_FRAGMENT_LEGNTH = 25
MAX_SEQUENCES_PER_VARIANT = 5


def translate_compatible_reading_frames(
        dna_sequence_prefix,
        dna_sequence_variant,
        dna_sequence_suffix,
        variant_is_insertion,
        base1_variant_start_location,
        base1_variant_end_location,
        transcripts,
        max_transcript_mismatches=MAX_TRANSCRIPT_MISMATCHES,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH):
    """
    Use the given reference transcripts to attempt to establish the ORF
    of a sequence context extract from RNA reads. It's expected that the
    sequence has been aligned to the reference genome potentially using its
    reverse-complement, and thus the sequence we are given is always from the
    positive strand of DNA.

    Parameters
    ----------
    dna_sequence_prefix : str
        Nucleotides before the variant.

    dna_sequence_variant : str
        Mutated nucleotides, should be empty for a deletion.

    dna_sequence_suffix : str
        Nucleotides after the variant

    variant_is_insertion : bool

    base1_variant_start_location : int
        For deletions and substitutions, this position is the first modified
        nucleotide. For insertions, this is the position before any inserted
        nucleotides.

    base1_variant_end_location : int
        For deletions and substitutions, this is the positions of the last
        affected reference nucleotide. For insertions, this is the location of
        the base after the insertion.

    transcripts : list of pyensembl.Transcript
        List of candidate reference transcripts from which we try to determine
        the ORF of a variant RNA sequence.

    max_transcript_mismatches : int
        Ignore transcripts with more than this number of mismatching nucleotides

    min_transcript_prefix_length : int
        Don't consider transcripts with less than this number of nucleotides
        before the variant position (setting this value to 0 will enable use
        of transcripts without 5' UTR)

    Returns list of TranslationFromReferenceORF objects.
    """
    results = []

    for transcript in transcripts:
        if transcript.strand == "+":
            cdna_prefix = dna_sequence_prefix
            cdna_suffix = dna_sequence_suffix
            cdna_variant = dna_sequence_variant
            variant_in_transcript_idx = transcript.spliced_offset(
                base1_variant_start_location)
        else:
            # if the transcript is on the reverse strand then we have to
            # take the sequence PREFIX|VARIANT|SUFFIX
            # and take the complement of XIFFUS|TNAIRAV|XIFERP
            cdna_prefix = str(DNA(dna_sequence_suffix).reverse_complement())
            cdna_suffix = str(DNA(dna_sequence_prefix).reverse_complement())
            cdna_variant = str(
                DNA(dna_sequence_variant).reverse_complement())
            variant_in_transcript_idx = transcript.spliced_offset(
                base1_variant_end_location)

        if variant_is_insertion:
            # insertions don't actually affect the base referred to
            # by the start position of the variant, but rather the
            # variant gets inserted *after* that position
            query_sequence_end_idx = variant_in_transcript_idx + 1
        else:
            query_sequence_end_idx = variant_in_transcript_idx

        start_codon_idx = min(transcript.start_codon_spliced_offsets)

        if variant_in_transcript_idx < start_codon_idx + 3:
            logging.info(
                "Skipping %s because variant appears in 5' UTR" % (
                    transcript,))

        query_sequence_start_idx = query_sequence_end_idx - len(cdna_prefix)

        if query_sequence_start_idx < 0:
            logging.warn("Transcript %s not long enough for observed sequence" % (
                transcript,))
            continue
        transcript_sequence_before_variant = transcript.sequence[
            query_sequence_start_idx:query_sequence_end_idx]

        assert len(transcript_sequence_before_variant) == len(cdna_prefix)
        if len(transcript_sequence_before_variant) < min_transcript_prefix_length:
            continue

        n_mismatch = sum(
            xi != yi
            for (xi, yi) in zip(
                transcript_sequence_before_variant, cdna_prefix))
        if n_mismatch > max_transcript_mismatches:
            continue

        fraction_mismatch = float(n_mismatch) / len(cdna_prefix)

        reading_frame = (query_sequence_start_idx - start_codon_idx) % 3

        if reading_frame == 1:
            # if we're 1 nucleotide into the codon then we need to shift
            # over two more to restore the ORF
            orf_offset = 2
        elif reading_frame == 2:
            orf_offset = 1
        else:
            orf_offset = 0

        combined_variant_cdna_sequence = DNA(
            cdna_prefix + cdna_variant + cdna_suffix)

        in_frame_combined_variant_sequence = combined_variant_cdna_sequence[
            orf_offset:]
        variant_protein_fragment = in_frame_combined_variant_sequence.translate()

        combined_transcript_sequence = DNA(
            transcript.sequence[
                query_sequence_start_idx:
                query_sequence_start_idx + len(combined_variant_cdna_sequence)])
        in_frame_transcript_sequence = combined_transcript_sequence[
            orf_offset:]
        transcript_protein_sequence = in_frame_transcript_sequence.translate()
        protein_fragment_start_offset = (
            query_sequence_start_idx - start_codon_idx) // 3
        protein_variant_offset = (
            query_sequence_end_idx - start_codon_idx) // 3
        results.append(
            TranslationFromReferenceORF(
                cdna_prefix=cdna_prefix,
                cdna_variant=cdna_variant,
                cdna_suffix=cdna_suffix,
                transcript_id=transcript.id,
                transcript_name=transcript.name,
                number_transcript_sequence_mismatches=n_mismatch,
                fraction_transcript_sequence_mismatches=fraction_mismatch,
                reading_frame_at_start_of_cdna_sequence=reading_frame,
                transcript_sequence_before_variant=transcript_sequence_before_variant,
                variant_protein_sequence=str(variant_protein_fragment),
                reference_protein_sequence=str(transcript_protein_sequence),
                protein_fragment_start_offset=protein_fragment_start_offset,
                protein_variant_offset=protein_variant_offset // 3))
    return results

def rna_sequence_key(fragment_info):
    """
    Get the cDNA (prefix, variant, suffix) fields from a
    TranslationFromReferenceORF object.
    """
    return (
        fragment_info.cdna_prefix,
        fragment_info.cdna_variant,
        fragment_info.cdna_suffix
    )

def group_protein_fragments(
        protein_fragment_and_read_names_list,
        desired_protein_length=PROTEIN_FRAGMENT_LEGNTH):
    """
    If we end up with multiple equivalent protein fragments from distinct
    cDNA sequences then we want to group them since ultimately it's
    not relevant which codon gives us a particular amino acid.

    Parameters
    ----------
    fragments_and_read_names : list
        List of tuples containing (1) a TranslationFromReferenceORF object and
        (2) a set of read names from spanning RNA reads for that
        unique sequence.

    Returns list of ProteinFragment objects.
    """
    # map every distinct protein sequence to a tuple with the following fields:
    # - transcript ID
    # - transcript name
    # - cDNA sequence
    # - reading frame

    protein_sequence_dict = defaultdict(list)
    cdna_to_read_names = {}
    for (translation_info, read_names) in protein_fragment_and_read_names_list:
        cnda_tuple = rna_sequence_key(translation_info)
        if cnda_tuple in cdna_to_read_names:
            assert cdna_to_read_names[cnda_tuple] == read_names
        else:
            cdna_to_read_names[cnda_tuple] = read_names

        protein_sequence = translation_info.variant_protein_sequence
        # if we can chop up the translated protein fragment into shorter pieces
        # then go for it!

        variant_protein_sequence = translation_info.variant_protein_sequence
        ref_protein_sequence = translation_info.reference_protein_sequence
        n_amino_acids = len(variant_protein_sequence)

        for start, end in enumerate(range(desired_protein_length, n_amino_acids)):

            # all protein fragments must overlap the variant
            if start > translation_info.base0_variant_amino_acid_start_offset:
                break
            if end < translation_info.base0_variant_amino_acid_start_offset:
                continue

            variant_protein_subsequence = protein_sequence[start:end]
            ref_protein_subsequence = ref_protein_sequence[start:end]
            info_tuple = (
                variant_protein_subsequence,
                ref_protein_subsequence,
                translation_info.transcript_id,
                translation_info.transcript_name,
                # TODO: actually add these fields!
                translation_info.base0_variant_amino_acid_start_offset,
                translation_info.base0_variant_amino_acid_end_offset,
            )
            protein_sequence_dict[protein_sequence].append(info_tuple)

    results = []
    for (protein_sequence, info_objs) in protein_sequence_dict.items():
        transcript_names = list(set([x.transcript_name for x in info_objs]))
        transcript_ids = list(set([x.transcript_id for x in info_objs]))
        cdna_sequence_keys = list(set([rna_sequence_key(x) for x in info_objs]))

        reference_protein_sequences = []

        total_read_count = sum(
            len(cdna_to_read_names[cnda_tuple])
            for cnda_tuple in cdna_sequence_keys)
        transcript_ids = list(set([x.transcript_id for x in info_objs]))
        transcript_names = list(set([x.transcript_name for x in info_objs]))
        reference_protein_sequences = list(set([
            x.reference_protein_sequence
            for x in info_objs]))
        cdna_sequences_to_read_names = {
            rna_sequence_key(x): cdna_to_read_names[rna_sequence_key(x)]
            for x in info_objs
        }
        list(set(
            [rna_sequence_key(x) for x in info_objs]))
        results.append(
            ProteinFragment(
                number_supporting_reads=total_read_count,
                variant_protein_sequence=protein_sequence,
                reference_transcript_ids=transcript_ids,
                reference_transcript_names=transcript_names,
                reference_protein_sequences=reference_protein_sequences,
                cdna_sequences_to_read_names=cdna_sequences_to_read_names))
    return results

def variant_protein_fragments_with_read_counts(
        variant,
        samfile,
        reference_transcripts,
        protein_fragment_length=PROTEIN_FRAGMENT_LEGNTH,
        min_reads_supporting_rna_sequence=MIN_READS_SUPPORTING_RNA_SEQUENCE,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
        max_transcript_mismatches=MAX_TRANSCRIPT_MISMATCHES,
        max_sequences=MAX_SEQUENCES_PER_VARIANT,
        chromosome_name=None):
    """
    Returns list of tuples containing:
        1) ProteinFragment object
        2) set of read names supporting the protein fragment
    """

    if len(reference_transcripts) == 0:
        return []

    variant_reads = gather_variant_reads(
        samfile=samfile,
        chromosome=chromosome_name if chromosome_name else variant.contig,
        base1_location=variant.start,
        ref=variant.ref,
        alt=variant.alt)
    if len(variant_reads) < min_reads_supporting_rna_sequence:
        return []

    # adding 2nt to total RNA sequence length  in case we need to clip 1 or 2
    # bases of the sequence to match a reference ORF but still want to end up
    # with the desired number of amino acids
    rna_sequence_length = protein_fragment_length * 3 + 2

    # the number of context nucleotides on either side of the variant
    # is half the desired length (minus the number of variant nucleotides)
    n_surrounding_nucleotides = rna_sequence_length - len(variant.alt)

    flanking_context_size = int(np.ceil(n_surrounding_nucleotides / 2.0))
    logging.info(
        "Looking or %dnt RNA sequence around %s for %d codons)" % (
            rna_sequence_length,
            variant,
            protein_fragment_length))

    sequence_count_info = sequence_counts(
        variant_reads,
        context_size=flanking_context_size)

    sequences_and_read_names = sequence_count_info.full_read_names

    if len(sequences_and_read_names) == 0:
        return []

    variant_seq = sequence_count_info.variant_nucleotides

    result_list = []

    for i, ((prefix, suffix), read_names) in enumerate(sorted(
            sequences_and_read_names.items(),
            key=lambda x: -len(x[1]))):
        if i >= max_sequences:
            break

        count = len(read_names)

        if count < min_reads_supporting_rna_sequence:
            break

        for protein_fragment in translate_compatible_reading_frames(
                dna_sequence_prefix=prefix,
                dna_sequence_variant=variant_seq,
                dna_sequence_suffix=suffix,
                base1_variant_start_location=variant.start,
                base1_variant_end_location=variant.end,
                variant_is_insertion=len(variant.ref) == 0,
                transcripts=reference_transcripts,
                max_transcript_mismatches=max_transcript_mismatches,
                min_transcript_prefix_length=max_transcript_mismatches):
            result_list.append((protein_fragment, count))
    return group_protein_fragments(result_list)

def translate_variant_collection(
        variants,
        samfile,
        transcript_id_whitelist=None,
        protein_fragment_length=PROTEIN_FRAGMENT_LEGNTH,
        min_reads_supporting_rna_sequence=MIN_READS_SUPPORTING_RNA_SEQUENCE,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
        max_transcript_mismatches=MAX_TRANSCRIPT_MISMATCHES,
        max_sequences_per_variant=MAX_SEQUENCES_PER_VARIANT):
    """
    Returns a dictionary mapping each variant to a list of protein fragment,
    read count pairs.
    """
    variant_to_protein_fragments = {}

    chromosome_names = set(samfile.references)
    for variant in variants:
        chromosome = variant.contig

        # I imagine the conversation went like this:
        # A: "Hey, I have an awesome idea"
        # B: "What's up?"
        # A: "Let's make two nearly identical reference genomes"
        # B: "But...that sounds like it might confuse people."
        # A: "Nah, it's cool, we'll give the chromosomes different prefixes!"
        # B: "OK, sounds like a good idea."
        if chromosome not in chromosome_names:
            if "chr" + chromosome in chromosome_names:
                chromosome = "chr" + chromosome
            else:
                logging.warn(
                    "Chromosome '%s' from variant %s not in alignment file %s" % (
                        chromosome, variant, samfile))
                continue

        start = variant.start
        end = variant.end
        # because of the "1" vs "chr1" mess the contig we use to search
        # Ensembl may be different than the one we use in the BAM file
        ensembl_contig = variant.contig
        reference_transcripts = [
            transcript
            for transcript in variant.transcripts
            if transcript.contains(ensembl_contig, start, end) and transcript.complete
        ]

        if transcript_id_whitelist:
            n_transcripts_before_filtering = len(reference_transcripts)
            reference_transcripts = [
                transcript for transcript in reference_transcripts
                if transcript.id in transcript_id_whitelist
            ]
            n_transcripts_after_filtering = len(reference_transcripts)
            n_dropped = (
                n_transcripts_before_filtering - n_transcripts_after_filtering)
            if n_dropped > 0:
                logging.info("Dropped %d/%d candidate transcripts for %s" % (
                    n_dropped,
                    n_transcripts_before_filtering,
                    variant))

        if len(reference_transcripts) == 0:
            continue

        variant_to_protein_fragments[variant] = \
            variant_protein_fragments_with_read_counts(
                chromosome_name=chromosome,
                variant=variant,
                samfile=samfile,
                reference_transcripts=reference_transcripts,
                protein_fragment_length=protein_fragment_length,
                min_reads_supporting_rna_sequence=min_reads_supporting_rna_sequence,
                min_transcript_prefix_length=min_transcript_prefix_length,
                max_transcript_mismatches=max_transcript_mismatches)
    return variant_to_protein_fragments


def variant_protein_fragments_dataframe(
        variants,
        samfile,
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

        base1_start_pos : int
            First reference nucleotide affected by variant (or position before
            insertion)

        base1_end_pos : int
            Last reference nucleotide affect by variant (or position after
            insertion)

        ref : str
            Reference nucleotides
        alt : str
            Variant nucleotides

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
    variant_to_proteins = translate_variant_collection(
        variants,
        samfile,
        protein_fragment_length=protein_fragment_length,
        min_reads_supporting_rna_sequence=min_reads_supporting_rna_sequence,
        min_transcript_prefix_length=min_transcript_prefix_length,
        max_transcript_mismatches=max_transcript_mismatches,
        max_sequences_per_variant=max_sequences_per_variant)

    # construct a dictionary incrementally which we'll turn into a
    # DataFrame
    column_dict = OrderedDict([
        ("chr", []),
        ("base1_start_pos", []),
        ("base1_end_pos", []),
        ("ref", []),
        ("alt", []),
        ("variant_protein_sequence", []),
        ("variant_protein_sequence_length", []),
        ("reference_transcript_ids", []),
        ("reference_transcript_names", []),
        ("reference_protein_sequences", []),
        ("cdna_sequences_to_read_names", []),
        ("total_supporting_read_count", [])
    ])
    for (variant, protein_fragments) in variant_to_proteins.items():
        for protein_fragment in protein_fragments:
            column_dict["chr"].append(variant.contig)
            column_dict["base1_start_pos"].append(variant.start)
            column_dict["base1_end_pos"].append(variant.end)
            column_dict["ref"].append(variant.ref)
            column_dict["alt"].append(variant.alt)
            column_dict["variant_protein_sequence"].append(
                protein_fragment.variant_protein_sequence)
            column_dict["variant_protein_sequence_length"].append(
                len(protein_fragment.variant_protein_sequence))
            column_dict["reference_transcript_ids"].append(
                ";".join(protein_fragment.reference_transcript_ids))
            column_dict["reference_transcript_names"].append(
                ";".join(protein_fragment.reference_transcript_names))
            column_dict["reference_protein_sequences"].append(
                ";".join(protein_fragment.reference_protein_sequences))
            column_dict["cdna_sequences_to_read_names"].append(
                protein_fragment.cdna_sequences_to_read_names)
            column_dict["total_supporting_read_count"].append(
                protein_fragment.number_supporting_reads)
    df = pd.DataFrame(column_dict)
    return df
