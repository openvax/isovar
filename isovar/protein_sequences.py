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
from collections import namedtuple

from skbio import DNA
import numpy as np

from .variant_reads import gather_variant_reads
from .sequence_counts import sequence_counts

ProteinFragment = namedtuple(
    "ProteinFragment",
    [
        "cdna_prefix",
        "cdna_variant",
        "cdna_suffix",
        "transcript_id",
        "number_transcript_sequence_mismatches",
        "fraction_transcript_sequence_mismatches",
        "reading_frame_at_start_of_cdna_sequence",
        "transcript_sequence_before_variant",
        "variant_protein_sequence",
        "transcript_protein_sequence",
        "protein_sequence_offset",
    ])

MIN_READS_SUPPORTING_RNA_SEQUENCE = 3
MIN_TRANSCRIPT_PREFIX_LENGTH = 15
MAX_TRANSCRIPT_MISMATCHES = 2
PROTEIN_FRAGMENT_LEGNTH = 25
MAX_SEQUENCES_PER_VARIANT = 5


def matching_isoform_fragment_protein_sequences(
        dna_sequence_prefix,
        dna_sequence_variant,
        dna_sequence_suffix,
        variant_is_insertion,
        base1_variant_start_location,
        base1_variant_end_location,
        transcripts,
        max_transcript_mismatches=MAX_TRANSCRIPT_MISMATCHES,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH):
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
            subseq_end_idx = variant_in_transcript_idx + 1
        else:
            subseq_end_idx = variant_in_transcript_idx

        start_codon_idx = min(transcript.start_codon_spliced_offsets)
        if variant_in_transcript_idx < start_codon_idx + 3:
            logging.info(
                "Skipping %s because variant is not after start codon" % (
                    transcript,))

        subseq_start_idx = subseq_end_idx - len(cdna_prefix)
        if subseq_start_idx < 0:
            continue

        transcript_sequence_before_variant = transcript.sequence[
            subseq_start_idx:subseq_end_idx]

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

        reading_frame = (subseq_start_idx - start_codon_idx) % 3

        combined_variant_cdna_sequence = DNA(
            cdna_prefix + cdna_variant + cdna_suffix)
        in_frame_combined_variant_sequence = combined_variant_cdna_sequence[
            reading_frame:]
        variant_protein_fragment = in_frame_combined_variant_sequence.translate()

        combined_transcript_sequence = DNA(
            transcript.sequence[
                subseq_start_idx:
                subseq_start_idx + len(combined_variant_cdna_sequence)])
        in_frame_transcript_sequence = combined_transcript_sequence[
            reading_frame:]
        transcript_protein_sequence = in_frame_transcript_sequence.translate()
        results.append(
            ProteinFragment(
                cdna_prefix=cdna_prefix,
                cdna_variant=cdna_variant,
                cdna_suffix=cdna_suffix,
                transcript_id=transcript.id,
                number_transcript_sequence_mismatches=n_mismatch,
                fraction_transcript_sequence_mismatches=fraction_mismatch,
                reading_frame_at_start_of_cdna_sequence=reading_frame,
                transcript_sequence_before_variant=transcript_sequence_before_variant,
                variant_protein_sequence=variant_protein_fragment,
                transcript_protein_sequence=transcript_protein_sequence,
                protein_sequence_offset=subseq_start_idx - start_codon_idx % 3))
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

    rna_sequence_length = protein_fragment_length * 3

    # the number of context nucleotides on either side of the variant
    # is half the desired length (minus the number of variant nucleotides)
    n_surrounding_nucleotides = rna_sequence_length - len(variant.alt)
    flanking_context_size = int(np.ceil(n_surrounding_nucleotides / 2.0))

    sequence_count_info = sequence_counts(
        variant_reads,
        context_size=flanking_context_size)

    sequence_count_dict = sequence_count_info.full_read_counts

    if len(sequence_count_dict) == 0:
        return []

    variant_seq = sequence_count_info.variant_nucleotides

    result_list = []

    for i, ((prefix, suffix), count) in enumerate(sorted(
            sequence_count_dict.items(),
            key=lambda x: -x[1])):
        if i >= max_sequences:
            break

        if count < min_reads_supporting_rna_sequence:
            break

        for protein_fragment in matching_isoform_fragment_protein_sequences(
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
    return result_list

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
