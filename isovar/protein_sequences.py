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

MIN_READS_SUPPORTING_CONTEXT = 3
MIN_TRANSCRIPT_PREFIX_LENGTH = 15
MAX_TRANSCRIPT_MISMATCHES = 2


def matching_isoform_fragment_protein_sequences(
        variant,
        dna_sequence_prefix,
        dna_sequence_variant,
        dna_sequence_suffix,
        transcripts,
        max_transcript_mismatches=MAX_TRANSCRIPT_MISMATCHES,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH):
    results = []

    for transcript in transcripts:
        if not transcript.contains(variant.contig, variant.start, variant.end):
            logging.info(
                "Skipping transcript %s because it does not overlap %s" % (
                    transcript, variant))
            continue

        if not transcript.complete:
            continue

        if transcript.strand == "+":
            cdna_prefix = dna_sequence_prefix
            cdna_suffix = dna_sequence_suffix
            cdna_variant = dna_sequence_variant
            variant_in_transcript_idx = transcript.spliced_offset(variant.start)
        else:
            # if the transcript is on the reverse strand then we have to
            # take the sequence PREFIX|VARIANT|SUFFIX
            # and take the complement of XIFFUS|TNAIRAV|XIFERP
            cdna_prefix = str(DNA(dna_sequence_suffix).reverse_complement())
            cdna_suffix = str(DNA(dna_sequence_prefix).reverse_complement())
            cdna_variant = str(
                DNA(dna_sequence_variant).reverse_complement())
            variant_in_transcript_idx = transcript.spliced_offset(variant.end)

        if len(variant.ref) == 0:
            # insertions don't actually affect the base referred to
            # by the start position of the variant, but rather the
            # variant gets inserted *after* that position
            subseq_end_idx = variant_in_transcript_idx + 1
        else:
            subseq_end_idx = variant_in_transcript_idx

        start_codon_idx = min(transcript.start_codon_spliced_offsets)
        if variant_in_transcript_idx < start_codon_idx + 3:
            logging.info(
                "Skipping %s because variant %s is not after start codon" % (
                    transcript, variant))

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

def variant_to_protein_sequences(
        variant,
        samfile,
        transcript_id_whitelist=None,
        sequence_context_size=45,
        min_reads_supporting_context=MIN_READS_SUPPORTING_CONTEXT,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
        max_transcript_mismatches=MAX_TRANSCRIPT_MISMATCHES):
    variant_reads = gather_variant_reads(
        samfile=samfile,
        chromosome="chr" + variant.contig,
        base1_location=variant.start,
        ref=variant.ref,
        alt=variant.alt)
    if len(variant_reads) == 0:
        return {}
    collapsed_sequences = sequence_counts(
        variant_reads, context_size=sequence_context_size)
    sequence_to_count_dict = collapsed_sequences.full_read_counts

    if len(sequence_to_count_dict) == 0:
        return {}

    variant_seq = collapsed_sequences.variant_nucleotides

    candidate_transcripts = variant.transcripts

    if transcript_id_whitelist:
        candidate_transcripts = [
            transcript
            for transcript in candidate_transcripts
            if transcript.id in transcript_id_whitelist
        ]

    if len(candidate_transcripts) == 0:
        return {}

    result_list = []
    for ((prefix, suffix), read_count) in sequence_to_count_dict.items():
        if read_count < min_reads_supporting_context:
            continue
        for protein_fragment in matching_isoform_fragment_protein_sequences(
                variant=variant,
                dna_sequence_prefix=prefix,
                dna_sequence_variant=variant_seq,
                dna_sequence_suffix=suffix,
                transcripts=candidate_transcripts,
                max_transcript_mismatches=max_transcript_mismatches,
                min_transcript_prefix_length=max_transcript_mismatches):
            result_list.append((protein_fragment, read_count))
    return result_list

def variants_to_protein_sequences(
        variants,
        samfile,
        transcript_id_whitelist=None,
        sequence_context_size=45,
        min_reads_supporting_context=MIN_READS_SUPPORTING_CONTEXT,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
        max_transcript_mismatches=MAX_TRANSCRIPT_MISMATCHES):
    """
    Returns a dictionary mapping each variant to a list of protein fragment,
    read count pairs.
    """
    variant_to_protein_fragments = {}
    for variant in variants:
        variant_to_protein_fragments[variant] = \
            variant_to_protein_sequences(
                variant=variant,
                samfile=samfile,
                transcript_id_whitelist=transcript_id_whitelist,
                sequence_context_size=sequence_context_size,
                min_reads_supporting_context=min_reads_supporting_context,
                min_transcript_prefix_length=min_transcript_prefix_length,
                max_transcript_mismatches=max_transcript_mismatches)
    return variant_to_protein_fragments
