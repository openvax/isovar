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

"""
Helper functions for building ArgumentParser instances. Keeps us from
having to copy the same arguments across many scripts.
"""

from __future__ import print_function, division, absolute_import

from .default_parameters import (
    MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
    MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
    MIN_TRANSCRIPT_PREFIX_LENGTH,
    MIN_READ_MAPPING_QUALITY,
    PROTEIN_SEQUENCE_LENGTH,
    MAX_PROTEIN_SEQUENCES_PER_VARIANT,
)

def add_somatic_vcf_args(
        parser,
        vcf_arg="--vcf",
        vcf_required=True,
        genome_arg="--genome",
        genome_required=False):
    """
    Extends an ArgumentParser instance with the following commandline arguments:
        --vcf
        --genome
    """
    parser.add_argument(
        vcf_arg,
        required=vcf_required,
        help="Path to VCF file containing somatic variants")

    parser.add_argument(
        genome_arg,
        default=None,
        required=genome_required,
        help="Name of reference genome for VCF of somatic variants")
    return parser

def add_rna_args(
        parser,
        rna_bam_arg="--bam",
        rna_bam_required=True,
        min_rna_reads_arg="--min-reads",
        min_rna_reads_default=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
        min_mapping_quality_arg="--min-mapping-quality",
        min_mapping_quality_default=MIN_READ_MAPPING_QUALITY):
    """
    Extends an ArgumentParser instance with the following commandline arguments:
        --bam
        --min-reads
        --min-mapping-quality
    """
    parser.add_argument(
        rna_bam_arg,
        required=rna_bam_required,
        help="BAM file containing RNAseq reads")

    parser.add_argument(
        min_rna_reads_arg,
        type=int,
        default=min_rna_reads_default,
        help="Minimum number of reads supporting a variant sequence")

    parser.add_argument(
        "--min-mapping-quality",
        type=int,
        default=min_mapping_quality_default,
        help="Minimum MAPQ value to allow for a read")
    return parser

def add_reference_context_args(
        parser,
        max_reference_transcript_mismatches_arg="--max-reference-transcript-mismatches",
        max_reference_transcript_mismatches_default=MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
        min_transcript_prefix_length_arg="--min-transcript-prefix-length",
        min_transcript_prefix_length_default=MIN_TRANSCRIPT_PREFIX_LENGTH):
    """
    Extends an ArgumentParser instance with the following commandline arguments:
        --max-reference-transcript-mismatches
        --min-transcript-prefix-length
    """
    parser.add_argument(
        max_reference_transcript_mismatches_arg,
        type=int,
        default=max_reference_transcript_mismatches_default,
        help=(
            "Maximum number of mismatches between variant sequence"
            " reference sequence before a candidate reading frame is ignored."))

    parser.add_argument(
        min_transcript_prefix_length_arg,
        type=int,
        default=min_transcript_prefix_length_default,
        help=(
            "Number of nucleotides before the variant we try to match against "
            "a reference transcript. Values greater than zero exclude variants "
            "near the start codon of transcrPROTEIN_SEQUENCE_LEGNTHipts without 5' UTRs."))
    return parser

def add_protein_sequence_args(
        parser,
        protein_sequence_length_arg="--protein-sequence-length",
        protein_sequence_length_default=PROTEIN_SEQUENCE_LENGTH,
        max_protein_sequences_per_variant_arg="--max-protein-sequences-per-variant",
        max_protein_sequences_per_variant_default=MAX_PROTEIN_SEQUENCES_PER_VARIANT):
    """
    Extends an ArgumentParser instance with the following args:
        --protein-sequence-length
        --max-protein-sequences-per-variant
    """
    parser.add_argument(
        protein_sequence_length_arg,
        default=protein_sequence_length_default,
        type=int)
    parser.add_argument(
        max_protein_sequences_per_variant_arg,
        type=int,
        default=max_protein_sequences_per_variant_default)
    return parser
