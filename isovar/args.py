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
)

def extend_parser_with_somatic_vcf_args(
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

def extend_parser_with_rna_args(
        parser,
        rna_bam_arg="--bam",
        rna_bam_required=True,
        min_rna_reads_arg="--min-reads",
        min_rna_reads_default=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE):
    """
    Extends an ArgumentParser instance with the following commandline arguments:
        --bam
        --min-reads
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

    return parser

def extend_parser_with_reference_context_args(
        parser,
        max_reference_transcript_mismatches_arg="--max-reference-transcript-mismatches",
        max_reference_transcript_mismatches_default=MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
        min_transcript_prefix_length_arg="--min-transcript-prefix-length",
        min_transcript_prefix_length_default=MIN_TRANSCRIPT_PREFIX_LENGTH):
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
            "near the start codon of transcripts without 5' UTRs."))
    return parser
