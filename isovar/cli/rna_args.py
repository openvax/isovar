# Copyright (c) 2016-2019. Mount Sinai School of Medicine
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
Common command-line arguments for all Isovar commands which use RNA
"""

from __future__ import print_function, division, absolute_import

from pysam import AlignmentFile

from varcode.cli import make_variants_parser, variant_collection_from_args

from ..default_parameters import (
    MIN_READ_MAPPING_QUALITY,
    MIN_ALT_RNA_READS,
    MIN_ALT_RNA_FRAGMENTS,
    MIN_RNA_VAF,
    MIN_RATIO_ALT_TO_OTHER_NONREF_RNA_FRAGMENTS
)
from ..read_creator import  ReadCreator
from ..dataframe_helpers import allele_reads_to_dataframe


def add_rna_args(
        parser,
        min_mapping_quality_default=MIN_READ_MAPPING_QUALITY):
    """
    Extends an ArgumentParser instance with the following commandline arguments:
        --bam
        --min-reads
        --min-mapping-quality
        --use-duplicate-reads
        --drop-secondary-alignments
    """
    rna_group = parser.add_argument_group("RNA")
    rna_group.add_argument(
        "--bam",
        required=True,
        help="BAM file containing RNAseq reads")

    rna_group.add_argument(
        "--min-mapping-quality",
        type=int,
        default=min_mapping_quality_default,
        help="Minimum MAPQ value to allow for a read (default %(default)s)")

    rna_group.add_argument(
        "--use-duplicate-reads",
        default=False,
        action="store_true",
        help=(
            "By default, reads which have been marked as duplicates are excluded."
            "Use this option to include duplicate reads."))

    rna_group.add_argument(
        "--drop-secondary-alignments",
        default=False,
        action="store_true",
        help=(
            "By default, secondary alignments are included in reads, "
            "use this option to instead only use primary alignments."))

    rna_group.add_argument(
        "--min-alt-rna-reads",
        type=int,
        default=MIN_ALT_RNA_READS,
        help="Minimum number of reads supporting variant allele (default %(default)s)")

    rna_group.add_argument(
        "--min-alt-rna-fragments",
        type=int,
        default=MIN_ALT_RNA_FRAGMENTS,
        help=(
            "Minimum number of fragments supporting variant allele (default %(default)s). "
            "Note that this option is the same as --min-alt-rna-reads for single-end "
            "sequencing."))

    rna_group.add_argument(
        "--min-rna_vaf",
        type=float,
        default=MIN_RNA_VAF,
        help=(
            "Minimum ratio of fragments supporting variant allele to total RNA fragments "
            "(default %(default)s)."))

    rna_group.add_argument(
        "--min-ratio-alt-to-other-nonref-fragments",
        type=float,
        default=MIN_RATIO_ALT_TO_OTHER_NONREF_RNA_FRAGMENTS,
        help=(
            "At loci where alleles other than the ref and a single alt are supported, "
            "this parameter controls how many more times fragments supporting "
            "the variant allele are required relative to other non-reference "
            "alleles (default %(default)s)."))

    return rna_group


def make_rna_reads_arg_parser(**kwargs):
    """
    Parameters
    ----------
    **kwargs : dict
        Passed directly to argparse.ArgumentParser

    Creates argparse.ArgumentParser instance with all of the options
    needed to load a set of variants and their supporting RNA reads.
    """
    parser = make_variants_parser(**kwargs)
    add_rna_args(parser)
    return parser


def samfile_from_args(args):
    return AlignmentFile(args.bam)

def read_creator_from_args(args):
    return ReadCreator(
        min_mapping_quality=args.min_mapping_quality,
        use_duplicate_reads=args.use_duplicate_reads,
        use_secondary_alignments=not args.drop_secondary_alignments)

def overlapping_reads_generator_from_args(args):
    variants = variant_collection_from_args(args)
    samfile = samfile_from_args(args)
    read_creator = read_creator_from_args(args)
    return read_creator.allele_reads_overlapping_variants(
        variants=variants,
        alignments=samfile)

def allele_reads_dataframe_from_args(args):
    return allele_reads_to_dataframe(overlapping_reads_generator_from_args(args))

def supporting_reads_generator_from_args(args):
    variants = variant_collection_from_args(args)
    samfile = samfile_from_args(args)
    read_creator = read_creator_from_args(args)
    return read_creator.allele_reads_supporting_variants(
        variants=variants,
        alignments=samfile)


def variant_reads_dataframe_from_args(args):
    return allele_reads_to_dataframe(supporting_reads_generator_from_args(args))
