# Copyright (c) 2016-2018. Mount Sinai School of Medicine
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

from pysam import AlignmentFile

from varcode.cli import make_variants_parser, variant_collection_from_args

from ..default_parameters import MIN_READ_MAPPING_QUALITY
from ..allele_reads import reads_overlapping_variants, reads_to_dataframe
from ..variant_reads import reads_supporting_variants


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
        action="store_true")

    rna_group.add_argument(
        "--drop-secondary-alignments",
        default=False,
        action="store_true")
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


def allele_reads_generator_from_args(args):
    variants = variant_collection_from_args(args)
    samfile = samfile_from_args(args)
    return reads_overlapping_variants(
        variants=variants,
        samfile=samfile,
        use_duplicate_reads=args.use_duplicate_reads,
        use_secondary_alignments=not args.drop_secondary_alignments,
        min_mapping_quality=args.min_mapping_quality)


def allele_reads_dataframe_from_args(args):
    return reads_to_dataframe(allele_reads_generator_from_args(args))


def variant_reads_generator_from_args(args):
    variants = variant_collection_from_args(args)
    samfile = samfile_from_args(args)
    return reads_supporting_variants(
        variants=variants,
        samfile=samfile,
        use_duplicate_reads=args.use_duplicate_reads,
        use_secondary_alignments=not args.drop_secondary_alignments,
        min_mapping_quality=args.min_mapping_quality)


def variant_reads_dataframe_from_args(args):
    return reads_to_dataframe(variant_reads_generator_from_args(args))
