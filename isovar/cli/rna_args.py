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

from pysam import AlignmentFile

from varcode.cli import make_variants_parser, variant_collection_from_args

from ..default_parameters import (
    MIN_READ_MAPPING_QUALITY,
    USE_DUPLICATE_READS,
    USE_SECONDARY_ALIGNMENTS,
    USE_SOFT_CLIPPED_BASES,
    MERGE_OVERLAPPING_FRAGMENTS,
    NUM_RNA_DECOMPRESSION_THREADS,
)

from ..read_collector import ReadCollector
from ..dataframe_helpers import (
    allele_counts_dataframe,
    allele_reads_to_dataframe,
    read_evidence_generator_to_dataframe,
)


def add_rna_args(
        parser,
        min_mapping_quality_default=MIN_READ_MAPPING_QUALITY):
    """
    Extends an ArgumentParser instance with the following commandline arguments:
        --bam
        --min-mapping-quality
        --use-duplicate-reads
        --drop-secondary-alignments
        --use-soft-clipped-bases
        --no-merge-overlapping-fragments
        --num-rna-decompression-threads
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
        default=USE_DUPLICATE_READS,
        action="store_true",
        help=(
            "Include reads marked as duplicates (default %(default)s)."))

    rna_group.add_argument(
        "--drop-secondary-alignments",
        default=not USE_SECONDARY_ALIGNMENTS,
        action="store_true",
        help=(
            "Drop secondary alignments (default %(default)s)."))

    rna_group.add_argument(
        "--use-soft-clipped-bases",
        default=USE_SOFT_CLIPPED_BASES,
        action="store_true",
        help=(
            "Include soft-clipped bases at the ends of reads (default %(default)s)."))

    rna_group.add_argument(
        "--no-merge-overlapping-fragments",
        dest="merge_overlapping_fragments",
        action="store_false",
        default=MERGE_OVERLAPPING_FRAGMENTS,
        help="Keep overlapping mates as separate reads (merge by default: %(default)s).")

    rna_group.add_argument(
        "--num-rna-decompression-threads",
        type=int,
        help=(
            "Number of threads to use for decompression of BAM/CRAM files "
            "(default %(default)s)."),
        default=NUM_RNA_DECOMPRESSION_THREADS)

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


def alignment_file_from_args(args):
    """
    Use parsed arguments to load a file of aligned RNA reads.
    """
    return AlignmentFile(
        args.bam,
        threads=args.num_rna_decompression_threads)


def read_collector_from_args(args):
    """
    Use parsed arguments to create a ReadCollector object
    """
    return ReadCollector(
        min_mapping_quality=args.min_mapping_quality,
        use_duplicate_reads=args.use_duplicate_reads,
        use_secondary_alignments=not args.drop_secondary_alignments,
        use_soft_clipped_bases=args.use_soft_clipped_bases,
        merge_overlapping_fragments=getattr(
            args, "merge_overlapping_fragments", MERGE_OVERLAPPING_FRAGMENTS))


def read_evidence_generator_from_args(args):
    """
    Creates a generator of (Variant, ReadEvidence) pairs from parsed
    arguments.
    """
    variants = variant_collection_from_args(args)
    samfile = alignment_file_from_args(args)
    read_creator = read_collector_from_args(args)
    return read_creator.read_evidence_generator(
        variants=variants,
        alignment_file=samfile)

def variant_reads_generator_from_args(args):
    """
    Creates a generator of (Variant, list of AlleleRead) from parsed
    arguments, where all AlleleRead objects must have alleles matching
    the variant.
    """
    for variant, read_evidence in read_evidence_generator_from_args(args):
        yield variant, read_evidence.alt_reads


def allele_reads_generator_from_args(args):
    """
    Creates a generator of (Variant, list of AlleleRead) from parsed
    arguments, including all reads at each variant locus.
    """
    for variant, read_evidence in read_evidence_generator_from_args(args):
        yield variant, (
            list(read_evidence.ref_reads)
            + list(read_evidence.alt_reads)
            + list(read_evidence.other_reads)
        )


def read_evidence_dataframe_from_args(args):
    """
    Collect ReadEvidence for each variant and turn them into a DataFrame
    """
    return read_evidence_generator_to_dataframe(
        read_evidence_generator_from_args(args))


def allele_counts_dataframe_from_args(args):
    """
    Collect read and fragment counts for each variant and turn them into a
    DataFrame.
    """
    return allele_counts_dataframe(
        read_evidence_generator_from_args(args))


def allele_reads_dataframe_from_args(args):
    """
    Collect all allele reads for each variant and turn them into a DataFrame.
    """
    return allele_reads_to_dataframe(
        allele_reads_generator_from_args(args))


def variant_reads_dataframe_from_args(args):
    """
    Collect alt-supporting reads for each variant and turn them into a
    DataFrame.
    """
    return allele_reads_to_dataframe(
        variant_reads_generator_from_args(args))


def variants_reads_dataframe_from_args(args):
    """
    Collect variant reads for each variant and turn them into a DataFrame
    """
    return variant_reads_dataframe_from_args(args)
