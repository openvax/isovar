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
from argparse import ArgumentTypeError

from varcode.cli import make_variants_parser

from ..default_parameters import (
    MIN_READ_MAPPING_QUALITY,
    USE_DUPLICATE_READS,
    USE_SECONDARY_ALIGNMENTS,
    USE_SOFT_CLIPPED_BASES,
    USE_READS_WITHOUT_BASE_QUALITIES,
    MERGE_OVERLAPPING_FRAGMENTS,
    NUM_RNA_DECOMPRESSION_THREADS,
    INFER_READ_ENDS, TRIM_ADAPTERS, TRIM_POLY_A, READ_END_PROFILE,
)

from ..read_collector import ReadCollector
from ..read_end_inference import read_end_profiles_from_json
from ..dataframe_helpers import allele_counts_dataframe, allele_reads_to_dataframe
from .validation import CommandInputError, non_negative_int, positive_int, variant_collection_from_args


def read_end_profile_argument(filename):
    """Report explicit profile errors as argument errors, not tracebacks."""
    try:
        return read_end_profiles_from_json(filename)
    except (OSError, ValueError, TypeError, AttributeError) as error:
        raise ArgumentTypeError(str(error)) from error


def add_rna_args(parser):
    """
    Extends an ArgumentParser instance with the following commandline arguments:
        --bam
        --min-mapping-quality
        --use-duplicate-reads
        --drop-secondary-alignments
        --use-soft-clipped-bases
        --require-base-qualities
        --infer-read-ends
        --read-end-profile
        --trim-adapters
        --trim-poly-a
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
        type=non_negative_int,
        default=MIN_READ_MAPPING_QUALITY,
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
        "--require-base-qualities",
        dest="use_reads_without_base_qualities",
        default=USE_READS_WITHOUT_BASE_QUALITIES,
        action="store_false",
        help="Discard reads with missing QUAL instead of retaining sequence/alignment evidence with unknown base quality.")

    rna_group.add_argument(
        "--infer-read-ends", action="store_true", default=INFER_READ_ENDS,
        help="Annotate candidate adapter/poly-A/T ends without changing read sequence.")
    rna_group.add_argument(
        "--read-end-profile", type=read_end_profile_argument, default=READ_END_PROFILE,
        help="JSON adapter/kit profile, optionally keyed by read_groups; also enables end annotation.")
    rna_group.add_argument(
        "--trim-adapters", action="store_true", default=TRIM_ADAPTERS,
        help="Trim inferred adapters from terminal soft clips (requires --read-end-profile; raw BAM unchanged).")
    rna_group.add_argument(
        "--trim-poly-a", action="store_true", default=TRIM_POLY_A,
        help="Trim candidate poly-A/T tails from terminal soft clips; enables inference, preserves aligned bases.")

    rna_group.add_argument(
        "--no-merge-overlapping-fragments",
        dest="merge_overlapping_fragments",
        action="store_false",
        default=MERGE_OVERLAPPING_FRAGMENTS,
        help="Keep overlapping mates as separate reads (merge by default: %(default)s).")

    rna_group.add_argument(
        "--num-rna-decompression-threads",
        type=positive_int,
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
    Use parsed arguments to open an indexed BAM/CRAM file of aligned RNA reads.

    Raises CommandInputError if it is missing, unreadable, SAM or unindexed.
    """
    try:
        alignment_file = AlignmentFile(
            args.bam,
            threads=args.num_rna_decompression_threads)
    except (OSError, ValueError) as error:
        raise CommandInputError("Cannot open alignment file %s: %s" % (args.bam, error)) from error
    if not (alignment_file.is_bam or alignment_file.is_cram):
        alignment_file.close()
        raise CommandInputError(
            "%s is not BAM or CRAM; region queries need an indexed BAM/CRAM file" % args.bam)
    if not alignment_file.has_index():
        alignment_file.close()
        raise CommandInputError("%s has no index; create one with 'samtools index'" % args.bam)
    return alignment_file


def read_collector_from_args(args):
    """
    Use parsed arguments to create a ReadCollector object
    """
    try:
        return _read_collector_from_args(args)
    except (TypeError, ValueError) as error:
        raise CommandInputError(str(error)) from error


def _read_collector_from_args(args):
    return ReadCollector(
        min_mapping_quality=args.min_mapping_quality,
        use_duplicate_reads=args.use_duplicate_reads,
        use_secondary_alignments=not args.drop_secondary_alignments,
        use_soft_clipped_bases=args.use_soft_clipped_bases,
        infer_read_ends=getattr(args, "infer_read_ends", INFER_READ_ENDS),
        read_end_profile=getattr(args, "read_end_profile", READ_END_PROFILE),
        trim_adapters=getattr(args, "trim_adapters", TRIM_ADAPTERS),
        trim_poly_a=getattr(args, "trim_poly_a", TRIM_POLY_A),
        use_reads_without_base_qualities=getattr(
            args, "use_reads_without_base_qualities", USE_READS_WITHOUT_BASE_QUALITIES),
        merge_overlapping_fragments=getattr(
            args, "merge_overlapping_fragments", MERGE_OVERLAPPING_FRAGMENTS))


def read_evidence_generator_from_args(args):
    """
    Creates a generator of (Variant, ReadEvidence) pairs from parsed
    arguments.
    """
    read_creator = read_collector_from_args(args)
    variants = variant_collection_from_args(args)
    samfile = alignment_file_from_args(args)
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


def allele_counts_dataframe_from_args(args):
    """
    Collect read and fragment counts for each variant and turn them into a
    DataFrame, with cell/UMI label counts when ``--cell-umi-labels`` is given.
    """
    pairs = list(read_evidence_generator_from_args(args))
    df = allele_counts_dataframe(pairs)
    if getattr(args, "cell_umi_labels", False):
        from ..cell_evidence import ALLELES, CellUmiAlleles
        with alignment_file_from_args(args) as alignment_file:
            labels = CellUmiAlleles(alignment_file, sample_id=args.sample_id, source=args.source or args.bam,
                                    read_collector=read_collector_from_args(args))
            evidence = [labels.evidence(variant, read_evidence) for variant, read_evidence in pairs]
        for allele in ALLELES:
            for column, field in (("cells", "observed_cells"), ("cell_umi_labels", "observed_labels"),
                                  ("unlabeled_segments", "unresolved_segments")):
                df["num_%s_%s" % (allele, column)] = [e["alleles"][allele][field] for e in evidence]
        df["num_cells_with_ref_and_alt"] = [e["cells_with_ref_and_alt"] for e in evidence]
    return df


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

