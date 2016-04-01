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
Collect reads containing a variant and split them into prefix, variant, and
suffix portions
"""

from __future__ import print_function, division, absolute_import

from collections import namedtuple, OrderedDict

from pandas import DataFrame

from .reads_at_locus import (
    gather_reads_at_locus,
    DEFAULT_MIN_MAPPING_QUALITY,
    DEFAULT_USE_SECONDARY_ALIGNMENTS,
    DEFAULT_USE_DUPLICATE_READS
)
from .logging import create_logger
from .variant_helpers import trim_variant

logger = create_logger(__name__)

VariantRead = namedtuple(
    "VariantRead", "prefix variant suffix name")

def variant_reads_from_overlapping_reads(reads, ref, alt):
    """
    Given a collection of pysam.AlignedSegment objects, generates a
    sequence of VariantRead objects (which are split into prefix/variant/suffix
    nucleotides).

    Parameters
    ----------
    reads : sequence of ReadAtLocus records

    ref : str
        Reference sequence of the variant (empty for insertions)

    alt : str
        Alternate sequence of the variant (empty for deletions)

    Returns a sequence of VariantRead objects.
    """
    for read in reads:
        reference_positions = read.reference_positions
        offset_before = read.offset_before_variant
        offset_after = read.offset_after_variant

        sequence = read.sequence

        #################
        #
        # INSERTIONS
        #
        #################
        if len(ref) == 0:

            if len(reference_positions) < end_offset + 1:
                # if we can't check the base after the inserted nucleotides
                logger.debug(
                    "Skipping read, can't get reference position after insertion")
                continue
            elif start_offset == 0:
                # if we can't check the base before the variant, also skip the read
                logger.debug(
                    "Skipping read, can't get reference position before insertion")
                continue

            # insertions require a sequence of non-aligned bases
            # followed by the subsequence reference position
            insert_positions = reference_positions[start_offset:end_offset + 1]
            if any(insert_pos is not None for insert_pos in insert_positions):
                # all these inserted nucleotides should *not* align to the
                # reference
                logger.debug(
                    "Skipping read, inserted nucleotides shouldn't map to reference")
                continue
            ref_pos_before_insertion = reference_positions[start_offset - 1]
            ref_pos_after_insertion = reference_positions[end_offset + 1]
            if ref_pos_after_insertion != ref_pos_before_insertion + 1:
                logger.debug(
                    "Skipping read, positions before+after insertion should be adjacent")
                continue
            prefix = sequence[:start_offset]
            suffix = sequence[end_offset + 1:]
        ###################
        #
        # DELETIONS
        #
        ###################
        elif len(alt) == 0:
            # the meaning of the locus is different for deletions, since
            # these offsets point to the bases before and after the deleted
            # nucleotides

            ref_pos_before = reference_positions[start_offset]
            ref_pos_after = reference_positions[end_offset]

            if ref_pos_after - ref_pos_before != len(ref):
                # if the number of nucleotides skipped isn't the same
                # as the number deleted in the variant then
                # don't use this read
                logger.debug("Positions before and after deletion should be adjacent")
                continue
            # capture the prefix sequence including the base before the deletion
            prefix = sequence[:start_offset + 1]
            # capture the suffix sequence including the base after the deletion
            suffix = sequence[end_offset:]
        ###################
        #
        # SUBSTITUTIONS
        #
        ###################
        else:
            # deletions and substitutions work similarly, we just need
            # all the reference bases to be adjacently aligned
            ref_pos_start = reference_positions[start_offset]
            ref_pos_end = reference_positions[end_offset]
            if ref_pos_end - ref_pos_start + 1 != len(ref):
                logger.debug(
                    "Positions before and after substitution should be adjacent: %d:%d" % (
                        ref_pos_start, ref_pos_end))
                continue
            prefix = sequence[:start_offset]
            suffix = sequence[end_offset + 1:]
        if isinstance(prefix, bytes):
            prefix = str(prefix, "ascii")
        if isinstance(suffix, bytes):
            suffix = str(suffix, "ascii")
        yield VariantRead(prefix, alt, suffix, name=read.name)


def gather_reads_for_single_variant(
        samfile,
        chromosome,
        variant,
        use_duplicate_reads=DEFAULT_USE_DUPLICATE_READS,
        use_secondary_alignments=DEFAULT_USE_SECONDARY_ALIGNMENTS,
        min_mapping_quality=DEFAULT_MIN_MAPPING_QUALITY):
    """
    Find reads in the given SAM/BAM file which overlap the given variant, filter
    to only include those which agree with the variant's nucleotide(s), and turn
    them into a list of VariantRead objects.

    Parameters
    ----------
    samfile : pysam.AlignmentFile

    chromosome : str

    variant : varcode.Variant

    use_duplicate_reads : bool
        Should we use reads that have been marked as PCR duplicates

    use_secondary_alignments : bool
        Should we use reads at locations other than their best alignment

    min_mapping_quality : int
        Drop reads below this mapping quality

    Returns list of VariantRead objects.
    """
    logger.info("Gathering variant reads for variant %s (chromosome=%s)" % (
        variant,
        chromosome))

    base1_position, ref, alt = trim_variant(variant)
    if len(ref) == 0:
        # if the variant is an insertion
        base1_position_before_variant = base1_position
        base1_position_after_variant = base1_position + 1
    else:
        base1_position_before_variant = base1_position - 1
        base1_position_after_variant = base1_position + len(ref)

    reads_at_locus_generator = gather_reads_at_locus(
        samfile=samfile,
        chromosome=chromosome,
        base1_position_before_variant=base1_position_before_variant,
        base1_position_after_variant=base1_position_after_variant,
        use_duplicate_reads=use_duplicate_reads,
        use_secondary_alignments=use_secondary_alignments,
        min_mapping_quality=min_mapping_quality)
    variant_reads = list(
        variant_reads_from_overlapping_reads(
            overlapping_reads=overlapping_reads,
            ref=ref,
            alt=alt))
    logger.info("Variant reads: %s" % (variant_reads,))
    return variant_reads

def variant_reads_generator(
        variants,
        samfile,
        use_duplicate_reads=False,
        use_secondary_alignments=True,
        min_mapping_quality=5):
    """
    Generates sequence of tuples, each containing a variant paired with
    a list of VariantRead objects.

    Parameters
    ----------
    variants : varcode.VariantCollection

    samfile : pysam.AlignmentFile

    use_duplicate_reads : bool
        Should we use reads that have been marked as PCR duplicates

    use_secondary_alignments : bool
        Should we use reads at locations other than their best alignment

    min_mapping_quality : int
        Drop reads below this mapping quality
    """
    chromosome_names = set(samfile.references)
    for variant in variants:
        # I imagine the conversation went like this:
        # A: "Hey, I have an awesome idea"
        # B: "What's up?"
        # A: "Let's make two nearly identical reference genomes"
        # B: "But...that sounds like it might confuse people."
        # A: "Nah, it's cool, we'll give the chromosomes different prefixes!"
        # B: "OK, sounds like a good idea."
        if variant.contig in chromosome_names:
            chromosome = variant.contig
        elif "chr" + variant.contig in chromosome_names:
            chromosome = "chr" + variant.contig
        else:
            logger.warn(
                "Chromosome '%s' from variant %s not in alignment file %s" % (
                    chromosome, variant, samfile.filename))
            continue

        variant_reads = gather_reads_for_single_variant(
            samfile=samfile,
            chromosome=chromosome,
            variant=variant,
            use_duplicate_reads=use_duplicate_reads,
            use_secondary_alignments=use_secondary_alignments,
            min_mapping_quality=min_mapping_quality)
        logger.info("%s => %s" % (
            variant,
            variant_reads))
        yield variant, variant_reads


def variant_reads_dataframe(variants, samfile):
    columns = OrderedDict([
        ("chr", []),
        ("pos", []),
        ("ref", []),
        ("alt", []),
        ("read_name", []),
        ("read_prefix", []),
        ("read_variant", []),
        ("read_suffix", []),
    ])
    for variant, variant_reads in variant_reads_generator(variants, samfile):
        for variant_read in variant_reads:
            columns["chr"].append(variant.contig)
            columns["pos"].append(variant.start)
            columns["ref"].append(variant.ref)
            columns["alt"].append(variant.alt)
            columns["read_name"].append(variant_read.name)
            columns["read_prefix"].append(variant_read.prefix)
            columns["read_variant"].append(variant_read.variant)
            columns["read_suffix"].append(variant_read.suffix)
    return DataFrame(columns)
