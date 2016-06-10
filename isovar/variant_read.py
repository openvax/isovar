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
from collections import namedtuple
import logging

from .read_at_locus import read_at_locus_generator
from .default_parameters import (
    MIN_READ_MAPPING_QUALITY,
    USE_SECONDARY_ALIGNMENTS,
    USE_DUPLICATE_READS,
)
from .variant_helpers import trim_variant
from .dataframe_builder import DataFrameBuilder


VariantRead = namedtuple(
    "VariantRead", "prefix alt suffix name")


def trim_N_nucleotides(prefix, suffix):
    """
    Drop all occurrences of 'N' from prefix and suffix nucleotide strings
    by trimming.
    """
    if 'N' in prefix:
        # trim prefix to exclude all occurrences of N
        rightmost_index = prefix.rfind('N')
        logging.debug(
            "Trimming %d nucleotides from read prefix '%s'" % (
                rightmost_index + 1, prefix))
        prefix = prefix[rightmost_index + 1:]

    if 'N' in suffix:
        leftmost_index = suffix.find('N')
        logging.debug(
            "Trimming %d nucleotides from read suffix '%s'" % (
                len(suffix) - leftmost_index,
                suffix))
        suffix = suffix[:leftmost_index]

    return prefix, suffix

def convert_from_bytes_if_necessary(prefix, suffix):
    """
    Depending on how we extract data from pysam we may end up with either
    a string or a byte array of nucleotides. For consistency and simplicity,
    we want to only use strings in the rest of our code.
    """
    if isinstance(prefix, bytes):
        prefix = prefix.decode('ascii')

    if isinstance(suffix, bytes):
        suffix = suffix.decode('ascii')

    return prefix, suffix

def variant_read_from_single_read_at_locus(read, ref, alt):
    """
    Given a single ReadAtLocus object, return either a VariantRead or None
    (if the read's sequence didn't contain the variant nucleotides).

    Parameters
    ----------
    read : ReadAtLocus
        Read which may possibly contain the alternate nucleotides

    ref : str
        Reference sequence of the variant (empty for insertions)

    alt : str
        Alternate sequence of the variant (empty for deletions)

    """
    sequence = read.sequence
    reference_positions = read.reference_positions

    # positions of the nucleotides before and after the variant within
    # the read sequence
    read_pos_before = read.base0_read_position_before_variant
    read_pos_after = read.base0_read_position_after_variant

    # positions of the nucleotides before and after the variant on the
    # reference genome
    ref_pos_before = reference_positions[read_pos_before]

    if ref_pos_before is None:
        logging.warn(
            "Missing reference pos for nucleotide before variant on read: %s" % (
                read,))
        return None

    ref_pos_after = reference_positions[read_pos_after]

    if ref_pos_after is None:
        logging.warn(
            "Missing reference pos for nucleotide after variant on read: %s" % (
                read,))
        return None

    if len(ref) == 0:
        if ref_pos_after - ref_pos_before != 1:
            # if the number of nucleotides skipped isn't the same
            # as the number of reference nucleotides in the variant then
            # don't use this read
            logging.debug(
                "Positions before (%d) and after (%d) variant should be adjacent on read %s" % (
                    ref_pos_before,
                    ref_pos_after,
                    read))
            return None

        # insertions require a sequence of non-aligned bases
        # followed by the subsequence reference position
        ref_positions_for_inserted = reference_positions[
            read_pos_before + 1:read_pos_after]
        if any(insert_pos is not None for insert_pos in ref_positions_for_inserted):
            # all these inserted nucleotides should *not* align to the
            # reference
            logging.debug(
                "Skipping read, inserted nucleotides shouldn't map to reference")
            return None
    else:
        # substitutions and deletions
        if ref_pos_after - ref_pos_before != len(ref) + 1:
            # if the number of nucleotides skipped isn't the same
            # as the number of reference nucleotides in the variant then
            # don't use this read
            logging.debug(
                "Positions before (%d) and after (%d) variant should be adjacent on read %s" % (
                    ref_pos_before,
                    ref_pos_after,
                    read))
            return None

    nucleotides_at_variant_locus = sequence[read_pos_before + 1:read_pos_after]
    if nucleotides_at_variant_locus != alt:
        # read sequence doesn't match variant
        logging.debug("Read '%s' has '%s' at variant locus (%d:%d), required '%s'" % (
            read.name,
            nucleotides_at_variant_locus,
            read_pos_before + 1,
            read_pos_after,
            alt))
        return None

    prefix = sequence[:read_pos_before + 1]
    suffix = sequence[read_pos_after:]

    prefix, suffix = convert_from_bytes_if_necessary(prefix, suffix)
    prefix, suffix = trim_N_nucleotides(prefix, suffix)

    return VariantRead(
        prefix,
        nucleotides_at_variant_locus,
        suffix,
        name=read.name)

def variant_reads_from_reads_at_locus(reads, ref, alt):
    """
    Given a collection of ReadAtLocus objects, returns a
    list of VariantRead objects (which are split into prefix/variant/suffix
    nucleotides).

    Parameters
    ----------
    reads : sequence of ReadAtLocus records

    ref : str
        Reference sequence of the variant (empty for insertions)

    alt : str
        Alternate sequence of the variant (empty for deletions)

    Returns a list of VariantRead objects.
    """
    variant_reads = []
    for read in reads:
        variant_read = variant_read_from_single_read_at_locus(read, ref, alt)
        if variant_read is not None:
            variant_reads.append(variant_read)
    return variant_reads

def gather_reads_for_single_variant(
        samfile,
        variant,
        chromosome=None,
        use_duplicate_reads=USE_DUPLICATE_READS,
        use_secondary_alignments=USE_SECONDARY_ALIGNMENTS,
        min_mapping_quality=MIN_READ_MAPPING_QUALITY):
    """
    Find reads in the given SAM/BAM file which overlap the given variant, filter
    to only include those which agree with the variant's nucleotide(s), and turn
    them into a list of VariantRead objects.

    Parameters
    ----------
    samfile : pysam.AlignmentFile

    variant : varcode.Variant

    chromosome : str

    use_duplicate_reads : bool
        Should we use reads that have been marked as PCR duplicates

    use_secondary_alignments : bool
        Should we use reads at locations other than their best alignment

    min_mapping_quality : int
        Drop reads below this mapping quality

    Returns list of VariantRead objects.
    """
    logging.info("Gathering reads for %s" % variant)
    if chromosome is None:
        chromosome = variant.contig

    logging.info("Gathering variant reads for variant %s (chromosome=%s)" % (
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

    reads = list(read_at_locus_generator(
        samfile=samfile,
        chromosome=chromosome,
        base1_position_before_variant=base1_position_before_variant,
        base1_position_after_variant=base1_position_after_variant,
        use_duplicate_reads=use_duplicate_reads,
        use_secondary_alignments=use_secondary_alignments,
        min_mapping_quality=min_mapping_quality))

    logging.info(
        "Found %d reads at locus overlapping %s" % (len(reads), variant))

    variant_reads = list(
        variant_reads_from_reads_at_locus(
            reads=reads,
            ref=ref,
            alt=alt))
    logging.info("Found %d VariantReads for variant %s" % (
        len(variant_reads),
        variant))
    return variant_reads

def variant_reads_generator(
        variants,
        samfile,
        use_duplicate_reads=USE_DUPLICATE_READS,
        use_secondary_alignments=USE_SECONDARY_ALIGNMENTS,
        min_mapping_quality=MIN_READ_MAPPING_QUALITY):
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
            logging.warn(
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
        yield variant, variant_reads


def variant_reads_dataframe(
        variants,
        samfile,
        use_duplicate_reads=USE_DUPLICATE_READS,
        use_secondary_alignments=USE_SECONDARY_ALIGNMENTS,
        min_mapping_quality=MIN_READ_MAPPING_QUALITY):
    """
    Creates a DataFrame containing sequences around variant from variant
    collection and BAM/SAM file.

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
    df_builder = DataFrameBuilder(VariantRead)
    for variant, variant_reads in variant_reads_generator(
            variants,
            samfile,
            use_duplicate_reads=use_duplicate_reads,
            use_secondary_alignments=use_secondary_alignments,
            min_mapping_quality=min_mapping_quality):
        for variant_read in variant_reads:
            df_builder.add(variant, variant_read)
    return df_builder.to_dataframe()
