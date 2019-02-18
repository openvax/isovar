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
Functions for filtering, grouping, and summarizing collections of
AlleleRead objects.
"""

from collections import defaultdict

from .allele_read import AlleleRead
from .common import groupby
from .default_parameters import (
    MIN_READ_MAPPING_QUALITY,
    USE_SECONDARY_ALIGNMENTS,
    USE_DUPLICATE_READS,
)
from .locus_read import locus_read_generator
from .logging import get_logger
from .variant_helpers import trim_variant

logger = get_logger(__name__)


def allele_reads_from_locus_reads(locus_reads, n_ref):
    """
    Given a collection of LocusRead objects, returns a
    list of AlleleRead objects
    (which are split into prefix/allele/suffix nucleotide strings).

    Parameters
    ----------
    locus_reads : sequence of LocusRead records

    n_ref : int
        Number of reference nucleotides affected by variant.

    Generates AlleleRead objects.
    """

    for locus_read in locus_reads:
        allele_read = AlleleRead.from_locus_read(locus_read, n_ref)
        if allele_read is None:
            continue
        else:
            yield allele_read


def reads_overlapping_variant(
        samfile,
        variant,
        chromosome=None,
        use_duplicate_reads=USE_DUPLICATE_READS,
        use_secondary_alignments=USE_SECONDARY_ALIGNMENTS,
        min_mapping_quality=MIN_READ_MAPPING_QUALITY):
    """
    Find reads in the given SAM/BAM file which overlap the given variant and
    return them as a list of AlleleRead objects.

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

    Returns sequence of AlleleRead objects.
    """
    logger.info("Gathering reads for %s", variant)
    if chromosome is None:
        chromosome = variant.contig

    logger.info(
        "Gathering variant reads for variant %s (chromosome = %s, gene names = %s)",
        variant,
        chromosome,
        variant.gene_names)

    base1_position, ref, alt = trim_variant(variant)

    if len(ref) == 0:
        # if the variant is an insertion
        base1_position_before_variant = base1_position
        base1_position_after_variant = base1_position + 1
    else:
        base1_position_before_variant = base1_position - 1
        base1_position_after_variant = base1_position + len(ref)

    locus_reads = locus_read_generator(
        samfile=samfile,
        chromosome=chromosome,
        base1_position_before_variant=base1_position_before_variant,
        base1_position_after_variant=base1_position_after_variant,
        use_duplicate_reads=use_duplicate_reads,
        use_secondary_alignments=use_secondary_alignments,
        min_mapping_quality=min_mapping_quality)

    allele_reads = allele_reads_from_locus_reads(
        locus_reads=locus_reads,
        n_ref=len(ref))

    return allele_reads


def reads_overlapping_variants(variants, samfile, **kwargs):
    """
    Generates sequence of tuples, each containing a variant paired with
    a list of AlleleRead objects.

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
            logger.warning(
                "Chromosome '%s' from variant %s not in alignment file %s",
                variant.contig,
                variant,
                samfile.filename)
            yield variant, []
            continue
        allele_reads = reads_overlapping_variant(
            samfile=samfile,
            chromosome=chromosome,
            variant=variant,
            **kwargs)
        yield variant, allele_reads


def group_reads_by_allele(allele_reads):
    """
    Returns dictionary mapping each allele's nucleotide sequence to a list of
    supporting AlleleRead objects.
    """
    allele_to_reads_dict = defaultdict(list)
    for allele_read in allele_reads:
        allele_to_reads_dict[allele_read.allele].append(allele_read)
    return allele_to_reads_dict


def filter_non_alt_reads_for_variant(variant, allele_reads):
    """
    Given a variant and an unfiltered collection of AlleleRead objects,
    return only the AlleleRead object which match the alt allele of the variant.
    """
    _, _, alt = trim_variant(variant)
    return [read for read in allele_reads if read.allele == alt]


def filter_non_alt_reads_for_variants(variants_and_allele_reads_sequence):
    """
    Given a sequence of variants paired with all of their overlapping reads,
    yields a sequence of variants paired only with reads which contain their
    mutated nucleotide sequence.
    """
    for variant, allele_reads in variants_and_allele_reads_sequence:
        yield variant, filter_non_alt_reads_for_variant(variant, allele_reads)


def reads_supporting_variant(variant, samfile, **kwargs):
    """
    Parameters
    ----------
    variant: varcode.Variant

    samfile:  pysam.AlignmentFile

    Given a variant and a SAM/BAM file, finds all AlleleRead objects overlapping a variant and returns
    those that support the variant's alt allele.
    """
    allele_reads = reads_overlapping_variant(
        variant=variant,
        samfile=samfile,
        **kwargs)
    return filter_non_alt_reads_for_variant(
        variant=variant,
        allele_reads=allele_reads)


def reads_supporting_variants(variants, samfile, **kwargs):
    """
    Given a SAM/BAM file and a collection of variants, generates a sequence
    of variants paired with reads which support each variant.
    """
    for variant, allele_reads in reads_overlapping_variants(
            variants=variants,
            samfile=samfile,
            **kwargs):
        yield variant, filter_non_alt_reads_for_variant(variant, allele_reads)


def get_single_allele_from_reads(allele_reads):
    """
    Given a sequence of AlleleRead objects, which are expected to all have
    the same allele, return that allele.
    """
    allele_reads = list(allele_reads)

    if len(allele_reads) == 0:
        raise ValueError("Expected non-empty list of AlleleRead objects")

    seq = allele_reads[0].allele
    if any(read.allele != seq for read in allele_reads):
        raise ValueError("Expected all AlleleRead objects to have same allele '%s', got %s" % (
            seq, allele_reads))
    return seq


def group_unique_sequences(
        allele_reads,
        max_prefix_size=None,
        max_suffix_size=None):
    """
    Given a list of AlleleRead objects, extracts all unique
    (prefix, allele, suffix) sequences and associate each with a list
    of reads that contained that sequence.
    """
    groups = defaultdict(set)
    for r in allele_reads:
        prefix = r.prefix
        allele = r.allele
        suffix = r.suffix
        if max_prefix_size and len(prefix) > max_prefix_size:
            prefix = prefix[-max_prefix_size:]
        if max_suffix_size and len(suffix) > max_suffix_size:
            suffix = suffix[:max_suffix_size]
        key = (prefix, allele, suffix)
        groups[key].add(r)
    return groups


def count_unique_sequences(
        allele_reads,
        max_prefix_size=None,
        max_suffix_size=None):
    """
    Given a list of AlleleRead objects, extracts all unique
    (prefix, allele, suffix) sequences and associate each with the number
    of reads that contain that sequence.
    """
    groups = group_unique_sequences(
        allele_reads,
        max_prefix_size=max_prefix_size,
        max_suffix_size=max_suffix_size)
    return {
        seq_tuple: len(read_names)
        for (seq_tuple, read_names) in groups.items()
    }


def group_reads_by_allele(allele_reads):
    """
    Returns dictionary mapping each allele's nucleotide sequence to a list of
    supporting AlleleRead objects.
    """
    return groupby(allele_reads, lambda read: read.allele)
