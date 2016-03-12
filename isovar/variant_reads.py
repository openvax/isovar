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

import logging
from collections import namedtuple, OrderedDict

from pandas import DataFrame

from .overlapping_reads import gather_overlapping_reads

VariantRead = namedtuple(
    "VariantRead", "prefix variant suffix name")

def trim_variant(location, ref, alt):
    """
    Trims common prefixes from the ref and alt sequences

    Parameters
    ----------
    location : int
        Position (starting from 1) on some chromosome

    ref : str
        Reference nucleotides

    alt : str
        Alternate (mutant) nucleotide

    Returns adjusted triplet (location, ref, alt)
    """
    if len(alt) > 0 and ref.startswith(alt):
        # if alt is a prefix of the ref sequence then we actually have a
        # deletion like:
        #   g.10 GTT > GT
        # which can be trimmed to
        #   g.12 'T'>''
        ref = ref[len(alt):]
        location += len(alt)
        alt = ""
    if len(ref) > 0 and alt.startswith(ref):
        # if ref sequence is a prefix of the alt sequence then we actually have
        # an insertion like:
        #   g.10 GT>GTT
        # which can be trimmed to
        #   g.11 ''>'T'
        # Note that we are selecting the position *before* the insertion
        # (as an arbitrary convention)
        alt = alt[len(ref):]
        location += len(ref) - 1
        ref = ""
    return location, ref, alt

def get_variant_base0_interval(base1_location, ref, alt):
    if len(ref) == 0:
        # if the variant is an insertion then we need to check to make sure
        # both sides of the insertion are matches
        base0_start = base1_location - 1
        base0_end = base1_location + 1
    elif len(alt) == 0:
        # if we're deleting from the sequence, then move the interval
        # to the base behind the deletion so that we have an aligned
        # nucleotide to use from string slicing
        base0_start = base1_location - 2
        base0_end = base1_location + 1
    else:
        base0_start = base1_location - 1
        base0_end = base0_start + len(ref)
    return base0_start, base0_end


def variant_reads_from_overlapping_reads(overlapping_reads, ref, alt):
    """
    Given a collection of pysam.AlignedSegment objects, generates a
    sequence of VariantRead objects (which are split into prefix/variant/suffix
    nucleotides).

    Parameters
    ----------
    overlapping_reads : list of pysam.AlignedSegment

    ref : str
        Reference sequence of the variant (empty for insertions)

    alt : str
        Alternate sequence of the variant (empty for deletions)

    Returns a sequence of VariantRead objects.
    """
    for read in overlapping_reads:
        reference_positions = read.reference_positions
        offset = read.locus_offset
        sequence = read.sequence
        if len(ref) == 0:
            # insertions require a sequence of non-aligned bases
            # followed by the subsequence reference position
            if len(reference_positions) < offset + len(alt) + 1:
                continue
            ref_pos = reference_positions[offset]
            insert_positions = reference_positions[offset + 1:offset + len(alt)]
            if any(insert_pos is not None for insert_pos in insert_positions):
                # all these inserted nucleotides should *not* align to the
                # reference
                continue
            next_ref_pos = reference_positions[offset + 1]
            if next_ref_pos != ref_pos + 1:
                continue
            prefix = sequence[:offset + 1]
            suffix = sequence[offset + len(alt) + 1:]
        elif len(alt) == 0:
            if len(reference_positions) < offset + 2:
                # if we're missing the position after the deletion then
                # skip this read
                continue
            ref_pos_before = reference_positions[offset]
            ref_pos_after = reference_positions[offset + 1]
            if ref_pos_after - ref_pos_before - 1 != len(ref):
                # if the number of nucleotides skipped isn't the same
                # as the number deleted in the variant then
                # don't use this read
                continue
            prefix = sequence[:offset + 1]
            suffix = sequence[offset + 1:]
        else:
            # deletions and substitutions work similarly, we just need
            # all the reference bases to be adjacently aligned
            ref_pos_start = reference_positions[offset]
            ref_pos_end = reference_positions[offset + len(ref)]
            if ref_pos_end - ref_pos_start != len(ref):
                continue
            prefix = sequence[:offset]
            suffix = sequence[offset + len(ref):]
        if isinstance(prefix, bytes):
            prefix = str(prefix, "ascii")
        if isinstance(suffix, bytes):
            suffix = str(suffix, "ascii")
        yield VariantRead(prefix, alt, suffix, name=read.name)

def gather_reads_for_single_variant(
        samfile,
        chromosome,
        base1_location,
        ref,
        alt):
    """
    Find reads in the given SAM/BAM file which overlap the given variant, filter
    to only include those which agree with the variant's nucleotide(s), and turn
    them into a list of VariantRead objects.

    Parameters
    ----------
    samfile : pysam.AlignmentFile

    chromosome : str

    base1_location : int

    ref : str
        Reference nucleotides

    alt : str
        Variant nucleotides

    Returns list of VariantRead objects.
    """
    base1_location, ref, alt = trim_variant(base1_location, ref, alt)

    logging.info("Gathering variant reads for variant %s:%s '%s'>'%s'" % (
        chromosome,
        base1_location,
        ref,
        alt))
    base0_start, base0_end = get_variant_base0_interval(
        base1_location=base1_location,
        ref=ref,
        alt=alt)
    overlapping_reads = gather_overlapping_reads(
        samfile=samfile,
        chromosome=chromosome,
        base0_start=base0_start,
        base0_end=base0_end,
        is_del=len(alt) == 0)
    return list(
        variant_reads_from_overlapping_reads(
            overlapping_reads=overlapping_reads,
            ref=ref,
            alt=alt))


def variant_reads_generator(variants, samfile):
    """
    Generates sequence of tuples, each containing a variant paired with
    a list of VariantRead objects.

    Parameters
    ----------
    variants : varcode.VariantCollection

    samfile : pysam.AlignmentFile
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
                    chromosome, variant, samfile))
            continue

        variant_reads = gather_reads_for_single_variant(
            samfile=samfile,
            chromosome=chromosome,
            base1_location=variant.start,
            ref=variant.ref,
            alt=variant.alt)
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
