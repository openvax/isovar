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

from collections import namedtuple

from .overlapping_reads import gather_overlapping_reads

VariantRead = namedtuple(
    "VariantRead", "prefix variant suffix name")

def trim_variant(location, ref, alt):
    """Trims common prefixes from the ref and alt sequences"""
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
    Parameters
    ----------
    read_tuples : sequence
        Each tuple contains the following elements:
        1) list of aligned reference positions
        2) offset into the reference (where the variant aligned to)
        3) nucleotide sequence of the read

    ref : str
        Reference sequence of the variant (empty for insertions)

    alt : str
        Alternate sequence of the variant (empty for deletions)

    Returns a sequence of (prefix, alt, suffix) string tuples.
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

def gather_variant_reads(samfile, chromosome, base1_location, ref, alt):
    base1_location, ref, alt = trim_variant(base1_location, ref, alt)
    base0_start, base0_end = get_variant_base0_interval(
        base1_location=base1_location,
        ref=ref,
        alt=alt)
    reads = gather_overlapping_reads(
        samfile=samfile,
        chromosome=chromosome,
        base0_start=base0_start,
        base0_end=base0_end,
        is_del=len(alt) == 0)
    return list(variant_reads_from_overlapping_reads(reads, ref=ref, alt=alt))
