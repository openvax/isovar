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


def overlapping_read_tuple_generator(
        samfile,
        chromosome,
        base0_start,
        base0_end,
        is_del):
    """
    Generator that yields a sequence of tuples from reads which overlap the
    given locus and have a matching alignment on the first position.

    The tuples contain the following information:
        - nucleotide sequence of the read
        - list of reference position alignments for each nucleotide
            (where a non-aligned nucleotide is given by None)

    """
    base1_start = base0_start + 1
    # Let pysam pileup the reads covering our location of interest for us
    for column in samfile.pileup(
            chromosome,
            base0_start,
            base0_end):
        if column.pos != base1_start:
            continue
        for i, pileup_element in enumerate(column.pileups):
            if pileup_element.is_refskip:
                # if read sequence doesn't actually align here, skip it
                continue
            elif pileup_element.is_del and not is_del:
                # if read has a deletion at this location and variant isn't a
                # deletion
                continue
            read = pileup_element.alignment
            reference_positions = read.get_reference_positions(
                full_length=False)
            if base0_start in reference_positions:
                offset = reference_positions.index(base0_start)
                yield (reference_positions, offset, read.query_sequence)


def get_variant_base0_interval(base1_location, ref, alt):
    base0_start = base1_location - 1
    if len(ref) == 0:
        # if the variant is an insertion then we need to check to make sure
        # both sides of the insertion are matches
        base0_end = base1_location + 1
    else:
        base0_end = base0_start + len(ref)
    return base0_start, base0_end


def partitioned_read_sequences_from_tuples(read_tuples, ref, alt):
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
    for (reference_positions, offset, sequence) in read_tuples:
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
            suffix = sequence[offset + 1:]
        else:
            # deletions and substitutions work similarly, we just need
            # all the reference bases to be adjacently aligned
            if any(
                    (
                        reference_positions[offset + i] -
                        reference_positions[offset + i - 1]
                    ) != 1
                    for i in range(1, len(ref))):
                continue
            prefix = sequence[:offset]
            suffix = sequence[offset + len(ref):]
        yield str(prefix, "ascii"), alt, str(suffix, "ascii")

def partition_variant_reads(samfile, chromosome, base1_location, ref, alt):
    base1_location, ref, alt = trim_variant(base1_location, ref, alt)
    base0_start, base0_end = get_variant_base0_interval(
        base1_location=base1_location,
        ref=ref,
        alt=alt)
    read_tuples = overlapping_read_tuple_generator(
        samfile=samfile,
        chromosome=chromosome,
        base0_start=base0_start,
        base0_end=base0_end,
        is_del=len(alt) == 0)
    return list(
        partitioned_read_sequences_from_tuples(read_tuples, ref=ref, alt=alt))
