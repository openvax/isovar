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
This module wraps pysam and gives us a view of any reads overlapping
a variant locus which includes offsets into the read sequence & qualities
for extracting variant nucleotides.
"""

from __future__ import print_function, division, absolute_import

from .value_object import ValueObject


class LocusRead(ValueObject):
    """
    Minimal set of information extracted from SAM/BAM alignment file at a particular
    locus to later figure out the allele at this locus.
    """
    __slots__ = [
        "name",
        "sequence",
        "reference_positions",
        "quality_scores",
        "reference_base0_start_inclusive",
        "reference_base0_end_exclusive",
        "read_base0_start_inclusive",
        "read_base0_end_exclusive"
    ]

    def __init__(
            self,
            name,
            sequence,
            reference_positions,
            quality_scores,
            reference_base0_start_inclusive,
            reference_base0_end_exclusive,
            read_base0_start_inclusive,
            read_base0_end_exclusive):
        """
        Parameters
        ----------
        name : str
            Fragment name, paired reads from the same fragment will have the same name

        sequence : str
            cDNA sequence

        reference_positions : list of (int or None)
            For every base in the sequence, which base-1 reference position
            does it map to, or None if the read base is an insertion or soft-clipped

        quality_scores : array of int
            Base qualities for every character in the sequence

        reference_base0_start_inclusive : int
            Start index of reference locus which is overlapped
            by this read (base 0, inclusive)

        reference_base0_end_exclusive : int
            End index of reference locus which is overlapped
            by this read (base 0, exclusive)

        read_base0_start_inclusive : int or None
            Start index of base in read which corresponds to
            start of reference locus (if it's mapped)

        read_base0_end_exclusive : int or None
            End index after last base in sequence which
            corresponds to reference locus (if it's mapped)
        """
        ######################################################################
        # When can the start or end of the read interval be None?
        # ---------------------------------------------------------
        # If a locus goes from  [x:y) on the reference chromosome then there are
        # a few possibilities for whether x and y will have corresponding mapped
        # positions on a read.
        #
        # If x=y, then we expect the read to either match the reference or have some bases
        # inserted between x and y. In this case, the only reason why x and y wouldn't be mapped
        # is if they occur after the end of a read, but the read still overlaps some
        # part of the interval.
        #
        # If y > x, then we're selecting some non-zero reference bases and either matching them
        # or deleting them. If they CIGAR operation is M (match), then x and y should have
        # corresponding positions on the read (unless, like previously, only part of the
        # interval is covered by a read).
        #
        # In the case of a deletion, however, the selected reference bases do not have
        # corresponding positions on the read.
        #
        # Diagram of a deletion:
        #
        #    REFERENCE
        #    19124 19125 19126 19127 19128 19129 19130 19131 19132 19133 19134
        #      A     C     T     G     G     C     A     T      T    T     T
        #
        # If we're interested in checking whether the sequence 'GCA' is deleted
        # between 19128:19131 then we'll look at an RNA read.
        #
        #    RNA WHICH SUPPORTS REFERENCE
        #    00024 00025 00026 00027 00028 00029 00030 00031 00032 00033 00034
        #      A     C     T     G     G     C     A     T      T    T     T
        #
        #    RNA WHICH SUPPORTS MUTATION
        #    00024 00025 00026 00027 00028 00029 00030 00031 00032 00033 00034
        #      A     C     T     G     T      T    T      T    A     A     A
        #
        # In the RNA read which does not support the mutation (matches reference) the
        # reference index 19128 is mapped to 28 and 19131 is mapped to 31.
        #
        # On the mutant RNA, however, the position 19128, 19129, 19130 are unmapped.
        #
        #    19124 19125 19126 19127 19128 19129 19130 19131 19132 19133 19134
        #      A     C     T     G     G     C     A     T      T    T     T
        #      |     |     |     |                       |      |    |     |
        #      |     |     |     |     *-----------------*      |    |     |
        #      |     |     |     |     |      *-----------------*    |     |
        #      |     |     |     |     |      |    |-----------------*     |
        #      |     |     |     |     |      |    |      *----------------*
        #      |     |     |     |     |      |    |      |
        #    00024 00025 00026 00027 00028 00029 00030 00031 00032 00033 00034
        #      A     C     T     G     T      T    T      T    A     A     A

        self.name = name
        self.sequence = sequence
        self.reference_positions = reference_positions
        self.quality_scores = quality_scores
        self.reference_base0_start_inclusive = reference_base0_start_inclusive
        self.reference_base0_end_exclusive = reference_base0_end_exclusive
        self.read_base0_start_inclusive = read_base0_start_inclusive
        self.read_base0_end_exclusive = read_base0_end_exclusive


