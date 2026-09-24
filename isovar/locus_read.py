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

from .value_object import ValueObject


class LocusRead(ValueObject):
    """
    Minimal set of information extracted from SAM/BAM alignment file at a particular
    locus to later figure out the allele at this locus.

    Overlapping paired-end reads from the same fragment may be merged into a
    single LocusRead. In that case `source_read_count` records how many
    sequenced segments were collapsed into this fragment-level view.
    """
    __slots__ = [
        "name",
        "sequence",
        "reference_positions",
        "quality_scores",
        "source_read_count",
        "reference_base0_start_inclusive",
        "reference_base0_end_exclusive",
        "read_base0_start_inclusive",
        "read_base0_end_exclusive",
        "splice_junctions",
        "source_alignments",
        "is_primary",
        "source_alignment_paths",
        "source_read_views",
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
            read_base0_end_exclusive,
            source_read_count=1,
            splice_junctions=(),
            source_alignments=(),
            is_primary=False,
            source_alignment_paths=(),
            source_read_views=()):
        """
        Parameters
        ----------
        name : str
            Fragment name, paired reads from the same fragment will have the same name

        sequence : str
            cDNA sequence

        reference_positions : list of (int or None)
            For every base in the sequence, which base-0 reference position
            does it map to, or None if the read base is an insertion or soft-clipped

        quality_scores : sequence of int or None
            Base qualities for every character in the sequence. None marks
            unavailable quality, not a fabricated low/high Phred score.

        source_read_count : int
            Number of raw reads represented by this LocusRead. Usually 1, but
            overlapping paired-end mates from the same fragment may be merged
            into a single LocusRead while keeping this count.

        splice_junctions : sequence of (int, int)
            Reference skips from CIGAR ``N`` operations, represented as
            0-based half-open genomic intervals ``(intron_start, intron_end)``.

        source_alignments : tuple
            Immutable (segment identity, alignment identity) pairs collected
            from SAM; see ``read_identity``. Empty for legacy caller-created
            objects. Alternative placements do not identify extra reads.

        source_alignment_paths : tuple
            Optional immutable SAM chimeric-path declarations, keyed by source
            segment. Used only for phasing; see ``chimeric_alignment``.

        source_read_views : tuple
            Optional (source alignment, ReadSequenceView) provenance for opt-in
            end inference. Views index original SAM SEQ, not merged coordinates.

        is_primary : bool
            Whether this view contains primary alignments. Mate merging also
            requires explicit, complementary segment flags and a common read
            group. Missing provenance is not evidence that two views are mates.

        reference_base0_start_inclusive : int
            Start index of reference locus which is overlapped
            by this read (base 0, inclusive)

        reference_base0_end_exclusive : int
            End index of reference locus which is overlapped
            by this read (base 0, exclusive)

        read_base0_start_inclusive, read_base0_end_exclusive : int or None
            Query interval for the allele. Both are None when the locus isn't
            mapped on this read, including an unassigned insertion boundary.
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
        # or deleting them. If the CIGAR operation is M (match), then x and y should have
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
        self.source_read_count = source_read_count
        self.reference_base0_start_inclusive = reference_base0_start_inclusive
        self.reference_base0_end_exclusive = reference_base0_end_exclusive
        self.read_base0_start_inclusive = read_base0_start_inclusive
        self.read_base0_end_exclusive = read_base0_end_exclusive
        self.splice_junctions = tuple(splice_junctions)
        self.source_alignments = tuple(source_alignments)
        self.is_primary = is_primary
        self.source_alignment_paths = tuple(source_alignment_paths)
        self.source_read_views = tuple(source_read_views)
