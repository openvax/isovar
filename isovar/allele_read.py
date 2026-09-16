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
Reads overlapping a locus of interest split into prefix,
allele (ref, alt, or otherwise), and suffix portions
"""

import logging

from .string_helpers import convert_from_bytes_if_necessary, trim_N_nucleotides
from .value_object import ValueObject

logger = logging.getLogger(__name__)


class AlleleRead(ValueObject):
    """
    Compact representation of a read at a locus. In addition to its sequence,
    it retains gapless genomic blocks and observed splice junctions so it can
    be classified against annotated transcript paths before assembly.

    When overlapping mates from the same fragment are merged upstream,
    `source_read_count` retains the number of raw reads represented by this
    fragment-level allele observation.

    ``source_alignments`` preserves the immutable segment/placement identities
    from SAM. It defaults to empty for existing manually constructed reads;
    public ``name`` remains the original QNAME string. Counts across collected
    observations must deduplicate segment IDs, not sum ``source_read_count``.
    """
    __slots__ = [
        "prefix",
        "allele",
        "suffix",
        "name",
        "sequence",
        "source_read_count",
        "reference_blocks",
        "splice_junctions",
        "compatible_transcript_ids",
        "source_alignments",
    ]

    def __init__(
            self,
            prefix,
            allele,
            suffix,
            name,
            source_read_count=1,
            reference_blocks=(),
            splice_junctions=(),
            compatible_transcript_ids=None,
            source_alignments=()):
        self.prefix = prefix
        self.allele = allele
        self.suffix = suffix
        self.name = name
        self.sequence = prefix + allele + suffix
        self.source_read_count = source_read_count
        self.source_alignments = tuple(source_alignments)
        self.reference_blocks = tuple(reference_blocks)
        self.splice_junctions = tuple(splice_junctions)
        self.compatible_transcript_ids = (
            None
            if compatible_transcript_ids is None
            else frozenset(compatible_transcript_ids))

    def __len__(self):
        return len(self.sequence)

    def with_compatible_transcript_ids(self, transcript_ids):
        """Return this observation classified against a transcript set."""
        return AlleleRead(
            prefix=self.prefix,
            allele=self.allele,
            suffix=self.suffix,
            name=self.name,
            source_read_count=self.source_read_count,
            reference_blocks=self.reference_blocks,
            splice_junctions=self.splice_junctions,
            compatible_transcript_ids=transcript_ids,
            source_alignments=self.source_alignments)

    @staticmethod
    def _reference_blocks(reference_positions):
        """Compress mappings into query/genomic 0-based half-open blocks."""
        blocks = []
        query_start = reference_start = previous = None
        for query_position, position in enumerate(reference_positions):
            if position is None:
                if query_start is not None:
                    blocks.append((
                        query_start,
                        query_position,
                        reference_start,
                        previous + 1))
                    query_start = reference_start = previous = None
            elif query_start is None:
                query_start = query_position
                reference_start = previous = position
            elif position == previous + 1:
                previous = position
            else:
                blocks.append((
                    query_start,
                    query_position,
                    reference_start,
                    previous + 1))
                query_start = query_position
                reference_start = previous = position
        if query_start is not None:
            blocks.append((
                query_start,
                len(reference_positions),
                reference_start,
                previous + 1))
        return tuple(blocks)

    @staticmethod
    def _slice_reference_blocks(reference_blocks, start, end):
        """Clip whole-read blocks to a query interval and rebase its offsets."""
        clipped = []
        for query_start, query_end, reference_start, _ in reference_blocks:
            clipped_start = max(query_start, start)
            clipped_end = min(query_end, end)
            if clipped_start >= clipped_end:
                continue
            clipped_reference_start = reference_start + clipped_start - query_start
            clipped.append((
                clipped_start - start,
                clipped_end - start,
                clipped_reference_start,
                clipped_reference_start + clipped_end - clipped_start,
            ))
        return tuple(clipped)

    @classmethod
    def from_locus_read(cls, locus_read):
        """
        Given a single LocusRead object, return either an AlleleRead or None

        Parameters
        ----------
        locus_read : LocusRead
            Read which overlaps a variant locus but doesn't necessarily contain the
            alternate nucleotides
        """
        sequence = locus_read.sequence
        read_name = locus_read.name

        reference_base0_start_inclusive = locus_read.reference_base0_start_inclusive
        reference_base0_end_exclusive = locus_read.reference_base0_end_exclusive

        read_base0_start_inclusive = locus_read.read_base0_start_inclusive
        read_base0_end_exclusive = locus_read.read_base0_end_exclusive

        if read_base0_start_inclusive is None or read_base0_end_exclusive is None:
            logger.debug(
                "Skipping read '%s' because required bases in reference interval %s:%s aren't mapped",
                read_name,
                reference_base0_start_inclusive,
                reference_base0_end_exclusive)
            return None

        n_ref_bases = reference_base0_end_exclusive - reference_base0_start_inclusive

        insertion = (n_ref_bases == 0)

        # ReadCollector may canonicalize equivalent indels to a logical query
        # interval even when the aligner placed the indel elsewhere in a repeat.
        nucleotides_at_variant_locus = convert_from_bytes_if_necessary(
            sequence[read_base0_start_inclusive:read_base0_end_exclusive])

        if "N" in nucleotides_at_variant_locus:
            logger.debug(
                "Skipping read '%s', found N nucleotides at variant locus",
                read_name)
            return None

        prefix = convert_from_bytes_if_necessary(sequence[:read_base0_start_inclusive])
        suffix = convert_from_bytes_if_necessary(sequence[read_base0_end_exclusive:])

        prefix, suffix = trim_N_nucleotides(prefix, suffix)

        retained_start = read_base0_start_inclusive - len(prefix)
        retained_end = read_base0_end_exclusive + len(suffix)
        if hasattr(locus_read, "reference_blocks"):
            reference_blocks = cls._slice_reference_blocks(
                locus_read.reference_blocks, retained_start, retained_end)
        else:
            reference_blocks = cls._reference_blocks(
                locus_read.reference_positions[retained_start:retained_end])
        observed_block_gaps = {
            (left[3], right[2])
            for left, right in zip(reference_blocks, reference_blocks[1:])
        }
        splice_junctions = tuple(
            junction
            for junction in locus_read.splice_junctions
            if junction in observed_block_gaps)

        return AlleleRead(
            prefix,
            nucleotides_at_variant_locus,
            suffix,
            name=read_name,
            source_read_count=locus_read.source_read_count,
            reference_blocks=reference_blocks,
            splice_junctions=splice_junctions,
            source_alignments=locus_read.source_alignments)


# Branch-local transcript classifications are not distinct observations.
# source_alignments DOES participate in identity: retain competing placements
# as hypotheses, then deduplicate support by segment rather than object count.
AlleleRead._fields = tuple(
    field
    for field in AlleleRead._fields
    if field not in {
        "reference_blocks",
        "splice_junctions",
        "compatible_transcript_ids",
    }
)
