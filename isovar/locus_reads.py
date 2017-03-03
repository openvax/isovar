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
This module wraps pysam and gives us a view of any reads overlapping
a variant locus which includes offsets into the read sequence & qualities
for extracting variant nucleotides.
"""

from __future__ import print_function, division, absolute_import
import logging

from .default_parameters import (
    MIN_READ_MAPPING_QUALITY,
    USE_DUPLICATE_READS,
    USE_SECONDARY_ALIGNMENTS,
)
from .common import list_to_string
from .dataframe_builder import DataFrameBuilder
from .value_object import ValueObject

logger = logging.getLogger(__name__)

class LocusRead(ValueObject):
    __slots__ = [
        "name",
        "sequence",
        "reference_positions",
        "quality_scores",
    ]

    def __init__(
            self,
            name,
            sequence,
            reference_positions,
            quality_scores):
        self.name = name
        self.sequence = sequence
        self.reference_positions = reference_positions
        self.quality_scores = quality_scores

    def translate_interval_to_read():
        pass

class LocusReadCollector(object):
    """
    Holds to an AlignmentFile object for a SAM/BAM file and then
    generates lists of LocusRead objects at the specified locus.
    """
    def __init__(
            self,
            samfile,
            use_duplicate_reads=USE_DUPLICATE_READS,
            use_secondary_alignments=USE_SECONDARY_ALIGNMENTS,
            min_mapping_quality=MIN_READ_MAPPING_QUALITY):
        self.samfile = samfile
        self.use_duplicate_reads = use_duplicate_reads
        self.use_secondary_alignments = use_secondary_alignments
        self.min_mapping_quality = min_mapping_quality

    def _convert_aligned_segment(
            self,
            aligned_segment):
        """
        Converts pysam.AlignedSegment object to LocusRead. Can also return None
        if aligned segment doesn't meet the filtering criteria such as
        minimum mapping quality.
        """
        # For future reference,  may get overlapping reads
        # which can be identified by having the same name
        name = aligned_segment.query_name

        if name is None:
            logger.warn(
                "Read missing name: %s", aligned_segment)
            return None

        if aligned_segment.is_unmapped:
            logger.warn(
                "Unmapped read: %s", aligned_segment)
            return None

        if aligned_segment.is_secondary and not self.use_secondary_alignments:
            logger.debug("Skipping secondary alignment of read '%s'", name)
            return None

        if aligned_segment.is_duplicate and not self.use_duplicate_reads:
            logger.debug("Skipping duplicate read '%s'", name)
            return None

        mapping_quality = aligned_segment.mapping_quality

        missing_mapping_quality = mapping_quality is None

        if self.min_mapping_quality > 0 and missing_mapping_quality:
            logger.debug("Skipping read '%s' due to missing MAPQ", name)
            return None
        elif mapping_quality < self.min_mapping_quality:
            logger.debug(
                "Skipping read '%s' due to low MAPQ %d (min=%d)",
                name,
                mapping_quality,
                self.min_mapping_quality)
            return None

        sequence = aligned_segment.query_sequence

        if sequence is None:
            logger.warn("Read '%s' missing sequence", name)
            return None

        base_qualities = aligned_segment.query_qualities

        if base_qualities is None:
            logger.warn("Read '%s' missing base qualities", name)
            return None

        # Documentation for pysam.AlignedSegment.get_reference_positions:
        # ------------------------------------------------------------------
        # By default, this method only returns positions in the reference
        # that are within the alignment. If full_length is set, None values
        # will be included for any soft-clipped or unaligned positions
        # within the read. The returned list will thus be of the same length
        # as the read.
        #
        # Source:
        # http://pysam.readthedocs.org/en/latest/
        # api.html#pysam.AlignedSegment.get_reference_positions
        #
        # We want a None value for every read position that does not have a
        # corresponding reference position.
        reference_positions = aligned_segment.get_reference_positions(
            full_length=True)
        """
        # pysam uses base-0 positions everywhere except region strings
        # Source:
        # http://pysam.readthedocs.org/en/latest/faq.html#pysam-coordinates-are-wrong
        if base0_locus_start not in reference_positions:
            logger.debug(
                "Skipping read '%s' because first position %d not mapped",
                name,
                base0_position_before_variant)
            return None
        else:
            base0_read_position_before_variant = reference_positions.index(
                base0_position_before_variant)

        if base0_position_after_variant not in reference_positions:
            logger.debug(
                "Skipping read '%s' because last position %d not mapped",
                name,
                base0_position_after_variant)
            return None
        else:
            base0_read_position_after_variant = reference_positions.index(
                base0_position_after_variant)
        """

        if isinstance(sequence, bytes):
            sequence = sequence.decode('ascii')

        return LocusRead(
            name=name,
            sequence=sequence,
            reference_positions=reference_positions,
            quality_scores=base_qualities)

    def get_locus_reads(
            self,
            chromosome,
            base0_locus_start,
            base0_locus_end):
        """
        Generator that yields a sequence of ReadAtLocus records for reads which
        overlap the requested locus. The actual work to figure out if what's at the
        locus matches a variant happens later.

        Parameters
        ----------
        samfile : pysam.AlignmentFile

        chromosome : str

        base0_locus_start : int
            Interbase position before first reference nucleotide

        base0_locus_end : int
            Interbase position after last reference nucleotide

        Returns list of LocusRead objects
        """
        logger.debug(
            "Gathering reads at locus chr=%s, interbase start=%d, interbase end=%d",
            chromosome,
            base0_locus_start,
            base0_locus_end)

        reads = []
        # iterate over any pysam.AlignedSegment objects which overlap this locus
        for aligned_segment in self.samfile.fetch(
                chromosome, base0_locus_start, base0_locus_end):
            locus_read = self._convert_aligned_segment(aligned_segment)
            reference_positions = locus_read.reference_positions

            aligned_positions = [
                pos in reference_positions
                for pos in range(base0_locus_start, base0_locus_end)]

            if len(aligned_positions) > 0 and all(aligned_positions):
                # if we're looking at positions in the reference and they
                # are all mapped in this read then we're good
                reads.append(locus_read)
            else:
                # if we're looking for an insertion between two bases or a
                # deletion then we want to at least make sure that one of the bases
                # to the left or right is mapped where we expect

                before = base0_locus_start - 1
                after = base0_locus_end - 1
                adjacent_position_aligned = (
                    before in reference_positions or
                    after in reference_positions)
                if not adjacent_position_aligned:
                    logger.info(
                        "Skipping %s since no adjacent positions are mapped",
                        locus_read)
                    continue

                if base0_locus_start == base0_locus_end:
                    # read which can be used to look for an insertion
                    reads.append(locus_read)
                elif not any(aligned_positions):
                    # read which is probably evidence of a deletion
                    reads.append(locus_read)
                else:
                    logger.info(
                        "Not all required positions were mapped in %s",
                        locus_read)

        logger.info(
            "Found %d reads overlapping locus %s[%d:%d]",
            len(reads),
            chromosome,
            base0_locus_start,
            base0_locus_end)
        # TODO: de-duplicate reads by name
        # since overlapping mate pairs will have the same name, then...
        #
        # TODO: combine overlapping mate pairs into single sequence before
        # de-duplication
        return reads

def locus_reads_dataframe(
        chromosome, base0_locus_start, base0_locus_end, **generator_kwargs):
    """
    Traverse a BAM file to find all the reads overlapping a specified locus.

    Parameters are the same as those for read_locus_generator.
    """
    df_builder = DataFrameBuilder(
        LocusRead,
        variant_columns=False,
        converters={
            "reference_positions": list_to_string,
            "quality_scores": list_to_string,
        })

    collector = LocusReadCollector(**generator_kwargs)
    reads = collector.get_locus_reads(
        chromosome=chromosome,
        base0_locus_start=base0_locus_start,
        base0_locus_end=base0_locus_end)
    for locus_read in reads:
        df_builder.add(variant=None, element=locus_read)
    return df_builder.to_dataframe()
