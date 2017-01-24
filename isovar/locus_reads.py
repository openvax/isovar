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
        "base0_locus_start",
        "base0_locus_end"
    ]

    def __init__(
            self,
            name,
            sequence,
            reference_positions,
            quality_scores,
            base0_locus_start,
            base0_locus_end):
        self.name = name
        self.sequence = sequence
        self.reference_positions = reference_positions
        self.quality_scores = quality_scores
        self.base0_locus_start = base0_locus_start
        self.base0_locus_end = base0_locus_end

    @classmethod
    def from_pysam_pileup_element(
            cls,
            pileup_element,
            base0_position_before_variant,
            base0_position_after_variant,
            use_secondary_alignments,
            use_duplicate_reads,
            min_mapping_quality):
        """
        Parameters
        ----------
        pileup_element : pysam.PileupRead

        base0_position_before_variant : int

        base0_position_after_variant : int

        use_secondary_alignments : bool

        use_duplicate_reads : bool

        min_mapping_quality : int

        Returns LocusRead or None
        """
        read = pileup_element.alignment

        # For future reference,  may get overlapping reads
        # which can be identified by having the same name
        name = read.query_name

        if name is None:
            logger.warn(
                "Read missing name at position %d",
                base0_position_before_variant + 1)
            return None

        if read.is_unmapped:
            logger.warn(
                "How did we get unmapped read '%s' in a pileup?", name)
            return None

        if pileup_element.is_refskip:
            # if read sequence doesn't actually align to the reference
            # base before a variant, skip it
            logger.debug("Skipping pileup element with CIGAR alignment N (intron)")
            return None
        elif pileup_element.is_del:
            logger.debug(
                "Skipping deletion at position %d (read name = %s)",
                base0_position_before_variant + 1,
                name)
            return None

        if read.is_secondary and not use_secondary_alignments:
            logger.debug("Skipping secondary alignment of read '%s'", name)
            return None

        if read.is_duplicate and not use_duplicate_reads:
            logger.debug("Skipping duplicate read '%s'", name)
            return None

        mapping_quality = read.mapping_quality

        missing_mapping_quality = mapping_quality is None

        if min_mapping_quality > 0 and missing_mapping_quality:
            logger.debug("Skipping read '%s' due to missing MAPQ", name)
            return None
        elif mapping_quality < min_mapping_quality:
            logger.debug(
                "Skipping read '%s' due to low MAPQ: %d < %d",
                read.mapping_quality,
                mapping_quality,
                min_mapping_quality)
            return None

        sequence = read.query_sequence

        if sequence is None:
            logger.warn("Read '%s' missing sequence", name)
            return None

        base_qualities = read.query_qualities

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
        reference_positions = read.get_reference_positions(
            full_length=True)

        # pysam uses base-0 positions everywhere except region strings
        # Source:
        # http://pysam.readthedocs.org/en/latest/faq.html#pysam-coordinates-are-wrong
        if base0_position_before_variant not in reference_positions:
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

        if isinstance(sequence, bytes):
            sequence = sequence.decode('ascii')

        return cls(
            name=name,
            sequence=sequence,
            reference_positions=reference_positions,
            quality_scores=base_qualities,
            base0_read_position_before_variant=base0_read_position_before_variant,
            base0_read_position_after_variant=base0_read_position_after_variant)

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

    def _convert_aligned_segment(self, aligned_segment):
        """
        Converts pysam.AlignedSegment object to LocusRead. Can also return None
        if aligned segment doesn't meet the filtering criteria such as
        minimum mapping quality.
        """
        pass

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
            pass

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
