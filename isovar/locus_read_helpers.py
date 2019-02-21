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
Helper functions for extracting LocusRead objects from read alignment files at
specified locations.
"""

from __future__ import print_function, division, absolute_import

from .default_parameters import (
    MIN_READ_MAPPING_QUALITY,
    USE_DUPLICATE_READS,
    USE_SECONDARY_ALIGNMENTS,
)
from .locus_read import  LocusRead
from .logging import get_logger

logger = get_logger(__name__)

def pileup_reads_at_position(samfile, chromosome, base0_position):
    """
    Returns a pileup column at the specified position. Unclear if a function
    like this is hiding somewhere in pysam API.
    """

    # TODO: I want to pass truncate=True, stepper="all"
    # but for some reason I get this error:
    #      pileup() got an unexpected keyword argument 'truncate'
    # ...even though these options are listed in the docs for pysam 0.9.0
    #
    for column in samfile.pileup(
            chromosome,
            start=base0_position,
            end=base0_position + 1):

        if column.pos != base0_position:
            # if this column isn't centered on the base before the
            # variant then keep going
            continue

        return column.pileups

    # if we get to this point then we never saw a pileup at the
    # desired position
    return []


def locus_read_generator(
        samfile,
        chromosome,
        base1_position_before_variant,
        base1_position_after_variant,
        use_duplicate_reads=USE_DUPLICATE_READS,
        use_secondary_alignments=USE_SECONDARY_ALIGNMENTS,
        min_mapping_quality=MIN_READ_MAPPING_QUALITY):
    """
    Generator that yields a sequence of ReadAtLocus records for reads which
    contain the positions before and after a variant. The actual work to figure
    out if what's between those positions matches a variant happens later in
    the `variant_reads` module.

    Parameters
    ----------
    samfile : pysam.AlignmentFile

    chromosome : str

    base1_position_before_variant : int
        Genomic position of reference nucleotide before a variant

    base1_position_after_variant : int
        Genomic position of reference nucleotide before a variant

    use_duplicate_reads : bool
        By default, we're ignoring any duplicate reads

    use_secondary_alignments : bool
        By default we are using secondary alignments, set this to False to
        only use primary alignments of reads.

    min_mapping_quality : int
        Drop reads below this mapping quality

    Yields ReadAtLocus objects
    """
    logger.debug(
        "Gathering reads at locus %s: %d-%d",
        chromosome,
        base1_position_before_variant,
        base1_position_after_variant)
    base0_position_before_variant = base1_position_before_variant - 1
    base0_position_after_variant = base1_position_after_variant - 1

    count = 0

    for aligned_segment in samfile.fetch(
                contig=chromosome,
                start=base0_position_before_variant,
                stop=base0_position_after_variant):
        read = LocusRead.from_pysam_aligned_segment(
            aligned_segment,
            base0_position_before_variant=base0_position_before_variant,
            base0_position_after_variant=base0_position_after_variant,
            use_secondary_alignments=use_secondary_alignments,
            use_duplicate_reads=use_duplicate_reads,
            min_mapping_quality=min_mapping_quality)

        if read is not None:
            count += 1
            yield read

    logger.info(
        "Found %d reads overlapping locus %s: %d-%d",
        count,
        chromosome,
        base1_position_before_variant,
        base1_position_after_variant)
