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

from __future__ import print_function, division, absolute_import

from collections import namedtuple

from .logging import create_logger

logger = create_logger(__name__)

ReadAtLocus = namedtuple(
    "ReadAtLocus",
    [
        "sequence",
        "reference_positions",
        "base_qualities",
        "locus_start",
        "locus_end",
        "locus_size",
        "is_deletion",
        "name"
    ])

def gather_overlapping_reads(
        samfile,
        chromosome,
        base1_position,
        n_bases,
        is_deletion,
        padding=1,
        use_duplicate_reads=False,
        use_secondary_alignments=True,
        min_mapping_quality=5):
    """
    Generator that yields a sequence of OverlappingRead records for reads which
    overlap the given locus and have a matching alignment on the first position.

    Parameters
    ----------
    samfile : pysam.AlignmentFile

    chromosome : str

    base1_position : int
        First variant position in a reference locus

    n_bases : int
        Number of bases to select

    is_deletion : bool
        Does this locus originate from a deletion variant?
        (if so, we want to include reads that have deletions at the locus)

    padding : int
        Number of bases to the left and right of a variant that we want a read
        to contain before considering it overlapping.

    use_duplicate_reads : bool
        By default, we're ignoring any duplicate reads

    use_secondary_alignments : bool
        By default we are using secondary alignments, set this to False to
        only use primary alignments of reads.

    min_mapping_quality : int
        Drop reads below this mapping quality

    Yields OverlappingRead objects contain the following information:
        - nucleotide sequence of the read
        - list of reference position alignments for each nucleotide
            (where a non-aligned nucleotide is given by None)
        - offset of the locus within the read

    """

    base0_position = base1_position - 1

    # Let pysam pileup the reads covering our location of interest for us
    #
    # Annoyingly this function takes base-0 intervals but returns
    # columns with base-1 positions
    for column in samfile.pileup(chromosome, base0_position, base1_position):
        if column.pos != base1_position:
            # if this column isn't centered on the first base of the
            # variant then keep going
            continue

        for i, pileup_element in enumerate(column.pileups):
            if pileup_element.is_refskip:
                # if read sequence doesn't actually align here, skip it
                continue
            elif pileup_element.is_del and not is_deletion:
                logger.debug("Skipping deletion")
                # if read has a deletion at this location and variant isn't a
                # deletion
                #
                # TODO: how to ensure that *all* of n_bases are deleted?
                continue

            read = pileup_element.alignment

            logger.debug(read)

            # For future reference,  may get overlapping reads
            # which can be identified by having the same name
            name = read.query_name

            if name is None:
                logger.warn("Read at locus %s:%d missing name" % (
                    chromosome,
                    base0_position + 1))
                continue

            if read.is_unmapped:
                logger.warn(
                    "How did we get unmapped read '%s' in a pileup?" % (name,))
                continue

            if read.is_secondary and not use_secondary_alignments:
                logger.debug("Skipping secondary alignment of read '%s'")
                continue

            if read.is_duplicate and not use_duplicate_reads:
                logger.debug("Skipping duplicate read '%s'" % name)
                continue

            if read.mapping_quality < min_mapping_quality:
                logger.debug("Skipping read '%s' due to low mapq: %d < %d" % (
                    read.mapping_quality, min_mapping_quality))
                continue

            #
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
            reference_positions = read.get_reference_positions(
                full_length=False)

            # if we're looking for deleted bases then inspect the bases to the
            # left and right of the deletion. Otherwise we should find
            # the reference bases directly on the read.
            genomic_locus_start = (
                base0_position - 1 if is_deletion else base0_position)
            genomic_locus_end = (
                base0_position + n_bases if is_deletion
                else base0_position + n_bases - 1)

            # pysam uses base-0 positions everywhere except region strings
            # Source:
            # http://pysam.readthedocs.org/en/latest/faq.html#pysam-coordinates-are-wrong
            if genomic_locus_start not in reference_positions:
                logger.debug(
                    "Skipping read '%s' because first position %d not mapped" % (
                        name,
                        genomic_locus_start))
                continue
            else:
                locus_start_in_read = reference_positions.index(genomic_locus_start)

            if genomic_locus_end not in reference_positions:
                logger.debug(
                    "Skipping read '%s' because last position %d not mapped" % (
                        name,
                        genomic_locus_end))
                continue
            else:
                locus_end_in_read = reference_positions.index(genomic_locus_end)

            sequence = read.query_sequence
            if sequence is None:
                logger.warn("Read '%s' missing sequence")
                continue

            base_qualities = read.query_qualities

            if base_qualities is None:
                logger.warn("Read '%s' missing base qualities" % (name,))
                continue

            yield ReadAtLocus(
                sequence=sequence,
                reference_positions=reference_positions,
                base_qualities=base_qualities,
                locus_start=locus_start_in_read,
                locus_end=locus_end_in_read,
                locus_size=n_bases,
                is_deletion=is_deletion,
                name=name)
