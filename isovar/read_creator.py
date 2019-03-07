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

from __future__ import print_function, division, absolute_import

from six import integer_types

from .allele_read import AlleleRead
from .allele_read_helpers import filter_non_alt_reads_for_variant
from .default_parameters import (
    USE_SECONDARY_ALIGNMENTS,
    USE_DUPLICATE_READS,
    MIN_READ_MAPPING_QUALITY,
    USE_SOFT_CLIPPED_BASES,
)
from .locus_read import LocusRead
from .logging import get_logger
from .string_helpers import convert_from_bytes_if_necessary, trim_N_nucleotides
from .variant_helpers import trim_variant

logger = get_logger(__name__)


class ReadCreator(object):
    """
    LocusReadCreator holds options related to extracting reads from SAM/BAM alignment files
    and provides methods for different ways to create LocusRead objects.
    """
    def __init__(
            self,
            use_secondary_alignments=USE_SECONDARY_ALIGNMENTS,
            use_duplicate_reads=USE_DUPLICATE_READS,
            min_mapping_quality=MIN_READ_MAPPING_QUALITY,
            use_soft_clipped_bases=USE_SOFT_CLIPPED_BASES):
        """
        Parameters
        ----------
        use_secondary_alignments : bool
            Use a read even when it's not the primary alignment at a locus

        use_duplicate_reads : bool
            Use a read even if it's been marked as a duplicate

        min_mapping_quality : int
            Minimum MAPQ (mapping quality) to use a read

        use_soft_clipped_bases : bool
            Include soft-clipped positions on a read which were ignored by the aligner
        """
        self.use_secondary_alignments = use_secondary_alignments
        self.use_duplicate_reads = use_duplicate_reads
        self.min_mapping_quality = min_mapping_quality
        self.use_soft_clipped_bases = use_soft_clipped_bases

    def locus_read_from_pysam_aligned_segment(
            self,
            pysam_aligned_segment,
            base0_start_inclusive,
            base0_end_exclusive):
        """
        Create LocusRead from pysam.AlignedSegment object and the start/end indices
        of the locus of interest. If any essential fields of the aligned segment
        are missing then None is returned instead.

        Parameters
        ----------
        pysam_aligned_segment : pysam.AlignedSegment
            AlignedSegment is expected to overlap the locus

        base0_start_inclusive : int

        base0_end_exclusive : int

        Returns
        -------
        LocusRead or None
        """
        if not isinstance(base0_start_inclusive, integer_types):
            raise TypeError("Expected base0_start_inclusive to be an integer but got %s" % (
                type(base0_start_inclusive),))
        if not isinstance(base0_end_exclusive, integer_types):
            raise TypeError("Expected base0_end_exclusive to be an integer but got %s" % (
                type(base0_end_exclusive),))

        name = pysam_aligned_segment.query_name
        if name is None:
            logger.warn(
                "Read missing name at position %d",
                base0_start_inclusive + 1)
            return None

        if pysam_aligned_segment.is_unmapped:
            logger.warn(
                "How did we get unmapped read '%s' in a pileup?", name)
            return None

        if pysam_aligned_segment.is_secondary and not self.use_secondary_alignments:
            logger.debug("Skipping secondary alignment of read '%s'", name)
            return None

        if pysam_aligned_segment.is_duplicate and not self.use_duplicate_reads:
            logger.debug("Skipping duplicate read '%s'", name)
            return None

        mapping_quality = pysam_aligned_segment.mapping_quality

        if self.min_mapping_quality > 0 and (mapping_quality is None):
            logger.debug("Skipping read '%s' due to missing MAPQ" % name)
            return None
        elif mapping_quality < self.min_mapping_quality:
            logger.debug(
                "Skipping read '%s' due to low MAPQ: %d < %d",
                name,
                mapping_quality,
                self.min_mapping_quality)
            return None

        sequence = pysam_aligned_segment.query_sequence
        if sequence is None:
            logger.warn("Skipping read '%s' due to missing sequence" % name)
            return None

        base_qualities = pysam_aligned_segment.query_qualities

        if base_qualities is None:
            logger.warn("Skipping read '%s' due to missing base qualities" % name)
            return None
        elif len(base_qualities) != len(sequence):
            logger.warn(
                "Skipping read '%s' due to mismatch in length of sequence (%d) and qualities (%d)" % (
                    name,
                    len(sequence),
                    len(base_qualities)))
            return None
        # By default, AlignedSegment.get_reference_positions only returns base-1 positions
        # from the reference that are within the alignment. If full_length is set,
        # None values will be included for any soft-clipped or unaligned positions
        # within the read. The returned list will thus be of the same
        # length as the read.

        base0_reference_positions = pysam_aligned_segment.get_reference_positions(full_length=True)

        if len(base0_reference_positions) != len(base_qualities):
            logger.warn(
                "Skipping read '%s' due to mismatch in length of positions (%d) and qualities (%d)" % (
                    name,
                    len(base0_reference_positions),
                    len(base_qualities)))
            return None

        base0_reference_positions_dict = {
            base0_reference_pos: base0_read_pos
            for (base0_read_pos, base0_reference_pos)
            in enumerate(base0_reference_positions)
            if base0_reference_pos is not None
        }

        reference_interval_size = base0_end_exclusive - base0_start_inclusive
        if reference_interval_size < 0:
            raise ValueError("Unexpected interval start after interval end")

        # TODO: Consider how to handle variants before splice sites, where
        # the bases before or after on the genome will not be mapped on the
        # read
        #
        # we have a dictionary mapping base-1 reference positions to base-0
        # read indices and we need to use that to convert the reference
        # half-open interval into a half-open interval on the read.
        if reference_interval_size == 0:
            # Reference interval is between two bases but read may contain
            # insertion.
            #
            # Reference:
            #   Insertion location:       *
            #   Reference position: 10000 | 10001 10002 10003 10004 10005 10006 10007
            #   Base sequence:        A   |   T     G     C     A     A     A     A
            #
            # Read with inserted nucleotide:
            #   Read position:      00000 00001 00002 00003 00004 00005 00006 00007
            #   Base sequence:        A    *A*    T     G     C     A     A     A
            #   Reference position: 10000 ----- 10001 10002 10003 10004 10005 10006
            #
            # The start/end of the reference interval may be mapped to a read position,
            # in this case reference:10000 -> read:00000, but it would be incorrect
            # to take this position as the start/end of the insertion on the read
            # since it does not cover the inserted bases. Instead, we look at the
            # read position of the next base in the reference and, if it's more than
            # 1 base away from the start, use that as the end of the interval. If it's
            # next to the start of the interval then we return the empty "between bases"
            # interval of [start, start).
            #
            # To deal with insertions at the beginning and end of a read we're
            # going to allow the start/end to be None.
            reference_position_before_insertion = base0_start_inclusive - 1
            reference_position_after_insertion = base0_start_inclusive
            read_base0_before_insertion = base0_reference_positions_dict.get(
                reference_position_before_insertion)
            read_base0_after_insertion = base0_reference_positions_dict.get(
                reference_position_after_insertion)
            print("[%d:%d)" % (base0_start_inclusive, base0_end_exclusive))
            print(sequence)
            print(base0_reference_positions)
            print("ref", reference_position_before_insertion, reference_position_after_insertion)

            print("read", read_base0_before_insertion, read_base0_after_insertion)

            if read_base0_before_insertion is None:
                logger.warning("Cannot use read '%s' because reference position %d is not mapped" % (
                    name,
                    reference_position_before_insertion))
                return None
            elif read_base0_after_insertion is None:
                logger.warning("Cannot use read '%s' because reference position %d is not mapped" % (
                    name,
                    reference_position_after_insertion))
                return None
            elif read_base0_after_insertion - read_base0_after_insertion == 1:
                read_base0_start_inclusive = read_base0_end_exclusive = read_base0_before_insertion + 1
            else:
                read_base0_start_inclusive = read_base0_before_insertion + 1
                read_base0_end_exclusive = read_base0_after_insertion
        else:
            # Reference bases are selected for match or deletion.
            #
            # What happens if the reference bases are interspersed with insertions?
            # Reference:
            #   10000 10001 10002 10003 10004 10005 10006 10007
            #     A     T     G     C     A     A     A     A
            #
            # Read:
            #   00000 00001 00002 00003 00004 00005 00006 00007
            #     A    *A*     T     G     C     A     A     A
            #
            # ...and our reference interval is base-1 inclusive 10000:10001
            # but the read has an inserted 'A' in between the two bases.
            #
            # In this case we need to figure out the first and last positions
            # which match the inclusive interval and then convert it to a half-open
            # interval. One possibly more obvious alternative is just to
            # figure out which read indices correspond to base0_start_inclusive and
            # base0_end_exclusive but this would fail if base0_end_exclusive is
            # after the end the end of the read.
            read_base0_start_inclusive = base0_reference_positions_dict.get(base0_start_inclusive)
            if read_base0_start_inclusive is None:
                # if first base of reference locus isn't mapped, try getting the base
                # before it and then adding one to its corresponding base index
                reference_base0_position_before_locus = base0_start_inclusive - 1
                if reference_base0_position_before_locus in base0_reference_positions_dict:
                    read_base0_position_before_locus = base0_reference_positions_dict[
                        reference_base0_position_before_locus]
                    read_base0_start_inclusive = read_base0_position_before_locus + 1
                else:
                    logger.warning(
                        "Cannot use read '%s' because neither reference positions %d or %d are not mapped" % (
                            name,
                            base0_start_inclusive,
                            reference_base0_position_before_locus))
                    return None

            read_base0_end_exclusive = base0_reference_positions_dict.get(base0_end_exclusive)
            if read_base0_end_exclusive is None:
                # if exclusve last index of reference interval doesn't have a corresponding
                # base position then try getting the base position of the reference
                # position before it and then adding one
                reference_base0_end_inclusive = base0_end_exclusive - 1
                if reference_base0_end_inclusive in base0_reference_positions_dict:
                    read_base0_end_inclusive = base0_reference_positions_dict[
                        reference_base0_end_inclusive]
                    read_base0_end_exclusive = read_base0_end_inclusive + 1
                else:
                    logger.warning(
                        "Cannot use read '%s' because neither reference positions %d or %d are not mapped" % (
                            name,
                            base0_end_exclusive,
                            reference_base0_end_inclusive))
                    return None


        if isinstance(sequence, bytes):
            sequence = sequence.decode('ascii')

        if not self.use_soft_clipped_bases:
            # if we're not allowing soft clipped based then
            # the fraction of the read which is usable may be smaller
            # than the sequence, qualities, and alignment positions
            # we've extracted, so slice through those to get rid of
            # soft-clipped ends of the read
            aligned_subsequence_start = pysam_aligned_segment.query_alignment_start
            aligned_subsequence_end = pysam_aligned_segment.query_alignment_end
            sequence = sequence[aligned_subsequence_start:aligned_subsequence_end]
            base0_reference_positions = base0_reference_positions[
                aligned_subsequence_start:aligned_subsequence_end]
            base_qualities = base_qualities[aligned_subsequence_start:aligned_subsequence_end]
            if read_base0_start_inclusive is not None:
                read_base0_start_inclusive -= aligned_subsequence_start
            if read_base0_end_exclusive is not None:
                read_base0_end_exclusive -= aligned_subsequence_start
        return LocusRead(
            name=name,
            sequence=sequence,
            reference_positions=base0_reference_positions,
            quality_scores=base_qualities,
            reference_base0_start_inclusive=base0_start_inclusive,
            reference_base0_end_exclusive=base0_end_exclusive,
            read_base0_start_inclusive=read_base0_start_inclusive,
            read_base0_end_exclusive=read_base0_end_exclusive)

    def get_locus_reads(
            self,
            alignments,
            chromosome,
            base0_start_inclusive,
            base0_end_exclusive):
        """
        Create LocusRead objects for reads which overlap the given chromosome,
        start, and end positions. The actual work to figure out if what's between
        those positions matches a variant happens later when LocusRead objects are
        converted to AlleleRead objects.

        Parameters
        ----------
        alignments : pysam.AlignmentFile

        chromosome : str

        base0_start_inclusive : int
            Start of genomic interval, base 0 and inclusive

        base0_end_exclusive : int
            End of genomic interval, base 0 and exclusive

        Returns a sequence of ReadAtLocus objects
        """
        logger.debug(
            "Gathering reads at locus %s:%d-%d",
            chromosome,
            base0_start_inclusive,
            base0_end_exclusive)
        reads = []
        for aligned_segment in alignments.fetch(
                contig=chromosome,
                start=base0_start_inclusive,
                stop=base0_end_exclusive):
            read = self.locus_read_from_pysam_aligned_segment(
                aligned_segment,
                base0_start_inclusive=base0_start_inclusive,
                base0_end_exclusive=base0_end_exclusive)
            if read is not None:
                reads.append(read)
        logger.info(
            "Found %d reads overlapping locus %s:%d-%d",
            len(reads),
            chromosome,
            base0_start_inclusive,
            base0_end_exclusive)
        return reads

    def allele_reads_overlapping_variant(
            self,
            alignments,
            variant,
            chromosome=None):
        """
        Find reads in the given SAM/BAM file which overlap the given variant and
        return them as a list of AlleleRead objects.

        Parameters
        ----------
        alignments : pysam.AlignmentFile

        variant : varcode.Variant

        chromosome : str or None

        Returns sequence of AlleleRead objects.
        """
        if chromosome is None:
            chromosome = variant.contig

        logger.info(
            "Gathering variant reads for variant %s (with gene names %s)",
            variant,
            variant.gene_names)

        base1_position, ref, alt = trim_variant(variant)

        if len(ref) == 0:
            # If there is no reference sequence in the variant
            # then it's an insertion and the base 0 coordinates
            # will select the space between two bases.
            #
            # For example, an insertion between base-1 positions chr1:3 and chr1:4
            #
            # Base 1 inclusive:
            #   |   1   |   2   |   3   |   4   |   5   |
            # Base 0 half-open:
            #   0       1       2       3       4       5
            #
            #  The space between chr1:3 and chr1:4 in base-0 coordinates is chr1 3:3
            #
            # So, to convert an insertion from base-1 inclusive to base-0 half-open we
            # keep the original position
            base0_start_inclusive = base1_position
            base0_end_exclusive = base1_position
        else:
            # if variant is SNV or deletion then some number of reference bases
            # are selected, so just get the interval for those.
            #
            # For example, if two bases at positions chr1:1000 and 1001 are deleted
            # then the base0 indices will be 9999:1001
            base0_start_inclusive = base1_position - 1
            base0_end_exclusive = base0_start_inclusive + len(ref)

        locus_reads = self.get_locus_reads(
            alignments=alignments,
            chromosome=chromosome,
            base0_start_inclusive=base0_start_inclusive,
            base0_end_exclusive=base0_end_exclusive)

        allele_reads = []

        for locus_read in locus_reads:
            allele_read = self.allele_read_from_locus_read(locus_read)
            if allele_read is None:
                continue
            else:
                allele_reads.append(allele_read)
        return allele_reads

    def allele_reads_overlapping_variants(self, variants, alignments):
        """
        Generates sequence of tuples, each containing a variant paired with
        a list of AlleleRead objects.

        Parameters
        ----------
        variants : varcode.VariantCollection

        alignments : pysam.AlignmentFile
        """
        chromosome_names = set(alignments.references)
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
                logger.warning(
                    "Chromosome '%s' from variant %s not in alignment file %s",
                    variant.contig,
                    variant,
                    alignments.filename)
                yield variant, []
                continue
            allele_reads = self.allele_reads_overlapping_variant(
                alignments=alignments,
                chromosome=chromosome,
                variant=variant)
            yield variant, allele_reads


    def allele_reads_supporting_variant(self, variant, alignments):
        """
        Parameters
        ----------
        variant: varcode.Variant
    
        alignments:  pysam.AlignmentFile
    
        Given a variant and a SAM/BAM file, finds all AlleleRead objects overlapping
        a variant and returns those that support the variant's alt allele.
        """
        allele_reads = self.allele_reads_overlapping_variant(
            variant=variant,
            alignments=alignments)
        return filter_non_alt_reads_for_variant(
            variant=variant,
            allele_reads=allele_reads)
    
    
    def allele_reads_supporting_variants(self, variants, alignments):
        """
        Given a SAM/BAM file and a collection of variants, generates a sequence
        of variants paired with reads which support each variant.

        Parameters
        ----------
        variants : varcode.VariantCollection

        alignments : pysam.AlignmentFile

        Generator of (varcode.Variant, list of AlleleRead) pairs.
        """
        for variant, allele_reads in self.reads_overlapping_variants(
                variants=variants,
                alignments=alignments):
            yield variant, filter_non_alt_reads_for_variant(variant, allele_reads)

    def allele_read_from_locus_read(self, locus_read):
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
                "Skipping read '%s' because some required bases in reference interval %s:%s aren't mapped",
                read_name,
                reference_base0_start_inclusive,
                reference_base0_end_exclusive)
            return None

        reference_positions = locus_read.reference_positions

        n_ref_bases = reference_base0_end_exclusive - reference_base0_start_inclusive

        insertion = (n_ref_bases == 0)

        if insertion:
            # insertions require a sequence of non-aligned bases
            # followed by the subsequence reference position
            for read_index in range(read_base0_start_inclusive, read_base0_end_exclusive):
                # all the inserted nucleotides should *not* align to the reference
                if reference_positions[read_index] is not None:
                    logger.debug(
                        "Skipping read '%s', inserted nucleotides shouldn't map to reference",
                        read_name)
                    return None

        nucleotides_at_variant_locus = convert_from_bytes_if_necessary(
            sequence[read_base0_start_inclusive:read_base0_end_exclusive])

        if "N" in nucleotides_at_variant_locus:
            logger.debug(
                "Skipping read '%s', found N nucleotides at variant locus",
                read_name)
        prefix = convert_from_bytes_if_necessary(sequence[:read_base0_start_inclusive])
        suffix = convert_from_bytes_if_necessary(sequence[read_base0_end_exclusive:])

        prefix, suffix = trim_N_nucleotides(prefix, suffix)

        return AlleleRead(
            prefix,
            nucleotides_at_variant_locus,
            suffix,
            name=read_name)
