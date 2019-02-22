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

from .allele_read import AlleleRead
from .allele_read_helpers import filter_non_alt_reads_for_variant
from .default_parameters import (
    USE_SECONDARY_ALIGNMENTS,
    USE_DUPLICATE_READS,
    MIN_READ_MAPPING_QUALITY
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
            use_soft_clipped_bases=False):
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

        # By default, AlignedSegment.get_reference_positions only returns base-1 positions
        # from the reference that are within the alignment. If full_length is set,
        # None values will be included for any soft-clipped or unaligned positions
        # within the read. The returned list will thus be of the same
        # length as the read.
        base1_reference_positions = pysam_aligned_segment.get_reference_positions(full_length=True)
        reference_positions_dict = {
            base1_reference_pos: base0_read_pos
            for (base0_read_pos, base1_reference_pos)
            in enumerate(base1_reference_positions)
            if base1_reference_pos is not None
        }

        reference_interval_size = base0_end_exclusive - base0_start_inclusive
        if reference_interval_size < 0:
            raise ValueError("Unexpected interval start after interval end")

        # we have a dictionary mapping base-1 reference positions to base-0
        # read indices and we need to use that to convert the reference
        # half-open interval into a half-open interval on the read.
        if reference_interval_size == 0:
            # To deal with insertions at the beginning and end of a read we're
            # # going to try two approaches to figure out the read interval
            # corresponding with a reference interval.
            #
            # First, try getting the read position of the base before the insertion
            # and if that's not mapped, we'll try getting the position of the
            # base after the insertion.
            position_before_insertion = base0_start_inclusive
            position_after_insertion = base0_start_inclusive + 1

            if position_before_insertion in reference_positions_dict:
                locus_on_read_base0_start = reference_positions_dict.get(position_before_insertion)
            elif position_after_insertion in reference_positions_dict:
                locus_on_read_base0_start = reference_positions_dict.get(position_after_insertion) - 1
            else:
                locus_on_read_base0_start = None
            locus_on_read_base0_end = locus_on_read_base0_start
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
            base1_start_inclusive = base0_start_inclusive + 1
            base1_end_inclusive = base1_start_inclusive + reference_interval_size
            locus_on_read_base0_start = reference_positions_dict.get(base1_start_inclusive)
            locus_on_read_base0_end_inclusive_index = reference_positions_dict.get(base1_end_inclusive)
            locus_on_read_base0_end = locus_on_read_base0_end_inclusive_index + 1

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
            base1_reference_positions = base1_reference_positions[
                aligned_subsequence_start:aligned_subsequence_end]
            base_qualities = base_qualities[aligned_subsequence_start:aligned_subsequence_end]
            if locus_on_read_base0_start is not None:
                locus_on_read_base0_start -= aligned_subsequence_start
            if locus_on_read_base0_end is not None:
                locus_on_read_base0_end -= aligned_subsequence_start

        return LocusRead(
            name=name,
            sequence=sequence,
            reference_positions=base1_reference_positions,
            quality_scores=base_qualities,
            reference_base0_start_inclusive=base0_start_inclusive,
            reference_base0_end_exclusive=base0_end_exclusive,
            read_base0_start_inclusive=locus_on_read_base0_start,
            read_base0_end_exclusive=locus_on_read_base0_end)

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
                base0_position_before_variant=base0_start_inclusive,
                base0_position_after_variant=base0_end_exclusive)
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
            allele_read = self.allele_read_from_locus_read(
                locus_read=locus_read,
                n_ref=len(ref),
                n_alt=len(alt))
            if allele_read is None:
                continue
            else:
                allele_reads.append(allele_read)
        return allele_reads

    def reads_overlapping_variants(self, variants, alignments):
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


    def reads_supporting_variant(self, variant, alignments):
        """
        Parameters
        ----------
        variant: varcode.Variant
    
        alignments:  pysam.AlignmentFile
    
        Given a variant and a SAM/BAM file, finds all AlleleRead objects overlapping
        a variant and returns those that support the variant's alt allele.
        """
        allele_reads = self.reads_overlapping_variant(
            variant=variant,
            alignments=alignments)
        return filter_non_alt_reads_for_variant(
            variant=variant,
            allele_reads=allele_reads)
    
    
    def reads_supporting_variants(self, variants, alignments):
        """
        Given a SAM/BAM file and a collection of variants, generates a sequence
        of variants paired with reads which support each variant.
        """
        for variant, allele_reads in self.reads_overlapping_variants(
                variants=variants,
                alignments=alignments):
            yield variant, filter_non_alt_reads_for_variant(variant, allele_reads)

    def allele_read_from_locus_read(self, locus_read, n_ref, n_alt):
        """
        Given a single LocusRead object, return either an AlleleRead or None

        Parameters
        ----------
        locus_read : LocusRead
            Read which overlaps a variant locus but doesn't necessarily contain the
            alternate nucleotides

        n_ref : int
            Number of reference positions we are expecting to be modified or
            deleted (for insertions this should be 0)

        n_alt : int
        """
        sequence = locus_read.sequence
        reference_positions = locus_read.reference_positions
        reference_base0_start_inclusive = locus_read.reference_base0_start_inclusive
        reference_base0_end_exclusive = locus_read.reference_base0_end_exclusive
        read_base0_start_inclusive = locus_read.read_base0_start_inclusive
        read_base0_end_exclusive = locus_read.read_base0_end_exclusive

        insertion = (n_ref == 0)
        deletion = (n_alt == 0)
        if insertion:
            if ref_pos_after - ref_pos_before != 1:
                # if the number of nucleotides skipped isn't the same
                # as the number of reference nucleotides in the variant then
                # don't use this read
                logger.debug(
                    "Positions before (%d) and after (%d) variant should be adjacent on read %s",
                    ref_pos_before,
                    ref_pos_after,
                    locus_read)
                return None

            # insertions require a sequence of non-aligned bases
            # followed by the subsequence reference position
            ref_positions_for_inserted = reference_positions[
                read_pos_before + 1:read_pos_after]
            if any(insert_pos is not None for insert_pos in ref_positions_for_inserted):
                # all these inserted nucleotides should *not* align to the
                # reference
                logger.debug(
                    "Skipping read, inserted nucleotides shouldn't map to reference")
                return None
        else:
            # substitutions and deletions
            if ref_pos_after - ref_pos_before != n_ref + 1:
                # if the number of nucleotides skipped isn't the same
                # as the number of reference nucleotides in the variant then
                # don't use this read
                logger.debug(
                    ("Positions before (%d) and after (%d) variant should be "
                     "adjacent on read %s"),
                    ref_pos_before,
                    ref_pos_after,
                    locus_read)
                return None

        nucleotides_at_variant_locus = sequence[read_pos_before + 1:read_pos_after]

        prefix = sequence[:read_pos_before + 1]
        suffix = sequence[read_pos_after:]

        prefix, suffix = convert_from_bytes_if_necessary(prefix, suffix)
        prefix, suffix = trim_N_nucleotides(prefix, suffix)

        return AlleleRead(
            prefix,
            nucleotides_at_variant_locus,
            suffix,
            name=locus_read.name)
