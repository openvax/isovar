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

from .default_parameters import (
    USE_SECONDARY_ALIGNMENTS,
    USE_DUPLICATE_READS,
    MIN_READ_MAPPING_QUALITY,
    USE_SOFT_CLIPPED_BASES,
)
from .locus_read import LocusRead
from .logging import get_logger
from .allele_read_helpers import allele_reads_from_locus_reads
from .variant_helpers import trim_variant
from .read_evidence import ReadEvidence

logger = get_logger(__name__)


class ReadCollector(object):
    """
    ReadCollector holds options related to extracting reads from SAM/BAM alignment files
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
        if pysam_aligned_segment.is_secondary and not self.use_secondary_alignments:
            return None

        if pysam_aligned_segment.is_duplicate and not self.use_duplicate_reads:
            return None

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
        base0_reference_positions = \
            pysam_aligned_segment.get_reference_positions(full_length=True)

        if len(base0_reference_positions) != len(base_qualities):
            logger.warn(
                "Skipping read '%s' due to mismatch in length of positions (%d) and qualities (%d)" % (
                    name,
                    len(base0_reference_positions),
                    len(base_qualities)))
            return None

        reference_interval_size = base0_end_exclusive - base0_start_inclusive
        if reference_interval_size < 0:
            raise ValueError("Unexpected interval start after interval end")

        # TODO:
        #  Consider how to handle variants before splice sites, where
        #  the bases before or after on the genome will not be mapped on the
        #  read
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
            if reference_position_before_insertion in base0_reference_positions:
                read_base0_before_insertion = \
                    base0_reference_positions.index(
                        reference_position_before_insertion)
            else:
                return None

            if reference_position_after_insertion in base0_reference_positions:
                read_base0_after_insertion = base0_reference_positions.index(
                    reference_position_after_insertion)
            else:
                return None

            if read_base0_after_insertion - read_base0_after_insertion == 1:
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
            if base0_start_inclusive in base0_reference_positions:
                read_base0_start_inclusive = base0_reference_positions.index(base0_start_inclusive)
            elif base0_start_inclusive - 1 in base0_reference_positions:
                # if first base of reference locus isn't mapped, try getting the base
                # before it and then adding one to its corresponding base index
                read_base0_position_before_locus = \
                    base0_reference_positions.index(base0_start_inclusive - 1)
                read_base0_start_inclusive = read_base0_position_before_locus + 1
            else:
                return None

            if base0_end_exclusive in base0_reference_positions:
                read_base0_end_exclusive = \
                    base0_reference_positions.index(base0_end_exclusive)
            elif (base0_end_exclusive - 1) in base0_reference_positions:
                # if exclusive last index of reference interval doesn't have a corresponding
                # base position then try getting the base position of the reference
                # position before it and then adding one
                read_base0_end_inclusive = \
                    base0_reference_positions.index(base0_end_exclusive - 1)
                read_base0_end_exclusive = read_base0_end_inclusive + 1
            else:
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
            alignment_file,
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
        alignment_file : pysam.AlignmentFile

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
        total_count = 0
        # check overlap against wider overlap to make sure we don't miss
        # any reads
        base0_pos_before_start = base0_start_inclusive - 1
        base0_pos_after_end = base0_end_exclusive + 1
        for aligned_segment in alignment_file.fetch(
                chromosome,
                base0_start_inclusive,
                base0_end_exclusive):
            total_count += 1
            # we get a significant speed up if we skip reads that have spliced
            # out the entire interval of interest. this is redundant with the
            # attempt to find mapping positions in
            #   self.locus_read_from_pysam_aligned_segment
            # but we do it here to skip the function call overhead for loci
            # where ~1M reads are mapped
            if aligned_segment.get_overlap(base0_pos_before_start, base0_pos_after_end) == 0:
                continue
            read = self.locus_read_from_pysam_aligned_segment(
                aligned_segment,
                base0_start_inclusive=base0_start_inclusive,
                base0_end_exclusive=base0_end_exclusive)
            if read is not None:
                reads.append(read)
        logger.info(
            "Kept %d/%d reads overlapping locus %s:%d-%d",
            len(reads),
            total_count,
            chromosome,
            base0_start_inclusive,
            base0_end_exclusive)
        return reads

    @staticmethod
    def _infer_chromosome_name(variant_chromosome_name, valid_chromosome_names):
        """
        In case the variant is using an hg19 reference name and the alignment
        was against b37 (or vice versa) we have to check whether adding or removing
        the prefix "chr" is necessary.
        Parameters
        ----------
        variant_chromosome_name : str

        valid_chromosome_names : set of str

        Returns
        -------
        str or None
        """
        # I imagine the conversation went like this:
        # A: "Hey, I have an awesome idea"
        # B: "What's up?"
        # A: "Let's make two nearly identical reference genomes"
        # B: "But...that sounds like it might confuse people."
        # A: "Nah, it's cool, we'll give the chromosomes different prefixes!"
        # B: "OK, sounds like a good idea."
        candidate_names = {variant_chromosome_name}
        if variant_chromosome_name.startswith("chr"):
            candidate_names.add(variant_chromosome_name[3:])
        else:
            candidate_names.add("chr" + variant_chromosome_name)
        for candidate in list(candidate_names):
            candidate_names.add(candidate.lower())
            candidate_names.add(candidate.upper())
        for candidate in candidate_names:
            if candidate in valid_chromosome_names:
                return candidate
        return None

    def locus_reads_overlapping_variant(
            self,
            alignment_file,
            variant,
            chromosome=None):
        """
        Find reads in the given SAM/BAM file which overlap the given variant and
        return them as a list of LocusRead objects.

        Parameters
        ----------
        alignment_file : pysam.AlignmentFile

        variant : varcode.Variant

        chromosome : str or None

        Returns sequence of LocusRead objects.
        """
        if chromosome is None:
            # if a chromosome name isn't manually specified then try
            # to figure out whether adding or removing "chr" is necessary
            # match chromosome names used for variant calling and those
            # found in read alignments
            chromosome = self._infer_chromosome_name(
                variant_chromosome_name=variant.contig,
                valid_chromosome_names=set(alignment_file.references))

        if chromosome is None:
            # failed to infer a chromsome name for this variant which
            # matches names used in SAM/BAM file
            logger.warning(
                "Chromosome '%s' from variant %s not in alignment file %s",
                variant.contig,
                variant,
                alignment_file.filename)
            return []

        logger.info(
            "Gathering reads for variant %s (with gene names %s)",
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

        return self.get_locus_reads(
            alignment_file=alignment_file,
            chromosome=chromosome,
            base0_start_inclusive=base0_start_inclusive,
            base0_end_exclusive=base0_end_exclusive)

    def allele_reads_overlapping_variant(
            self,
            variant,
            alignment_file):
        """
        Find reads in the given SAM/BAM file which overlap the given variant and
        return them as a list of AlleleRead objects.

        Parameters
        ----------
        variant : varcode.Variant

        alignment_file : pysam.AlignmentFile
            Aligned RNA reads

        Returns sequence of AlleleRead objects.
        """
        return allele_reads_from_locus_reads(
            self.locus_reads_overlapping_variant(
                alignment_file=alignment_file,
                variant=variant))

    def read_evidence_for_variant(
            self,
            variant,
            alignment_file):
        """
        Find reads in the given SAM/BAM file which overlap the given variant and
        return them as a ReadEvidence object, which splits the reads into
        ref/alt/other groups.

        Parameters
        ----------
        variant : varcode.Variant

        alignment_file : pysam.AlignmentFile
            Aligned RNA reads

        Returns ReadEvidence
        """
        allele_reads = self.allele_reads_overlapping_variant(
            variant=variant,
            alignment_file=alignment_file)
        return ReadEvidence.from_variant_and_allele_reads(
            variant,
            allele_reads)

    def allele_reads_supporting_variant(self, variant, alignment_file):
        """
        Gather AlleleRead objects which contain the same allele as the variant.

        Parameters
        ----------
        variant : varcode.VariantCollection
            Variants which will be the keys of the result

        alignment_file : pysam.AlignmentFile
            Aligned RNA reads

        Returns list of AlleleRead
        """
        read_evidence = self.read_evidence_for_variant(
            variant=variant,
            alignment_file=alignment_file)
        return read_evidence.alt_reads

    def read_evidence_generator(self, variants, alignment_file):
        """
        Consumes a generator of varcode.Variant objects, collects read evidence
        for each variant from the alignment_file, and generates a sequence
        of (Variant, ReadEvidence) pairs.

        Parameters
        ----------
        variants : varcode.VariantCollection
            Variants which will be the keys of the result

        alignment_file : pysam.AlignmentFile
            Aligned RNA reads

        Generates sequence of (varcode.Variant, ReadEvidence) pairs
        """
        for variant in variants:
            read_evidence = self.read_evidence_for_variant(
                variant=variant,
                alignment_file=alignment_file)
            yield variant, read_evidence
