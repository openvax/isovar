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

from collections import defaultdict

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
        use_soft_clipped_bases=USE_SOFT_CLIPPED_BASES,
        merge_overlapping_fragments=False,
    ):
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

        merge_overlapping_fragments : bool
            Merge overlapping paired-end reads from the same fragment into one
            fragment-level sequence while preserving the raw alignment count in
            `source_read_count`.
        """
        self.use_secondary_alignments = use_secondary_alignments
        self.use_duplicate_reads = use_duplicate_reads
        self.min_mapping_quality = min_mapping_quality
        self.use_soft_clipped_bases = use_soft_clipped_bases
        self.merge_overlapping_fragments = merge_overlapping_fragments

    @staticmethod
    def _previous_fully_aligned_pair_index(aligned_pairs, start_index):
        """
        Find the previous aligned-pair entry which maps both a query and
        reference position.
        """
        while start_index >= 0:
            query_pos, ref_pos, _ = aligned_pairs[start_index]
            if query_pos is not None and ref_pos is not None:
                return start_index
            start_index -= 1
        return None

    @classmethod
    def _iter_left_aligned_indel_events(
        cls,
        aligned_pairs,
        base1_start,
        ref,
        alt,
        read_start,
        read_end,
        previous_aligned_index,
    ):
        """
        Yield the queried indel and each successive one-base left shift within
        the local alignment.

        The query interval is shifted in lockstep with the indel representation
        so that downstream AlleleRead objects see the same canonical split for
        every equivalent alignment of the same indel.
        """
        ref = ref.upper()
        alt = alt.upper()

        while True:
            yield base1_start, ref, alt, read_start, read_end

            if previous_aligned_index is None or base1_start <= 1:
                break

            _, _, previous_ref_base = aligned_pairs[previous_aligned_index]
            if previous_ref_base is None:
                break

            previous_ref_base = previous_ref_base.upper()
            longer_allele = alt if len(alt) > len(ref) else ref
            if len(longer_allele) == 0 or previous_ref_base != longer_allele[-1]:
                break

            if len(alt) > len(ref):
                alt = previous_ref_base + alt[:-1]
            else:
                ref = previous_ref_base + ref[:-1]

            base1_start -= 1
            read_start -= 1
            read_end -= 1
            previous_aligned_index = cls._previous_fully_aligned_pair_index(
                aligned_pairs,
                previous_aligned_index - 1,
            )

    @classmethod
    def _left_aligned_indel_interval_for_variant(
        cls,
        pysam_aligned_segment,
        trimmed_base1_start,
        trimmed_ref,
        trimmed_alt,
    ):
        """
        If the read contains an equivalent indel aligned to the right of the
        queried variant locus, return the canonical read interval for the
        left-aligned representation of that indel.
        """
        is_insertion = len(trimmed_ref) == 0 and len(trimmed_alt) > 0
        is_deletion = len(trimmed_alt) == 0 and len(trimmed_ref) > 0
        if not (is_insertion or is_deletion):
            return None

        try:
            aligned_pairs = pysam_aligned_segment.get_aligned_pairs(
                matches_only=False,
                with_seq=True,
            )
        except ValueError:
            logger.debug(
                "Skipping indel left-alignment for read '%s' because the alignment "
                "does not expose reference bases (for example, missing MD tag)",
                pysam_aligned_segment.query_name,
            )
            return None

        query_sequence = pysam_aligned_segment.query_sequence
        if query_sequence is None:
            return None
        if isinstance(query_sequence, bytes):
            query_sequence = query_sequence.decode("ascii")

        i = 0
        while i < len(aligned_pairs):
            query_pos, ref_pos, ref_base = aligned_pairs[i]

            if is_insertion and query_pos is not None and ref_pos is None:
                j = i
                while (
                    j < len(aligned_pairs)
                    and aligned_pairs[j][0] is not None
                    and aligned_pairs[j][1] is None
                ):
                    j += 1

                previous_aligned_index = cls._previous_fully_aligned_pair_index(
                    aligned_pairs,
                    i - 1,
                )
                if previous_aligned_index is not None:
                    previous_query_pos, previous_ref_pos, _ = aligned_pairs[
                        previous_aligned_index
                    ]
                    inserted_sequence = query_sequence[query_pos:aligned_pairs[j - 1][0] + 1]
                    for (
                        normalized_base1_start,
                        normalized_ref,
                        normalized_alt,
                        normalized_read_start,
                        normalized_read_end,
                    ) in cls._iter_left_aligned_indel_events(
                        aligned_pairs=aligned_pairs,
                        base1_start=previous_ref_pos + 1,
                        ref="",
                        alt=inserted_sequence,
                        read_start=query_pos,
                        read_end=aligned_pairs[j - 1][0] + 1,
                        previous_aligned_index=previous_aligned_index,
                    ):
                        if (
                            normalized_base1_start == trimmed_base1_start
                            and normalized_ref == trimmed_ref
                            and normalized_alt == trimmed_alt.upper()
                        ):
                            return normalized_read_start, normalized_read_end
                i = j
                continue

            if is_deletion and query_pos is None and ref_pos is not None:
                j = i
                while (
                    j < len(aligned_pairs)
                    and aligned_pairs[j][0] is None
                    and aligned_pairs[j][1] is not None
                ):
                    j += 1

                previous_aligned_index = cls._previous_fully_aligned_pair_index(
                    aligned_pairs,
                    i - 1,
                )
                if previous_aligned_index is not None:
                    previous_query_pos, _, _ = aligned_pairs[previous_aligned_index]
                    deleted_sequence = "".join(
                        aligned_pairs[k][2].upper() for k in range(i, j)
                    )
                    for (
                        normalized_base1_start,
                        normalized_ref,
                        normalized_alt,
                        normalized_read_start,
                        normalized_read_end,
                    ) in cls._iter_left_aligned_indel_events(
                        aligned_pairs=aligned_pairs,
                        base1_start=ref_pos + 1,
                        ref=deleted_sequence,
                        alt="",
                        read_start=previous_query_pos + 1,
                        read_end=previous_query_pos + 1,
                        previous_aligned_index=previous_aligned_index,
                    ):
                        if (
                            normalized_base1_start == trimmed_base1_start
                            and normalized_ref == trimmed_ref.upper()
                            and normalized_alt == trimmed_alt
                        ):
                            return normalized_read_start, normalized_read_end
                i = j
                continue

            i += 1
        return None

    def locus_read_from_pysam_aligned_segment(
        self,
        pysam_aligned_segment,
        base0_start_inclusive,
        base0_end_exclusive,
        trimmed_base1_start=None,
        trimmed_ref=None,
        trimmed_alt=None,
    ):
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
            logger.warning("Read missing name at position %d", base0_start_inclusive + 1)
            return None

        if pysam_aligned_segment.is_unmapped:
            logger.warning("How did we get unmapped read '%s' in a pileup?", name)
            return None

        mapping_quality = pysam_aligned_segment.mapping_quality

        if mapping_quality is None:
            if self.min_mapping_quality > 0:
                logger.debug("Skipping read '%s' due to missing MAPQ" % name)
                return None
            else:
                mapping_quality = 0
        elif mapping_quality < self.min_mapping_quality:
            logger.debug(
                "Skipping read '%s' due to low MAPQ: %d < %d",
                name,
                mapping_quality,
                self.min_mapping_quality,
            )
            return None

        sequence = pysam_aligned_segment.query_sequence

        if sequence is None:
            logger.warning("Skipping read '%s' due to missing sequence" % name)
            return None

        base_qualities = pysam_aligned_segment.query_qualities

        if base_qualities is None:
            logger.warning("Skipping read '%s' due to missing base qualities" % name)
            return None
        elif len(base_qualities) != len(sequence):
            logger.warning(
                "Skipping read '%s' due to mismatch in length of sequence (%d) and qualities (%d)"
                % (name, len(sequence), len(base_qualities))
            )
            return None

        # By default, AlignedSegment.get_reference_positions only returns base-1 positions
        # from the reference that are within the alignment. If full_length is set,
        # None values will be included for any soft-clipped or unaligned positions
        # within the read. The returned list will thus be of the same
        # length as the read.
        base0_reference_positions = pysam_aligned_segment.get_reference_positions(
            full_length=True
        )

        if len(base0_reference_positions) != len(base_qualities):
            logger.warning(
                "Skipping read '%s' due to mismatch in length of positions (%d) and qualities (%d)"
                % (name, len(base0_reference_positions), len(base_qualities))
            )
            return None

        reference_position_to_read_position = {
            ref_pos: read_pos
            for read_pos, ref_pos in enumerate(base0_reference_positions)
            if ref_pos is not None
        }

        aligned_subsequence_start = pysam_aligned_segment.query_alignment_start
        aligned_subsequence_end = pysam_aligned_segment.query_alignment_end

        reference_interval_size = base0_end_exclusive - base0_start_inclusive
        if reference_interval_size < 0:
            raise ValueError("Unexpected interval start after interval end")

        normalized_indel_interval = None
        if (
            trimmed_base1_start is not None
            and trimmed_ref is not None
            and trimmed_alt is not None
        ):
            normalized_indel_interval = self._left_aligned_indel_interval_for_variant(
                pysam_aligned_segment=pysam_aligned_segment,
                trimmed_base1_start=trimmed_base1_start,
                trimmed_ref=trimmed_ref,
                trimmed_alt=trimmed_alt,
            )

        # TODO:
        #  Consider how to handle variants before splice sites, where
        #  the bases before or after on the genome will not be mapped on the
        #  read
        #
        # we have a dictionary mapping base-1 reference positions to base-0
        # read indices and we need to use that to convert the reference
        # half-open interval into a half-open interval on the read.
        if normalized_indel_interval is not None:
            read_base0_start_inclusive, read_base0_end_exclusive = (
                normalized_indel_interval
            )
        elif reference_interval_size == 0:
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
            read_base0_before_insertion = reference_position_to_read_position.get(
                reference_position_before_insertion
            )
            read_base0_after_insertion = reference_position_to_read_position.get(
                reference_position_after_insertion
            )

            if (
                read_base0_before_insertion is None
                and read_base0_after_insertion is None
            ):
                return None
            elif read_base0_before_insertion is None:
                read_base0_start_inclusive = aligned_subsequence_start
                read_base0_end_exclusive = read_base0_after_insertion
            elif read_base0_after_insertion is None:
                read_base0_start_inclusive = read_base0_before_insertion + 1
                read_base0_end_exclusive = aligned_subsequence_end
            elif read_base0_after_insertion - read_base0_before_insertion == 1:
                read_base0_start_inclusive = read_base0_end_exclusive = (
                    read_base0_before_insertion + 1
                )
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
            if base0_start_inclusive in reference_position_to_read_position:
                read_base0_start_inclusive = reference_position_to_read_position[
                    base0_start_inclusive
                ]
            elif base0_start_inclusive - 1 in reference_position_to_read_position:
                # if first base of reference locus isn't mapped, try getting the base
                # before it and then adding one to its corresponding base index
                read_base0_position_before_locus = reference_position_to_read_position[
                    base0_start_inclusive - 1
                ]
                read_base0_start_inclusive = read_base0_position_before_locus + 1
            else:
                return None

            if base0_end_exclusive in reference_position_to_read_position:
                read_base0_end_exclusive = reference_position_to_read_position[
                    base0_end_exclusive
                ]
            elif (base0_end_exclusive - 1) in reference_position_to_read_position:
                # if exclusive last index of reference interval doesn't have a corresponding
                # base position then try getting the base position of the reference
                # position before it and then adding one
                read_base0_end_inclusive = reference_position_to_read_position[
                    base0_end_exclusive - 1
                ]
                read_base0_end_exclusive = read_base0_end_inclusive + 1
            else:
                return None

        if isinstance(sequence, bytes):
            sequence = sequence.decode("ascii")

        if not self.use_soft_clipped_bases:
            # if we're not allowing soft clipped based then
            # the fraction of the read which is usable may be smaller
            # than the sequence, qualities, and alignment positions
            # we've extracted, so slice through those to get rid of
            # soft-clipped ends of the read
            sequence = sequence[aligned_subsequence_start:aligned_subsequence_end]
            base0_reference_positions = base0_reference_positions[
                aligned_subsequence_start:aligned_subsequence_end
            ]
            base_qualities = base_qualities[
                aligned_subsequence_start:aligned_subsequence_end
            ]
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
            read_base0_end_exclusive=read_base0_end_exclusive,
            source_read_count=1,
        )

    @staticmethod
    def _locus_read_sort_key(locus_read):
        mapped_reference_positions = [
            ref_pos for ref_pos in locus_read.reference_positions if ref_pos is not None
        ]
        if mapped_reference_positions:
            return (
                mapped_reference_positions[0],
                mapped_reference_positions[-1],
                len(locus_read.sequence),
            )
        return (float("inf"), float("inf"), len(locus_read.sequence))

    @staticmethod
    def _base_tokens_from_locus_read(locus_read):
        sequence = locus_read.sequence
        reference_positions = locus_read.reference_positions
        quality_scores = list(locus_read.quality_scores)
        index_to_key = {}
        tokens = []

        i = 0
        n = len(sequence)
        while i < n:
            reference_position = reference_positions[i]
            if reference_position is not None:
                key = (2 * reference_position, 0)
                tokens.append((key, reference_position, sequence[i], quality_scores[i]))
                index_to_key[i] = key
                i += 1
                continue

            j = i
            while j < n and reference_positions[j] is None:
                j += 1

            previous_reference_position = None
            for k in range(i - 1, -1, -1):
                if reference_positions[k] is not None:
                    previous_reference_position = reference_positions[k]
                    break

            next_reference_position = None
            for k in range(j, n):
                if reference_positions[k] is not None:
                    next_reference_position = reference_positions[k]
                    break

            if previous_reference_position is None and next_reference_position is None:
                return None, None

            if previous_reference_position is None:
                anchor = 2 * next_reference_position - 1
            else:
                anchor = 2 * previous_reference_position + 1

            for offset, k in enumerate(range(i, j)):
                key = (anchor, offset)
                tokens.append((key, None, sequence[k], quality_scores[k]))
                index_to_key[k] = key
            i = j

        return tokens, index_to_key

    @staticmethod
    def _translate_read_interval(locus_read, index_to_key, merged_index_by_key):
        start = locus_read.read_base0_start_inclusive
        end = locus_read.read_base0_end_exclusive

        if start is None or end is None:
            return None

        if start < end:
            merged_indices = [
                merged_index_by_key[index_to_key[i]]
                for i in range(start, end)
            ]
            return min(merged_indices), max(merged_indices) + 1

        if start == 0:
            boundary = 0
        elif start == len(locus_read.sequence):
            boundary = len(merged_index_by_key)
        else:
            boundary = merged_index_by_key[index_to_key[start]]
        return boundary, boundary

    @classmethod
    def _merge_locus_read_pair(cls, first, second):
        if (
            first.name != second.name
            or first.reference_base0_start_inclusive != second.reference_base0_start_inclusive
            or first.reference_base0_end_exclusive != second.reference_base0_end_exclusive
        ):
            return None

        first_tokens, first_index_to_key = cls._base_tokens_from_locus_read(first)
        second_tokens, second_index_to_key = cls._base_tokens_from_locus_read(second)
        if first_tokens is None or second_tokens is None:
            return None

        overlapping_keys = {key for key, _, _, _ in first_tokens}.intersection(
            key for key, _, _, _ in second_tokens
        )
        if not overlapping_keys:
            return None

        merged_tokens = {}
        for token in first_tokens + second_tokens:
            key, reference_position, base, quality_score = token
            existing = merged_tokens.get(key)
            if existing is None:
                merged_tokens[key] = token
                continue

            _, existing_reference_position, existing_base, existing_quality_score = existing
            if existing_reference_position != reference_position:
                return None

            if base == existing_base:
                merged_tokens[key] = (
                    key,
                    reference_position,
                    base,
                    max(existing_quality_score, quality_score),
                )
            elif quality_score > existing_quality_score:
                merged_tokens[key] = token
            elif quality_score == existing_quality_score:
                return None

        merged_sequence = []
        merged_reference_positions = []
        merged_quality_scores = []
        merged_index_by_key = {}

        for merged_index, key in enumerate(sorted(merged_tokens)):
            _, reference_position, base, quality_score = merged_tokens[key]
            merged_sequence.append(base)
            merged_reference_positions.append(reference_position)
            merged_quality_scores.append(quality_score)
            merged_index_by_key[key] = merged_index

        translated_intervals = [
            translated_interval
            for translated_interval in [
                cls._translate_read_interval(first, first_index_to_key, merged_index_by_key),
                cls._translate_read_interval(second, second_index_to_key, merged_index_by_key),
            ]
            if translated_interval is not None
        ]
        if not translated_intervals:
            return None

        non_empty_intervals = [
            interval for interval in translated_intervals if interval[0] != interval[1]
        ]
        if non_empty_intervals:
            read_base0_start_inclusive = min(start for start, _ in non_empty_intervals)
            read_base0_end_exclusive = max(end for _, end in non_empty_intervals)
        else:
            read_base0_start_inclusive = read_base0_end_exclusive = translated_intervals[0][0]

        return LocusRead(
            name=first.name,
            sequence="".join(merged_sequence),
            reference_positions=merged_reference_positions,
            quality_scores=merged_quality_scores,
            reference_base0_start_inclusive=first.reference_base0_start_inclusive,
            reference_base0_end_exclusive=first.reference_base0_end_exclusive,
            read_base0_start_inclusive=read_base0_start_inclusive,
            read_base0_end_exclusive=read_base0_end_exclusive,
            source_read_count=first.source_read_count + second.source_read_count,
        )

    @classmethod
    def _merge_overlapping_locus_reads(cls, reads):
        grouped_reads = defaultdict(list)
        for read in reads:
            grouped_reads[read.name].append(read)

        merged_reads = []
        for grouped_read_list in grouped_reads.values():
            pending_reads = sorted(grouped_read_list, key=cls._locus_read_sort_key)
            changed = True
            while changed and len(pending_reads) > 1:
                changed = False
                next_round = []
                used = [False] * len(pending_reads)
                for i, read in enumerate(pending_reads):
                    if used[i]:
                        continue
                    merged_read = read
                    for j in range(i + 1, len(pending_reads)):
                        if used[j]:
                            continue
                        candidate = cls._merge_locus_read_pair(merged_read, pending_reads[j])
                        if candidate is None:
                            continue
                        merged_read = candidate
                        used[j] = True
                        changed = True
                    used[i] = True
                    next_round.append(merged_read)
                pending_reads = sorted(next_round, key=cls._locus_read_sort_key)
            merged_reads.extend(pending_reads)
        return merged_reads

    def get_locus_reads(
        self,
        alignment_file,
        chromosome,
        base0_start_inclusive,
        base0_end_exclusive,
        trimmed_base1_start=None,
        trimmed_ref=None,
        trimmed_alt=None,
    ):
        """
        Create LocusRead objects for reads which overlap the given chromosome,
        start, and end positions. The actual work to figure out if what's between
        those positions matches a variant happens later when LocusRead objects are
        converted to AlleleRead objects.

        If `merge_overlapping_fragments` is enabled then overlapping paired-end
        reads from the same fragment are conservatively merged so downstream
        assembly sees one fragment-spanning sequence instead of double-counting
        the overlap. The raw number of source alignments is retained on each
        merged LocusRead via `source_read_count`.

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
            base0_end_exclusive,
        )
        reads = []
        total_count = 0
        # check overlap against wider overlap to make sure we don't miss
        # any reads
        base0_pos_before_start = max(0, base0_start_inclusive - 1)
        base0_pos_after_end = base0_end_exclusive + 1
        for aligned_segment in alignment_file.fetch(
            chromosome, base0_pos_before_start, base0_pos_after_end
        ):
            total_count += 1
            # we get a significant speed up if we skip reads that have spliced
            # out the entire interval of interest. this is redundant with the
            # attempt to find mapping positions in
            #   self.locus_read_from_pysam_aligned_segment
            # but we do it here to skip the function call overhead for loci
            # where ~1M reads are mapped
            if (
                aligned_segment.get_overlap(base0_pos_before_start, base0_pos_after_end)
                == 0
            ):
                continue
            read = self.locus_read_from_pysam_aligned_segment(
                aligned_segment,
                base0_start_inclusive=base0_start_inclusive,
                base0_end_exclusive=base0_end_exclusive,
                trimmed_base1_start=trimmed_base1_start,
                trimmed_ref=trimmed_ref,
                trimmed_alt=trimmed_alt,
            )
            if read is not None:
                reads.append(read)
        logger.info(
            "Kept %d/%d reads overlapping locus %s:%d-%d",
            len(reads),
            total_count,
            chromosome,
            base0_start_inclusive,
            base0_end_exclusive,
        )
        if not self.merge_overlapping_fragments:
            return reads

        merged_reads = self._merge_overlapping_locus_reads(reads)
        logger.info(
            "Merged overlapping paired reads into %d locus reads",
            len(merged_reads),
        )
        return merged_reads

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

    def locus_reads_overlapping_variant(self, alignment_file, variant, chromosome=None):
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
                valid_chromosome_names=set(alignment_file.references),
            )

        if chromosome is None:
            # failed to infer a chromsome name for this variant which
            # matches names used in SAM/BAM file
            logger.warning(
                "Chromosome '%s' from variant %s not in alignment file %s",
                variant.contig,
                variant,
                alignment_file.filename,
            )
            return []

        logger.info(
            "Gathering reads for variant %s (with gene names %s)",
            variant,
            variant.gene_names,
        )

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
            base0_end_exclusive=base0_end_exclusive,
            trimmed_base1_start=base1_position,
            trimmed_ref=ref,
            trimmed_alt=alt,
        )

    def allele_reads_overlapping_variant(self, variant, alignment_file):
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
                alignment_file=alignment_file, variant=variant
            )
        )

    def read_evidence_for_variant(self, variant, alignment_file):
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
            variant=variant, alignment_file=alignment_file
        )
        return ReadEvidence.from_variant_and_allele_reads(variant, allele_reads)

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
            variant=variant, alignment_file=alignment_file
        )
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
                variant=variant, alignment_file=alignment_file
            )
            yield variant, read_evidence
