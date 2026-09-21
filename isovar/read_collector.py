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
from collections.abc import Mapping
from itertools import groupby

from .default_parameters import (
    USE_SECONDARY_ALIGNMENTS,
    USE_DUPLICATE_READS,
    MIN_READ_MAPPING_QUALITY,
    USE_SOFT_CLIPPED_BASES,
    USE_READS_WITHOUT_BASE_QUALITIES,
    MERGE_OVERLAPPING_FRAGMENTS,
    INFER_READ_ENDS, TRIM_ADAPTERS, TRIM_POLY_A, READ_END_PROFILE,
)
from .locus_read import LocusRead
from .logging import get_logger
from .allele_read import AlleleRead
from .allele_read_helpers import allele_reads_from_locus_reads
from .variant_helpers import require_literal_variant, trim_variant
from .read_evidence import ReadEvidence
from .read_identity import source_alignments_from_pysam, source_read_ids
from .chimeric_alignment import source_alignment_paths_from_pysam
from .read_end_inference import ReadEndProfile, read_sequence_view_from_alignment

logger = get_logger(__name__)


class _CompactLocusRead(object):
    """Internal locus read with aligned blocks instead of per-base positions."""
    __slots__ = [
        "name",
        "sequence",
        "reference_blocks",
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
            reference_blocks,
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
        self.name = name
        self.sequence = sequence
        self.reference_blocks = tuple(reference_blocks)
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

    @classmethod
    def from_locus_read(cls, read):
        return cls(
            name=read.name,
            sequence=read.sequence,
            reference_blocks=AlleleRead._reference_blocks(read.reference_positions),
            quality_scores=read.quality_scores,
            reference_base0_start_inclusive=read.reference_base0_start_inclusive,
            reference_base0_end_exclusive=read.reference_base0_end_exclusive,
            read_base0_start_inclusive=read.read_base0_start_inclusive,
            read_base0_end_exclusive=read.read_base0_end_exclusive,
            source_read_count=read.source_read_count,
            splice_junctions=read.splice_junctions,
            source_alignments=read.source_alignments,
            is_primary=read.is_primary,
            source_alignment_paths=read.source_alignment_paths,
            source_read_views=read.source_read_views,
        )

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
        merge_overlapping_fragments=MERGE_OVERLAPPING_FRAGMENTS,
        use_reads_without_base_qualities=USE_READS_WITHOUT_BASE_QUALITIES,
        infer_read_ends=INFER_READ_ENDS,
        read_end_profile=READ_END_PROFILE,
        trim_adapters=TRIM_ADAPTERS,
        trim_poly_a=TRIM_POLY_A,
        read_filter=None,
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

        use_reads_without_base_qualities : bool
            Retain sequence/alignment evidence when QUAL is absent. Unknown
            qualities remain None, never inferred from MAPQ. Disagreeing mates
            cannot be quality-resolved when either base has unknown quality.

        infer_read_ends : bool
            Attach original-query adapter/poly-A/T candidates without trimming.
            Supplying a profile or enabling trimming also enables annotation.

        read_end_profile : ReadEndProfile or mapping of RG to ReadEndProfile
            Explicit kit configuration, optionally per read group. Missing groups
            remain unknown; no adapter profile is guessed from the instrument.

        trim_adapters, trim_poly_a : bool
            Opt-in removal from terminal soft clips only. Original alignments,
            aligned bases and insertion evidence are never rewritten.

        read_filter : callable or None
            Optional predicate on each eligible original pysam record. Return
            True to keep it. Used by both small-variant and SV reconstruction,
            before merging or segment building. No platform tags are decoded
            unless the predicate requests them. Missing-tag policy belongs to
            the caller; exceptions propagate rather than silently losing reads.
        """
        self.use_secondary_alignments = use_secondary_alignments
        self.use_duplicate_reads = use_duplicate_reads
        self.min_mapping_quality = min_mapping_quality
        self.use_soft_clipped_bases = use_soft_clipped_bases
        self.merge_overlapping_fragments = merge_overlapping_fragments
        self.use_reads_without_base_qualities = use_reads_without_base_qualities
        if read_end_profile is not None:
            profiles = read_end_profile.values() if isinstance(read_end_profile, Mapping) else (read_end_profile,)
            if any(not isinstance(profile, ReadEndProfile) for profile in profiles):
                raise TypeError("read_end_profile must contain ReadEndProfile objects")
        if trim_adapters and read_end_profile is None:
            raise ValueError("Adapter trimming requires an explicit read_end_profile")
        self.infer_read_ends = infer_read_ends
        self.read_end_profile = read_end_profile
        self.trim_adapters = trim_adapters
        self.trim_poly_a = trim_poly_a
        if read_filter is not None and not callable(read_filter):
            raise TypeError("read_filter must be callable or None")
        self.read_filter = read_filter

    def alignment_filter_reason(self, read):
        """Return a rejection reason, or None, for one original SAM record.

        MAPQ retains the historical numeric filter, including STAR's 255 for
        unique mappings. Metadata reports 255 as unavailable, never Q255.
        Unknown QUAL remains optional on every platform.
        """
        if read.is_unmapped:
            return "unmapped"
        if read.is_qcfail:
            return "qc_fail"
        if read.is_duplicate and not self.use_duplicate_reads:
            return "duplicate"
        if read.is_secondary and not self.use_secondary_alignments:
            return "secondary"
        if read.query_name is None:
            return "name_unavailable"
        mapq = read.mapping_quality
        if (mapq is None and self.min_mapping_quality > 0) or (
                mapq is not None and mapq < self.min_mapping_quality):
            return "mapping_quality"
        sequence = read.query_sequence
        if not sequence:
            return "sequence_unavailable"
        qualities = read.query_qualities
        if qualities is None:
            if not self.use_reads_without_base_qualities:
                return "qualities_unavailable"
        elif len(qualities) != len(sequence):
            return "qualities_length_mismatch"
        if self.read_filter is not None and not self.read_filter(read):
            return "read_filter"
        return None

    def read_sequence_view(self, read):
        """Infer original-query end structure independently of a variant locus."""
        profile = self.read_end_profile
        if isinstance(profile, Mapping):
            group = read.get_tag("RG") if read.has_tag("RG") else ""
            profile = profile.get(group)
        return read_sequence_view_from_alignment(
            read, profile, trim_adapters=self.trim_adapters, trim_poly_a=self.trim_poly_a)

    @staticmethod
    def _iter_left_aligned_indel_events(
        preceding_reference,
        base1_start,
        ref,
        alt,
        read_start,
        read_end,
    ):
        """
        Yield an indel and its left shifts within the immediately preceding
        contiguous aligned block. Never jump across another CIGAR event.

        The query interval is shifted in lockstep with the indel representation
        so that downstream AlleleRead objects see the same canonical split for
        every equivalent alignment of the same indel.
        """
        ref = ref.upper()
        alt = alt.upper()

        yield base1_start, ref, alt, read_start, read_end
        for previous_ref_base in reversed(preceding_reference):
            if base1_start <= 1:
                break
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
            yield base1_start, ref, alt, read_start, read_end

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

        read = pysam_aligned_segment
        query_sequence = read.query_sequence
        if query_sequence is None:
            return None
        if isinstance(query_sequence, bytes):
            query_sequence = query_sequence.decode("ascii")
        # MD describes M/=/X and D, but not N. Reconstructing this compact
        # reference avoids expanding potentially megabase introns into pairs
        # and keeps unavailable intronic bases out of deletion sequences.
        try:
            reference_sequence = read.get_reference_sequence().upper()
        except ValueError:
            logger.debug(
                "Skipping indel left-alignment for read '%s' because the alignment "
                "does not expose reference bases (for example, missing MD tag)",
                pysam_aligned_segment.query_name,
            )
            return None

        cigar = read.cigartuples or []
        if len(reference_sequence) != sum(n for op, n in cigar if op in (0, 2, 7, 8)):
            logger.debug("Skipping indel normalization for inconsistent MD/CIGAR on '%s'", read.query_name)
            return None
        query_pos, ref_pos, md_pos = 0, read.reference_start, 0
        preceding_reference = ""
        # Adjacent identical operations are legal and describe one event.
        for operation, group in groupby(cigar, key=lambda pair: pair[0]):
            length = sum(n for _, n in group)
            if operation in (0, 7, 8):
                preceding_reference += reference_sequence[md_pos:md_pos + length]
                query_pos += length
                ref_pos += length
                md_pos += length
                continue
            if (is_insertion and operation == 1) or (is_deletion and operation == 2):
                event_ref = reference_sequence[md_pos:md_pos + length] if operation == 2 else ""
                event_alt = query_sequence[query_pos:query_pos + length] if operation == 1 else ""
                for position, ref, alt, start, end in cls._iter_left_aligned_indel_events(
                        preceding_reference=preceding_reference,
                        base1_start=ref_pos + (operation == 2), ref=event_ref, alt=event_alt,
                        read_start=query_pos, read_end=query_pos + (length if operation == 1 else 0)):
                    if (position == trimmed_base1_start and ref == trimmed_ref.upper()
                            and alt == trimmed_alt.upper()):
                        return start, end
            # N, I, D, S, H and P all break the contiguous shifting context.
            preceding_reference = ""
            if operation in (1, 4):
                query_pos += length
            elif operation == 2:
                ref_pos += length
                md_pos += length
            elif operation == 3:
                ref_pos += length
        return None

    @staticmethod
    def _interval_overlaps_reference_skip(read, start, end):
        """Reject skipped reference bases, without rejecting unrelated splices.

        An insertion has an empty interval and needs the two neighboring
        reference positions; its boundary must not be in or touch an intron.
        Nonempty intervals may end immediately before or start after a skip.
        """
        pos = read.reference_start
        for operation, length in read.cigartuples or []:
            if operation == 3:
                if (start == end and pos <= start <= pos + length) or (
                        start < end and start < pos + length and end > pos):
                    return True
            if operation in (0, 2, 3, 7, 8):
                pos += length
        return False

    @staticmethod
    def _splice_junctions(read):
        """Return CIGAR ``N`` intervals in 0-based half-open coordinates."""
        reference_position = read.reference_start
        junctions = []
        for operation, length in read.cigartuples or ():
            if operation == 3:
                junctions.append((reference_position, reference_position + length))
            if operation in (0, 2, 3, 7, 8):
                reference_position += length
        return tuple(junctions)

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
        name = pysam_aligned_segment.query_name
        reason = self.alignment_filter_reason(pysam_aligned_segment)
        if reason is not None:
            logger.debug("Skipping read %r: %s", name, reason)
            return None

        if self._interval_overlaps_reference_skip(
                pysam_aligned_segment, base0_start_inclusive, base0_end_exclusive):
            return None

        sequence = pysam_aligned_segment.query_sequence
        base_qualities = pysam_aligned_segment.query_qualities
        if base_qualities is None:
            base_qualities = [None] * len(sequence)

        # By default, AlignedSegment.get_reference_positions returns only the
        # 0-based reference positions aligned to read bases. If full_length is
        # set, None values are included for soft-clipped or unaligned read
        # positions, so the returned list has the same length as the read.
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
            # Only CIGAR I supplies inserted bases. Missing anchors can be a
            # read end, clipping or D; none makes mapped flanks an insertion.
            # This also retains terminal/adjacent I operations without MD.
            read_base0_start_inclusive = read_base0_end_exclusive = None
            query_position = 0
            reference_position = pysam_aligned_segment.reference_start
            for operation, length in pysam_aligned_segment.cigartuples or ():
                if operation == 1 and reference_position == base0_start_inclusive:
                    if read_base0_start_inclusive is None:
                        read_base0_start_inclusive = query_position
                    read_base0_end_exclusive = query_position + length
                if operation in (0, 1, 4, 7, 8):
                    query_position += length
                if operation in (0, 2, 3, 7, 8):
                    reference_position += length
            if read_base0_start_inclusive is None:
                before = reference_position_to_read_position.get(base0_start_inclusive - 1)
                after = reference_position_to_read_position.get(base0_start_inclusive)
                if before is not None and after == before + 1:
                    read_base0_start_inclusive = read_base0_end_exclusive = after
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

        query_interval = (None if read_base0_start_inclusive is None else
                          (read_base0_start_inclusive, read_base0_end_exclusive))
        view = None
        retained_start, retained_end = 0, len(sequence)
        if (self.infer_read_ends or self.read_end_profile is not None
                or self.trim_adapters or self.trim_poly_a):
            view = self.read_sequence_view(pysam_aligned_segment)
            retained_start, retained_end = view.start, view.end
        if not self.use_soft_clipped_bases:
            retained_start = max(retained_start, aligned_subsequence_start)
            retained_end = min(retained_end, aligned_subsequence_end)
        # One slice and one coordinate shift for both policies. Allele and SA
        # intervals above came from the unchanged original CIGAR/SEQ.
        sequence = sequence[retained_start:retained_end]
        base0_reference_positions = base0_reference_positions[retained_start:retained_end]
        base_qualities = base_qualities[retained_start:retained_end]
        if read_base0_start_inclusive is not None:
            read_base0_start_inclusive -= retained_start
        if read_base0_end_exclusive is not None:
            read_base0_end_exclusive -= retained_start
        if view is not None:
            view = view._replace(start=retained_start, end=retained_end)
        source_alignments = source_alignments_from_pysam(pysam_aligned_segment, name)
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
            splice_junctions=self._splice_junctions(pysam_aligned_segment),
            source_alignments=source_alignments,
            source_alignment_paths=source_alignment_paths_from_pysam(
                pysam_aligned_segment, source_alignments, query_interval),
            source_read_views=() if view is None else ((source_alignments[0], view),),
            is_primary=not (pysam_aligned_segment.is_secondary
                            or pysam_aligned_segment.is_supplementary),
        )

    @staticmethod
    def _locus_read_sort_key(locus_read):
        if isinstance(locus_read, _CompactLocusRead):
            if locus_read.reference_blocks:
                return (
                    locus_read.reference_blocks[0][2],
                    locus_read.reference_blocks[-1][3] - 1,
                    len(locus_read.sequence),
                )
            return (float("inf"), float("inf"), len(locus_read.sequence))
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
        if isinstance(locus_read, _CompactLocusRead):
            return ReadCollector._base_tokens_from_compact_locus_read(locus_read)

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
    def _base_tokens_from_compact_locus_read(locus_read):
        """Tokenize aligned blocks without reconstructing a coordinate list."""
        sequence = locus_read.sequence
        quality_scores = locus_read.quality_scores
        tokens = []
        index_to_key = {}
        query_position = 0
        previous_reference_position = None

        for query_start, query_end, reference_start, reference_end in (
                locus_read.reference_blocks):
            if query_position < query_start:
                anchor = (
                    2 * previous_reference_position + 1
                    if previous_reference_position is not None
                    else 2 * reference_start - 1)
                for offset, index in enumerate(range(query_position, query_start)):
                    key = (anchor, offset)
                    tokens.append((
                        key, None, sequence[index], quality_scores[index]))
                    index_to_key[index] = key

            for index, reference_position in zip(
                    range(query_start, query_end),
                    range(reference_start, reference_end)):
                key = (2 * reference_position, 0)
                tokens.append((
                    key,
                    reference_position,
                    sequence[index],
                    quality_scores[index],
                ))
                index_to_key[index] = key
            query_position = query_end
            previous_reference_position = reference_end - 1

        if query_position < len(sequence):
            if previous_reference_position is None:
                return None, None
            anchor = 2 * previous_reference_position + 1
            for offset, index in enumerate(range(query_position, len(sequence))):
                key = (anchor, offset)
                tokens.append((key, None, sequence[index], quality_scores[index]))
                index_to_key[index] = key
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
        # QNAME alone also identifies secondary/supplementary placements of
        # ONE segment. Only complementary primary mates can be collapsed.
        first_ids, second_ids = source_read_ids(first), source_read_ids(second)
        if not (first.is_primary and second.is_primary
                and len(first_ids) == len(second_ids) == 1):
            return None
        first_id, second_id = next(iter(first_ids)), next(iter(second_ids))
        if first_id[:2] != second_id[:2] or {first_id[2], second_id[2]} != {64, 128}:
            return None
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

        first_keys = {key for key, _, _, _ in first_tokens}
        second_keys = {key for key, _, _, _ in second_tokens}
        if not first_keys.intersection(second_keys):
            return None

        # Both reads must describe the same alignment path where they overlap.
        # Unioning discordant paths can fill a deletion or invent an insertion.
        overlap_start = max(min(first_keys), min(second_keys))
        overlap_end = min(max(first_keys), max(second_keys))
        if any(overlap_start <= key <= overlap_end for key in first_keys ^ second_keys):
            return None
        # A splice skip and a deletion can omit exactly the same base tokens.
        for start, end in set(first.splice_junctions) ^ set(second.splice_junctions):
            if overlap_start[0] <= 2 * (start - 1) and 2 * end <= overlap_end[0]:
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
                    (quality_score if existing_quality_score is None else
                     existing_quality_score if quality_score is None else
                     max(existing_quality_score, quality_score)),
                )
            elif quality_score is None or existing_quality_score is None:
                return None
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

        read_kwargs = dict(
            name=first.name,
            sequence="".join(merged_sequence),
            quality_scores=merged_quality_scores,
            reference_base0_start_inclusive=first.reference_base0_start_inclusive,
            reference_base0_end_exclusive=first.reference_base0_end_exclusive,
            read_base0_start_inclusive=read_base0_start_inclusive,
            read_base0_end_exclusive=read_base0_end_exclusive,
            source_read_count=first.source_read_count + second.source_read_count,
            splice_junctions=tuple(sorted(set(
                first.splice_junctions + second.splice_junctions))),
            source_alignments=tuple(sorted(first.source_alignments + second.source_alignments)),
            source_alignment_paths=tuple(sorted(first.source_alignment_paths + second.source_alignment_paths)),
            source_read_views=first.source_read_views + second.source_read_views,
            is_primary=True,
        )
        if isinstance(first, _CompactLocusRead) and isinstance(second, _CompactLocusRead):
            return _CompactLocusRead(
                reference_blocks=AlleleRead._reference_blocks(
                    merged_reference_positions),
                **read_kwargs,
            )
        return LocusRead(
            reference_positions=merged_reference_positions,
            **read_kwargs,
        )

    @classmethod
    def _merge_overlapping_locus_reads(cls, reads):
        grouped_reads = defaultdict(list)
        for read in reads:
            grouped_reads[read.name].append(read)

        merged_reads = []
        for grouped_read_list in grouped_reads.values():
            if len(grouped_read_list) == 1:
                merged_reads.extend(grouped_read_list)
                continue
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
        """Return public, list-backed ``LocusRead`` objects for a locus.

        Equal built-in coordinate integers are shared within this call, but
        every read retains its own mutable ``reference_positions`` list. This
        representation and the subclass conversion/merge hooks are preserved
        for existing API consumers.
        """
        return self.collect_locus_reads_with_optional_compaction(
            alignment_file=alignment_file,
            chromosome=chromosome,
            base0_start_inclusive=base0_start_inclusive,
            base0_end_exclusive=base0_end_exclusive,
            trimmed_base1_start=trimmed_base1_start,
            trimmed_ref=trimmed_ref,
            trimmed_alt=trimmed_alt,
            compact=False,
        )

    def collect_locus_reads_with_optional_compaction(
        self,
        alignment_file,
        chromosome,
        base0_start_inclusive,
        base0_end_exclusive,
        trimmed_base1_start=None,
        trimmed_ref=None,
        trimmed_alt=None,
        compact=False,
    ):
        """
        Collect public LocusReads or compact internal equivalents.

        If `merge_overlapping_fragments` is enabled then overlapping paired-end
        reads from the same fragment are conservatively merged so downstream
        assembly sees one fragment-spanning sequence instead of double-counting
        the overlap. The number of sequenced segments is retained on each
        merged LocusRead via `source_read_count`.

        Parameters
        ----------
        alignment_file : pysam.AlignmentFile

        chromosome : str

        base0_start_inclusive : int
            Start of genomic interval, base 0 and inclusive

        base0_end_exclusive : int
            End of genomic interval, base 0 and exclusive

        ``compact=False`` returns the public list-backed LocusRead form;
        ``compact=True`` returns internal block-backed observations. The main
        pipeline selects compact storage only when collection hooks are
        unmodified. Compaction drops each per-base coordinate list after
        converting it to a handful of aligned blocks.
        """
        logger.debug(
            "Gathering reads at locus %s:%d-%d",
            chromosome,
            base0_start_inclusive,
            base0_end_exclusive,
        )
        reads = []
        total_count = 0
        # Deep reads repeat the same coordinates millions of times. Share only
        # immutable integers, never lists or read evidence, through a cache
        # local to this call. Every cached integer is already held by a read,
        # so the cache adds only table slots, and its size is bounded by the
        # reference span of reads at this locus. A fixed-size LRU cache shared
        # nothing once that span exceeded its bound (#234).
        shared_positions = {} if not compact else None
        share_position = shared_positions.setdefault if shared_positions is not None else None

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
                if compact:
                    read = _CompactLocusRead.from_locus_read(read)
                else:
                    # Hooks may return read-like objects without coordinates.
                    positions = getattr(read, "reference_positions", None)
                    if type(positions) is list:
                        # Keep the original list object, and leave custom sequence
                        # types / non-builtin coordinate objects from hooks alone.
                        positions[:] = [share_position(p, p) if type(p) is int else p
                                        for p in positions]
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
        Resolve a contig label, preferring an exact match to a unique alias.

        This does not establish reference-sequence or assembly equivalence.
        In particular, mitochondrial references can differ in both length
        and sequence despite using M/MT/chrM/chrMT labels.
        Parameters
        ----------
        variant_chromosome_name : str

        valid_chromosome_names : set of str

        Returns
        -------
        str or None
            Matching BAM name, or None if no candidate exists.

        Raises
        ------
        ValueError
            More than one alias exists and there is no exact match.
        """
        if variant_chromosome_name in valid_chromosome_names:
            return variant_chromosome_name
        unprefixed = variant_chromosome_name.lower().removeprefix("chr")
        candidate_names = {unprefixed, "chr" + unprefixed}
        if unprefixed in ("m", "mt"):
            candidate_names.update(("m", "mt", "chrm", "chrmt"))
        matches = {name for name in valid_chromosome_names if name.lower() in candidate_names}
        if len(matches) > 1:
            raise ValueError("Ambiguous alignment contig aliases for %r: %s" % (
                variant_chromosome_name, ", ".join(sorted(matches))))
        return next(iter(matches), None)

    def locus_reads_overlapping_variant(self, alignment_file, variant, chromosome=None):
        """Return public ``LocusRead`` objects overlapping a variant."""
        return self._locus_reads_overlapping_variant(
            alignment_file=alignment_file,
            variant=variant,
            chromosome=chromosome,
            compact=False,
        )

    def _locus_reads_overlapping_variant(
            self, alignment_file, variant, chromosome=None, compact=False):
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
        # Validate before consulting annotation or the alignment handle.
        require_literal_variant(variant)
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

        logger.info("Gathering reads for variant %s:%s %s>%s",
                    variant.contig, variant.start, variant.ref, variant.alt)

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

        kwargs = dict(
            alignment_file=alignment_file,
            chromosome=chromosome,
            base0_start_inclusive=base0_start_inclusive,
            base0_end_exclusive=base0_end_exclusive,
            trimmed_base1_start=base1_position,
            trimmed_ref=ref,
            trimmed_alt=alt,
        )
        if compact:
            return self.collect_locus_reads_with_optional_compaction(compact=True, **kwargs)
        return self.get_locus_reads(**kwargs)

    def _can_use_compact_evidence_path(self):
        """Keep every subclass on the historical public hook path."""
        return type(self) is ReadCollector and not any(
            name in self.__dict__
            for name in (
                "locus_reads_overlapping_variant",
                "get_locus_reads",
                "locus_read_from_pysam_aligned_segment",
                "_merge_overlapping_locus_reads",
            )
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
        if self._can_use_compact_evidence_path():
            locus_reads = self._locus_reads_overlapping_variant(
                alignment_file=alignment_file,
                variant=variant,
                compact=True,
            )
        else:
            locus_reads = self.locus_reads_overlapping_variant(
                alignment_file=alignment_file,
                variant=variant,
            )
        return allele_reads_from_locus_reads(locus_reads)

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
