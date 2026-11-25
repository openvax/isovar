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
Helper functions for constructing and filtering VariantSequence objects
from reads overlapping a variant locus.
"""

from bisect import bisect_left
from collections import defaultdict

import numpy as np

from .allele_read_helpers import group_unique_sequences
from .assembly import merge_identical_sequences
from .logging import get_logger
from .variant_sequence import VariantSequence
from .read_identity import observation_groups, partition_read_alignments, read_sort_key

logger = get_logger(__name__)


def _group_unique_sequences_and_transcript_paths(
        allele_reads,
        max_prefix_size=None,
        max_suffix_size=None):
    """Group identical cDNA observations without merging transcript paths."""
    sequence_groups = group_unique_sequences(
        allele_reads,
        max_prefix_size=max_prefix_size,
        max_suffix_size=max_suffix_size)
    groups = defaultdict(set)
    for (prefix, alt, suffix), reads in sequence_groups.items():
        for read in reads:
            groups[
                prefix,
                alt,
                suffix,
                read.compatible_transcript_ids,
            ].add(read)
    return groups


def initial_variant_sequences_from_reads(
        variant_reads,
        max_nucleotides_before_variant=None,
        max_nucleotides_after_variant=None):
    """
    Get all unique sequences from reads spanning a variant locus. This will
    include partial sequences due to reads starting in the middle of the
    sequence around a variant.

    Parameters
    ----------
    variant_reads : list of AlleleRead objects

    max_nucleotides_before_variant : int or None

    max_nucleotides_after_variant : int or None

    Returns
    -------
    list of VariantSequence
    """
    unique_sequence_groups = _group_unique_sequences_and_transcript_paths(
        variant_reads,
        max_prefix_size=max_nucleotides_before_variant,
        max_suffix_size=max_nucleotides_after_variant)

    return [
        VariantSequence(
            prefix=prefix,
            alt=alt,
            suffix=suffix,
            reads=compatible_reads)
        for ((prefix, alt, suffix, _compatible_transcript_ids), reads)
        in unique_sequence_groups.items()
        for compatible_reads in partition_read_alignments(reads)
    ]


class _FlankIndex:
    """Linear-space index of strings compatible with an anchored flank.

    Longer matches occupy one contiguous lexicographic range. Shorter matches
    are exact ancestors of the query. Reverse upstream flanks so that both
    sides use prefix matching. Store ranges, not a quadratic table of matches.
    """

    def __init__(self, groups, upstream):
        self.entries = sorted(
            ((prefix[::-1] if upstream else suffix, group_id)
             for group_id, (prefix, suffix, reads, transcript_ids) in enumerate(groups)))
        self.words = [word for word, _ in self.entries]
        self.exact = {}
        for i, word in enumerate(self.words):
            start = self.exact[word][0] if word in self.exact else i
            self.exact[word] = (start, i + 1)
        self.lengths = sorted({len(word) for word in self.exact})

    def compatible_ranges(self, query):
        start = bisect_left(self.words, query)
        low, high = start, len(self.words)
        # No alphabet/sentinel assumptions, including the empty query.
        while low < high:
            middle = (low + high) // 2
            if self.words[middle].startswith(query):
                low = middle + 1
            else:
                high = middle
        ranges = [(start, low)]
        for length in self.lengths:
            if length >= len(query):
                break
            interval = self.exact.get(query[:length])
            if interval is not None:
                ranges.append(interval)
        return sum(end - start for start, end in ranges), ranges


_CONFLICTING_PLACEMENT = object()


class SharedReadSupportIndex:
    """Reusable flank-and-alignment support index over a pool of AlleleReads.

    Precomputes linear-space flank indexes, per-flank bitmasks/read-sets,
    conflicting-segment sets, and multi-occurrence segment sets once so
    multiple candidates, context lengths, trimming passes, and transcript read
    groups can reconcile compatible read support without repeating read
    grouping or binary-search indexing.
    """

    def __init__(self, read_groups, max_prefix, max_suffix, *, extra_reads=()):
        self.max_prefix = max_prefix
        self.max_suffix = max_suffix
        normalized_read_groups = defaultdict(set)
        for key, reads in read_groups.items():
            if len(key) == 4:
                normalized_read_groups[key].update(reads)
                continue
            prefix, alt, suffix = key
            for read in reads:
                normalized_read_groups[
                    prefix,
                    alt,
                    suffix,
                    getattr(read, "compatible_transcript_ids", None),
                ].add(read)

        bounded_groups = defaultdict(set)
        for (prefix, alt, suffix, transcript_ids), reads in normalized_read_groups.items():
            # Only index the observable context; retain the ORIGINAL read objects
            # and their real bounds for coverage. Distinct overhangs outside every
            # candidate can safely share an index entry, never a read identity.
            key = (
                prefix[-max_prefix:] if max_prefix else "",
                alt,
                suffix[:max_suffix],
                transcript_ids,
            )
            bounded_groups[key].update(reads)

        # Scan all pool reads (and any extra candidate reads) once to identify
        # the small subset of segments with multiple placements (conflicts) or
        # multiple read entries (needing segment-aware observation grouping).
        segment_placement = {}
        conflicting_segments = set()
        multi_occurrence_segments = set()
        seen_reads = set()

        def _scan_read(read):
            if read in seen_reads:
                return
            seen_reads.add(read)
            for seg_key, placement in getattr(read, "source_alignments", ()):
                existing = segment_placement.get(seg_key)
                if existing is None:
                    segment_placement[seg_key] = placement
                else:
                    multi_occurrence_segments.add(seg_key)
                    if existing is not _CONFLICTING_PLACEMENT and existing != placement:
                        segment_placement[seg_key] = _CONFLICTING_PLACEMENT
                        conflicting_segments.add(seg_key)

        for reads in bounded_groups.values():
            for read in reads:
                _scan_read(read)
        for read in extra_reads:
            _scan_read(read)

        self.conflicting_segments = conflicting_segments
        self.multi_occurrence_segments = multi_occurrence_segments

        by_allele = defaultdict(list)

        def group_key(item):
            prefix, alt, suffix, ids = item[0]
            return prefix, alt, suffix, (0, ()) if ids is None else (1, tuple(sorted(ids)))

        for (prefix, alt, suffix, transcript_ids), reads in sorted(
                bounded_groups.items(), key=group_key):
            unique_reads = []
            multi_occ_reads = []
            has_conflict = False
            if conflicting_segments or multi_occurrence_segments:
                for read in reads:
                    placements = getattr(read, "source_alignments", ())
                    if placements:
                        read_conflict = False
                        read_multi = False
                        for seg_key, _ in placements:
                            if seg_key in conflicting_segments:
                                read_conflict = True
                                read_multi = True
                            elif seg_key in multi_occurrence_segments:
                                read_multi = True
                        if read_conflict:
                            has_conflict = True
                            multi_occ_reads.append(read)
                        elif read_multi:
                            multi_occ_reads.append(read)
                        else:
                            unique_reads.append(read)
                    else:
                        unique_reads.append(read)
            else:
                unique_reads = list(reads)

            if has_conflict:
                ordered_group = tuple(sorted(reads, key=read_sort_key))
                unique_tuple = ()
                multi_tuple = ordered_group
            else:
                ordered_group = tuple(unique_reads) + tuple(multi_occ_reads)
                unique_tuple = tuple(unique_reads)
                multi_tuple = tuple(multi_occ_reads)

            by_allele[alt].append((
                prefix,
                suffix,
                ordered_group,
                transcript_ids,
                len(prefix),
                len(suffix),
                has_conflict,
                unique_tuple,
                multi_tuple,
            ))

        self.by_allele = by_allele
        self.indexes = {
            alt: (
                _FlankIndex([(g[0], g[1], g[2], g[3]) for g in groups], True),
                _FlankIndex([(g[0], g[1], g[2], g[3]) for g in groups], False),
            )
            for alt, groups in by_allele.items()
        }
        self._sequence_matches_cache = {}
        # Fast-path bitmask / C-frozenset structures per alt when no placement
        # conflicts exist. Each non-conflicting single-occurrence read gets a
        # bit so coverage counts reduce to O(distinct_flank_lengths) bitwise
        # AND + bit_count() operations, and frozensets are built in pure C via
        # prefix_set & suffix_set & tx_set.
        self._bitmask_state = {}
        if not conflicting_segments:
            for alt, groups in by_allele.items():
                self._bitmask_state[alt] = self._build_allele_bitmask_state(groups)

    @staticmethod
    def _build_allele_bitmask_state(groups):
        read_bit = {}
        next_bit = 0
        group_read_masks = []
        group_multi = []
        prefix_len_masks = defaultdict(int)
        suffix_len_masks = defaultdict(int)
        nonempty_prefix_mask = 0
        nonempty_suffix_mask = 0
        nonempty_prefix_multi_mask = 0
        nonempty_suffix_multi_mask = 0
        prefix_to_reads = defaultdict(set)
        suffix_to_reads = defaultdict(set)
        prefix_to_mask = defaultdict(int)
        suffix_to_mask = defaultdict(int)
        prefix_to_multi_mask = defaultdict(int)
        suffix_to_multi_mask = defaultdict(int)
        tx_to_mask = defaultdict(int)
        tx_to_reads = defaultdict(set)
        tx_to_multi_mask = defaultdict(int)
        all_unique_reads = set()

        for g_idx, (
            prefix,
            suffix,
            _ordered_group,
            transcript_ids,
            prefix_len,
            suffix_len,
            _has_conflict,
            unique_tuple,
            multi_tuple,
        ) in enumerate(groups):
            g_mask = 0
            if unique_tuple:
                for read in unique_tuple:
                    bit = read_bit.get(read)
                    if bit is None:
                        bit = 1 << next_bit
                        read_bit[read] = bit
                        next_bit += 1
                    g_mask |= bit
                prefix_len_masks[prefix_len] |= g_mask
                suffix_len_masks[suffix_len] |= g_mask
                if prefix_len:
                    nonempty_prefix_mask |= g_mask
                if suffix_len:
                    nonempty_suffix_mask |= g_mask
                prefix_to_mask[prefix] |= g_mask
                suffix_to_mask[suffix] |= g_mask
                prefix_to_reads[prefix].update(unique_tuple)
                suffix_to_reads[suffix].update(unique_tuple)
                tx_to_mask[transcript_ids] |= g_mask
                tx_to_reads[transcript_ids].update(unique_tuple)
                all_unique_reads.update(unique_tuple)
            group_read_masks.append(g_mask)
            if multi_tuple:
                m_bit = 1 << len(group_multi)
                group_multi.append((multi_tuple, transcript_ids))
                prefix_to_multi_mask[prefix] |= m_bit
                suffix_to_multi_mask[suffix] |= m_bit
                tx_to_multi_mask[transcript_ids] |= m_bit
                if prefix_len:
                    nonempty_prefix_multi_mask |= m_bit
                if suffix_len:
                    nonempty_suffix_multi_mask |= m_bit

        # Build lightweight flank indexes over distinct prefixes/suffixes only.
        distinct_prefixes = sorted(prefix_to_mask.keys() | prefix_to_multi_mask.keys())
        distinct_suffixes = sorted(suffix_to_mask.keys() | suffix_to_multi_mask.keys())
        prefix_index = _FlankIndex(
            [(p, "", (), None) for p in distinct_prefixes], True
        )
        suffix_index = _FlankIndex(
            [("", s, (), None) for s in distinct_suffixes], False
        )
        return {
            "distinct_prefixes": distinct_prefixes,
            "distinct_suffixes": distinct_suffixes,
            "prefix_index": prefix_index,
            "suffix_index": suffix_index,
            "prefix_to_mask": prefix_to_mask,
            "suffix_to_mask": suffix_to_mask,
            "prefix_to_reads": prefix_to_reads,
            "suffix_to_reads": suffix_to_reads,
            "prefix_to_multi_mask": prefix_to_multi_mask,
            "suffix_to_multi_mask": suffix_to_multi_mask,
            "tx_to_mask": tx_to_mask,
            "tx_to_reads": tx_to_reads,
            "tx_to_multi_mask": tx_to_multi_mask,
            "tx_keys": tuple(tx_to_mask.keys() | tx_to_multi_mask.keys()),
            "group_multi": group_multi,
            "prefix_len_masks": tuple(sorted(prefix_len_masks.items())),
            "suffix_len_masks": tuple(sorted(suffix_len_masks.items())),
            "nonempty_prefix_mask": nonempty_prefix_mask,
            "nonempty_suffix_mask": nonempty_suffix_mask,
            "nonempty_prefix_multi_mask": nonempty_prefix_multi_mask,
            "nonempty_suffix_multi_mask": nonempty_suffix_multi_mask,
            "nonempty_prefix_reads": frozenset(
                r for p, rs in prefix_to_reads.items() if p for r in rs
            ),
            "nonempty_suffix_reads": frozenset(
                r for s, rs in suffix_to_reads.items() if s for r in rs
            ),
            "all_unique_reads": frozenset(all_unique_reads),
            # Per-query caches (amortized across all candidates & tx groups)
            "prefix_mask_cache": {},
            "suffix_mask_cache": {},
            "prefix_reads_cache": {},
            "suffix_reads_cache": {},
            "tx_compat_cache": {},
        }

    def _prefix_compat_masks(self, state, seq_prefix):
        cache = state["prefix_mask_cache"]
        cached = cache.get(seq_prefix)
        if cached is not None:
            return cached
        distinct_prefixes = state["distinct_prefixes"]
        p_to_mask = state["prefix_to_mask"]
        p_to_multi = state["prefix_to_multi_mask"]
        if len(distinct_prefixes) == 1:
            p = distinct_prefixes[0]
            if p.endswith(seq_prefix) or seq_prefix.endswith(p):
                res = (p_to_mask.get(p, 0), p_to_multi.get(p, 0), (p,))
            else:
                res = (0, 0, ())
            cache[seq_prefix] = res
            return res
        p_index = state["prefix_index"]
        _, ranges = p_index.compatible_ranges(seq_prefix[::-1])
        entries = p_index.entries
        mask = 0
        multi_mask = 0
        matched_prefixes = []
        for start, end in ranges:
            for i in range(start, end):
                p = distinct_prefixes[entries[i][1]]
                matched_prefixes.append(p)
                mask |= p_to_mask.get(p, 0)
                multi_mask |= p_to_multi.get(p, 0)
        res = (mask, multi_mask, tuple(matched_prefixes))
        cache[seq_prefix] = res
        return res

    def _suffix_compat_masks(self, state, seq_suffix):
        cache = state["suffix_mask_cache"]
        cached = cache.get(seq_suffix)
        if cached is not None:
            return cached
        distinct_suffixes = state["distinct_suffixes"]
        s_to_mask = state["suffix_to_mask"]
        s_to_multi = state["suffix_to_multi_mask"]
        if len(distinct_suffixes) == 1:
            s = distinct_suffixes[0]
            if s.startswith(seq_suffix) or seq_suffix.startswith(s):
                res = (s_to_mask.get(s, 0), s_to_multi.get(s, 0), (s,))
            else:
                res = (0, 0, ())
            cache[seq_suffix] = res
            return res
        s_index = state["suffix_index"]
        _, ranges = s_index.compatible_ranges(seq_suffix)
        entries = s_index.entries
        mask = 0
        multi_mask = 0
        matched_suffixes = []
        for start, end in ranges:
            for i in range(start, end):
                s = distinct_suffixes[entries[i][1]]
                matched_suffixes.append(s)
                mask |= s_to_mask.get(s, 0)
                multi_mask |= s_to_multi.get(s, 0)
        res = (mask, multi_mask, tuple(matched_suffixes))
        cache[seq_suffix] = res
        return res

    def _prefix_compat_reads(self, state, seq_prefix):
        cache = state["prefix_reads_cache"]
        cached = cache.get(seq_prefix)
        if cached is not None:
            return cached
        _, _, matched_prefixes = self._prefix_compat_masks(state, seq_prefix)
        p_to_reads = state["prefix_to_reads"]
        reads = set()
        for p in matched_prefixes:
            rs = p_to_reads.get(p)
            if rs:
                reads.update(rs)
        result = frozenset(reads)
        cache[seq_prefix] = result
        return result

    def _suffix_compat_reads(self, state, seq_suffix):
        cache = state["suffix_reads_cache"]
        cached = cache.get(seq_suffix)
        if cached is not None:
            return cached
        _, _, matched_suffixes = self._suffix_compat_masks(state, seq_suffix)
        s_to_reads = state["suffix_to_reads"]
        reads = set()
        for s in matched_suffixes:
            rs = s_to_reads.get(s)
            if rs:
                reads.update(rs)
        result = frozenset(reads)
        cache[seq_suffix] = result
        return result

    def _tx_compat(self, state, seq_tx_ids):
        cache = state["tx_compat_cache"]
        cached = cache.get(seq_tx_ids)
        if cached is not None:
            return cached
        tx_to_mask = state["tx_to_mask"]
        tx_to_reads = state["tx_to_reads"]
        tx_to_multi = state["tx_to_multi_mask"]
        mask = 0
        multi_mask = 0
        reads = set()
        matched_txs = []
        for tx_ids in state["tx_keys"]:
            if seq_tx_ids is None:
                ok = tx_ids is None
            else:
                ok = tx_ids is None or seq_tx_ids.issubset(tx_ids)
            if ok:
                mask |= tx_to_mask.get(tx_ids, 0)
                multi_mask |= tx_to_multi.get(tx_ids, 0)
                rs = tx_to_reads.get(tx_ids)
                if rs:
                    reads.update(rs)
                if tx_ids is not None:
                    matched_txs.append(tx_ids)
        result = (mask, multi_mask, frozenset(reads), tuple(matched_txs))
        cache[seq_tx_ids] = result
        return result

    def _fast_coverage_for_key(self, seq_prefix, seq_alt, seq_suffix, seq_tx_ids):
        """Compute (unique_mask, multi_mask, coverage_array) without building read sets."""
        seq_len = len(seq_prefix) + len(seq_alt) + len(seq_suffix)
        if seq_alt not in self._bitmask_state:
            return 0, 0, np.zeros(seq_len, dtype="int32")
        state = self._bitmask_state[seq_alt]
        p_mask, p_multi, _ = self._prefix_compat_masks(state, seq_prefix)
        s_mask, s_multi, _ = self._suffix_compat_masks(state, seq_suffix)
        t_mask, t_multi, _, _ = self._tx_compat(state, seq_tx_ids)
        mask = p_mask & s_mask & t_mask
        multi_mask = p_multi & s_multi & t_multi
        if not seq_alt:
            if not seq_prefix and not seq_suffix:
                return 0, 0, np.zeros(0, dtype="int32")
            if not seq_prefix:
                mask &= state["nonempty_suffix_mask"]
                multi_mask &= state["nonempty_suffix_multi_mask"]
            elif not seq_suffix:
                mask &= state["nonempty_prefix_mask"]
                multi_mask &= state["nonempty_prefix_multi_mask"]
            else:
                mask &= state["nonempty_prefix_mask"] | state["nonempty_suffix_mask"]
                multi_mask &= (
                    state["nonempty_prefix_multi_mask"]
                    | state["nonempty_suffix_multi_mask"]
                )

        variant_start = len(seq_prefix)
        variant_end = variant_start + len(seq_alt)
        boundaries = [0] * (seq_len + 1)
        if mask:
            for p_len, pm in state["prefix_len_masks"]:
                count = (mask & pm).bit_count()
                if count:
                    b_start = variant_start - p_len
                    boundaries[b_start if b_start > 0 else 0] += count
            for s_len, sm in state["suffix_len_masks"]:
                count = (mask & sm).bit_count()
                if count:
                    b_end = variant_end + s_len
                    boundaries[b_end if b_end < seq_len else seq_len] -= count
        if multi_mask:
            group_multi = state["group_multi"]
            multi_occ_reads = []
            m = multi_mask
            while m:
                lsb = m & -m
                multi_occ_reads.extend(group_multi[lsb.bit_length() - 1][0])
                m ^= lsb
            for obs_group in observation_groups(multi_occ_reads):
                b_start = max(
                    0, variant_start - max(len(read.prefix) for read in obs_group)
                )
                b_end = min(
                    seq_len, variant_end + max(len(read.suffix) for read in obs_group)
                )
                boundaries[b_start] += 1
                boundaries[b_end] -= 1
        return mask, multi_mask, np.cumsum(boundaries, dtype="int32")[:-1]

    def _fast_reads_for_key(
            self,
            seq_prefix,
            seq_alt,
            seq_suffix,
            seq_tx_ids,
            unique_mask,
            multi_mask,
            *,
            preserve_candidate_transcript_ids=False,
            frozenset_cache=None):
        """Materialize (reads_frozenset, final_tx_ids) via C-level frozenset intersection."""
        if seq_alt not in self._bitmask_state:
            return frozenset(), (seq_tx_ids if preserve_candidate_transcript_ids else None)
        state = self._bitmask_state[seq_alt]
        cache_key = (
            seq_alt,
            unique_mask,
            multi_mask,
            seq_tx_ids if preserve_candidate_transcript_ids else None,
        )
        if frozenset_cache is not None and cache_key in frozenset_cache:
            return frozenset_cache[cache_key]

        if unique_mask:
            p_reads = self._prefix_compat_reads(state, seq_prefix)
            s_reads = self._suffix_compat_reads(state, seq_suffix)
            _, _, t_reads, _ = self._tx_compat(state, seq_tx_ids)
            # Intersect smallest-first in pure C (same object pointers -> 0 __eq__ calls).
            sets = sorted((p_reads, s_reads, t_reads), key=len)
            unique_fs = sets[0] & sets[1]
            if len(unique_fs):
                unique_fs = unique_fs & sets[2]
            if not seq_alt and unique_fs:
                if not seq_prefix:
                    unique_fs = unique_fs & state["nonempty_suffix_reads"]
                elif not seq_suffix:
                    unique_fs = unique_fs & state["nonempty_prefix_reads"]
                else:
                    unique_fs = unique_fs & (
                        state["nonempty_prefix_reads"] | state["nonempty_suffix_reads"]
                    )
        else:
            unique_fs = frozenset()

        if multi_mask:
            group_multi = state["group_multi"]
            extra = []
            m = multi_mask
            while m:
                lsb = m & -m
                extra.extend(group_multi[lsb.bit_length() - 1][0])
                m ^= lsb
            reads_fs = unique_fs | frozenset(extra) if unique_fs else frozenset(extra)
        else:
            reads_fs = unique_fs

        if preserve_candidate_transcript_ids:
            final_tx_ids = seq_tx_ids
        else:
            final_tx_ids = None
            for read in reads_fs:
                ids = getattr(read, "compatible_transcript_ids", None)
                if ids is not None:
                    if final_tx_ids is None:
                        final_tx_ids = ids
                    elif ids is not final_tx_ids:
                        final_tx_ids = final_tx_ids.intersection(ids)

        res = (reads_fs, final_tx_ids)
        if frozenset_cache is not None:
            frozenset_cache[cache_key] = res
        return res

    def trim_and_support_no_conflicts(
            self,
            variant_sequences,
            min_variant_sequence_coverage,
            *,
            preserve_candidate_transcript_ids=False):
        """Fast conflict-free trim+support: trims via bitmask coverage before building sets."""
        # Deduplicate input candidates by canonical (prefix, alt, suffix, tx_ids).
        unique_keys = {}
        for seq in variant_sequences:
            if len(seq) == 0:
                continue
            tx_ids = seq.compatible_transcript_ids
            tx_key = (0, ()) if tx_ids is None else (1, tuple(sorted(tx_ids)))
            key = (seq.prefix, seq.alt, seq.suffix, tx_key)
            if key not in unique_keys:
                unique_keys[key] = (seq.prefix, seq.alt, seq.suffix, tx_ids)

        # Pre-trim coverage check via bitmasks only — zero frozensets allocated
        # for discarded or shortened candidates.
        surviving = {}
        for key in sorted(unique_keys):
            seq_prefix, seq_alt, seq_suffix, tx_ids = unique_keys[key]
            u_mask, m_mask, cov = self._fast_coverage_for_key(
                seq_prefix, seq_alt, seq_suffix, tx_ids
            )
            mask_ok = cov >= min_variant_sequence_coverage
            if not np.any(mask_ok):
                continue
            indices = np.flatnonzero(mask_ok)
            first_idx = int(indices[0])
            last_idx = int(indices[-1])
            v_start = len(seq_prefix)
            v_end = v_start + len(seq_alt)
            if first_idx > v_start or last_idx + 1 < v_end:
                continue
            new_prefix = seq_prefix[first_idx:]
            new_suffix = seq_suffix[: last_idx - v_end + 1]
            if not (new_prefix or seq_alt or new_suffix):
                continue
            trimmed_key = (new_prefix, seq_alt, new_suffix, key[3])
            was_untrimmed = (
                first_idx == 0 and last_idx + 1 == len(seq_prefix) + len(seq_alt) + len(seq_suffix)
            )
            if trimmed_key not in surviving:
                surviving[trimmed_key] = (
                    new_prefix,
                    seq_alt,
                    new_suffix,
                    tx_ids,
                    (u_mask, m_mask, cov) if was_untrimmed else None,
                )
            elif was_untrimmed and surviving[trimmed_key][4] is None:
                surviving[trimmed_key] = (
                    new_prefix,
                    seq_alt,
                    new_suffix,
                    tx_ids,
                    (u_mask, m_mask, cov),
                )

        # Materialize surviving candidates once in canonical _sequence_key order.
        frozenset_cache = {}
        result = []
        for trimmed_key in sorted(surviving):
            new_prefix, seq_alt, new_suffix, tx_ids, cached_cov = surviving[trimmed_key]
            if cached_cov is not None:
                u_mask, m_mask, cov = cached_cov
            else:
                u_mask, m_mask, cov = self._fast_coverage_for_key(
                    new_prefix, seq_alt, new_suffix, tx_ids
                )
            reads_fs, final_tx_ids = self._fast_reads_for_key(
                new_prefix,
                seq_alt,
                new_suffix,
                tx_ids,
                u_mask,
                m_mask,
                preserve_candidate_transcript_ids=preserve_candidate_transcript_ids,
                frozenset_cache=frozenset_cache,
            )
            result.append(
                VariantSequence._from_validated(
                    prefix=new_prefix,
                    alt=seq_alt,
                    suffix=new_suffix,
                    reads=reads_fs,
                    compatible_transcript_ids=final_tx_ids,
                    source_alignments=None,
                    coverage_cache=cov,
                )
            )
        return result

    def _compatible_group_indices(self, seq_prefix, seq_alt, seq_suffix):
        cache_key = (seq_prefix, seq_alt, seq_suffix)
        cached = self._sequence_matches_cache.get(cache_key)
        if cached is not None:
            return cached
        if seq_alt not in self.indexes:
            self._sequence_matches_cache[cache_key] = ()
            return ()
        upstream, downstream = self.indexes[seq_alt]
        n_left, left_ranges = upstream.compatible_ranges(seq_prefix[::-1])
        n_right, right_ranges = downstream.compatible_ranges(seq_suffix)
        use_left = n_left <= n_right
        index, ranges = (
            (upstream, left_ranges) if use_left else (downstream, right_ranges)
        )
        groups = self.by_allele[seq_alt]
        entries = index.entries
        matched = []
        if use_left:
            for start, end in ranges:
                for i in range(start, end):
                    group_idx = entries[i][1]
                    prefix, suffix = groups[group_idx][0], groups[group_idx][1]
                    if (
                        (seq_alt or (prefix and seq_prefix) or (suffix and seq_suffix))
                        and (suffix.startswith(seq_suffix) or seq_suffix.startswith(suffix))
                    ):
                        matched.append(group_idx)
        else:
            for start, end in ranges:
                for i in range(start, end):
                    group_idx = entries[i][1]
                    prefix, suffix = groups[group_idx][0], groups[group_idx][1]
                    if (
                        (seq_alt or (prefix and seq_prefix) or (suffix and seq_suffix))
                        and (prefix.endswith(seq_prefix) or seq_prefix.endswith(prefix))
                    ):
                        matched.append(group_idx)
        result = tuple(matched)
        self._sequence_matches_cache[cache_key] = result
        return result

    def support_sequences(self, variant_sequences, *, preserve_candidate_transcript_ids=False):
        """Reconcile compatible original reads and exact coverage for candidates."""
        variant_sequences = list(variant_sequences)
        if not variant_sequences:
            return []
        if not self.conflicting_segments:
            frozenset_cache = {}
            result = []
            for sequence in variant_sequences:
                seq_prefix = sequence.prefix
                seq_alt = sequence.alt
                seq_suffix = sequence.suffix
                seq_tx_ids = sequence.compatible_transcript_ids
                u_mask, m_mask, coverage_array = self._fast_coverage_for_key(
                    seq_prefix, seq_alt, seq_suffix, seq_tx_ids
                )
                reads_fs, final_tx_ids = self._fast_reads_for_key(
                    seq_prefix,
                    seq_alt,
                    seq_suffix,
                    seq_tx_ids,
                    u_mask,
                    m_mask,
                    preserve_candidate_transcript_ids=preserve_candidate_transcript_ids,
                    frozenset_cache=frozenset_cache,
                )
                if (
                    len(reads_fs) == len(sequence.reads)
                    and final_tx_ids == seq_tx_ids
                    and reads_fs == sequence.reads
                ):
                    supported = sequence
                    if supported._coverage_cache is None:
                        supported._coverage_cache = coverage_array
                else:
                    supported = VariantSequence._from_validated(
                        prefix=seq_prefix,
                        alt=seq_alt,
                        suffix=seq_suffix,
                        reads=reads_fs,
                        compatible_transcript_ids=final_tx_ids,
                        source_alignments=None,
                        coverage_cache=coverage_array,
                    )
                result.append(supported)
            return result

        conflicting_segments = self.conflicting_segments
        result = []
        for sequence in variant_sequences:
            seq_prefix = sequence.prefix
            seq_alt = sequence.alt
            seq_suffix = sequence.suffix
            seq_len = len(sequence)
            variant_start = len(seq_prefix)
            variant_end = variant_start + len(seq_alt)
            seq_tx_ids = sequence.compatible_transcript_ids

            supporting_reads = set()
            multi_occ_reads = []
            boundaries = [0] * (seq_len + 1)
            conflict_constraints = None
            result_tx_ids = None

            matched_indices = self._compatible_group_indices(seq_prefix, seq_alt, seq_suffix)
            if matched_indices:
                groups = self.by_allele[seq_alt]
                for group_idx in matched_indices:
                    (
                        _prefix,
                        _suffix,
                        ordered_group,
                        transcript_ids,
                        prefix_len,
                        suffix_len,
                        has_conflict,
                        unique_tuple,
                        multi_tuple,
                    ) = groups[group_idx]
                    if seq_tx_ids is None:
                        if transcript_ids is not None:
                            continue
                    elif transcript_ids is not None and not seq_tx_ids.issubset(transcript_ids):
                        continue

                    if not has_conflict:
                        if unique_tuple:
                            supporting_reads.update(unique_tuple)
                            n_unique = len(unique_tuple)
                            b_start = variant_start - prefix_len
                            if b_start < 0:
                                b_start = 0
                            b_end = variant_end + suffix_len
                            if b_end > seq_len:
                                b_end = seq_len
                            boundaries[b_start] += n_unique
                            boundaries[b_end] -= n_unique
                        if multi_tuple:
                            supporting_reads.update(multi_tuple)
                            multi_occ_reads.extend(multi_tuple)
                        if not preserve_candidate_transcript_ids and transcript_ids is not None:
                            if result_tx_ids is None:
                                result_tx_ids = transcript_ids
                            elif transcript_ids is not result_tx_ids:
                                result_tx_ids = result_tx_ids.intersection(transcript_ids)
                    else:
                        if conflict_constraints is None:
                            conflict_constraints = {
                                k: v
                                for k, v in sequence._source_alignments.items()
                                if k in conflicting_segments
                            }
                        added_any = False
                        for read in ordered_group:
                            placements = getattr(read, "source_alignments", ())
                            if all(
                                key not in conflict_constraints
                                or conflict_constraints[key] == value
                                for key, value in placements
                            ):
                                for key, value in placements:
                                    if key in conflicting_segments:
                                        conflict_constraints[key] = value
                                supporting_reads.add(read)
                                multi_occ_reads.append(read)
                                added_any = True
                        if (
                            added_any
                            and not preserve_candidate_transcript_ids
                            and transcript_ids is not None
                        ):
                            if result_tx_ids is None:
                                result_tx_ids = transcript_ids
                            elif transcript_ids is not result_tx_ids:
                                result_tx_ids = result_tx_ids.intersection(transcript_ids)

            if multi_occ_reads:
                for obs_group in observation_groups(multi_occ_reads):
                    b_start = max(
                        0, variant_start - max(len(read.prefix) for read in obs_group)
                    )
                    b_end = min(
                        seq_len, variant_end + max(len(read.suffix) for read in obs_group)
                    )
                    boundaries[b_start] += 1
                    boundaries[b_end] -= 1

            coverage_array = np.cumsum(boundaries, dtype="int32")[:-1]
            final_tx_ids = seq_tx_ids if preserve_candidate_transcript_ids else result_tx_ids

            if (
                len(supporting_reads) == len(sequence.reads)
                and final_tx_ids == seq_tx_ids
                and supporting_reads == sequence.reads
            ):
                supported = sequence
                if supported._coverage_cache is None:
                    supported._coverage_cache = coverage_array
            else:
                supported = VariantSequence._from_validated(
                    prefix=seq_prefix,
                    alt=seq_alt,
                    suffix=seq_suffix,
                    reads=frozenset(supporting_reads),
                    compatible_transcript_ids=final_tx_ids,
                    source_alignments=None,
                    coverage_cache=coverage_array,
                )
            result.append(supported)
        return result


def _variant_sequences_with_shared_read_support(variant_sequences, read_groups):
    """Assign compatible original reads to each retained candidate.

    Compatibility is checked against the reads themselves, not transitively
    through other candidates: a short core can match two conflicting longer
    branches without either branch supporting the other. Overhangs outside a
    candidate are irrelevant, and coverage still uses each read's real bounds.
    No sequence is extended or invented here. Reads may support more than one
    alternative; downstream aggregation deduplicates them by set union.
    ``read_groups`` is the canonical partition from ``group_unique_sequences``.
    """
    variant_sequences = list(variant_sequences)
    if not variant_sequences:
        return []
    max_prefix = max(len(s.prefix) for s in variant_sequences)
    max_suffix = max(len(s.suffix) for s in variant_sequences)
    extra_reads = [read for seq in variant_sequences for read in seq.reads]
    index = SharedReadSupportIndex(
        read_groups, max_prefix=max_prefix, max_suffix=max_suffix, extra_reads=extra_reads
    )
    return index.support_sequences(variant_sequences)


def trim_variant_sequences(
        variant_sequences,
        min_variant_sequence_coverage,
        *,
        support_index=None,
        preserve_candidate_transcript_ids=False):
    """
    Trim VariantSequences to desired coverage and combine exact duplicates.

    Contained but non-identical candidates remain separate because reference
    compatibility is not known at this stage. Reconcile their compatible read
    support before applying the coverage threshold, including when overlap
    assembly is disabled.

    Parameters
    ----------
    variant_sequences : list of VariantSequence

    min_variant_sequence_coverage : int

    support_index : SharedReadSupportIndex, optional
        Prebuilt read-support index covering all candidate flanks. When omitted,
        one is built from the unique reads present in ``variant_sequences`` and
        reused across both pre-trim and post-trim support passes.

    preserve_candidate_transcript_ids : bool, optional
        When True, preserve each candidate's ``compatible_transcript_ids``
        (used when multiple disjoint transcript read groups share one index
        over the original read pool).

    Returns list of VariantSequence
    """
    variant_sequences = [sequence for sequence in variant_sequences if len(sequence) > 0]
    n_total = len(variant_sequences)
    if not variant_sequences:
        logger.debug(
            "Kept 0/%d variant sequences after read coverage trimming to >=%dx",
            n_total,
            min_variant_sequence_coverage)
        return []
    if support_index is None:
        reads = {read for sequence in variant_sequences for read in sequence.reads}
        # Compare once per distinct read sequence, rather than once per fragment
        # or per appearance in the retained assembly history. Keep this original
        # pool even when some of its candidates are discarded below.
        read_groups = group_unique_sequences(reads)
        max_prefix = max(len(s.prefix) for s in variant_sequences)
        max_suffix = max(len(s.suffix) for s in variant_sequences)
        support_index = SharedReadSupportIndex(
            read_groups, max_prefix=max_prefix, max_suffix=max_suffix
        )
    if not support_index.conflicting_segments:
        trimmed_variant_sequences = support_index.trim_and_support_no_conflicts(
            variant_sequences,
            min_variant_sequence_coverage,
            preserve_candidate_transcript_ids=preserve_candidate_transcript_ids,
        )
        logger.debug(
            "Kept %d/%d variant sequences after read coverage trimming to >=%dx",
            len(trimmed_variant_sequences),
            n_total,
            min_variant_sequence_coverage)
        return trimmed_variant_sequences

    variant_sequences = merge_identical_sequences(variant_sequences)
    variant_sequences = support_index.support_sequences(
        variant_sequences,
        preserve_candidate_transcript_ids=preserve_candidate_transcript_ids,
    )
    # when no base of a sequence has sufficient coverage, trim_by_coverage
    # returns an empty sequence which still carries every one of the original
    # reads. Those sequences are useless downstream and their read counts are
    # actively misleading when comparing candidates by support, so drop them
    # here rather than relying on a later validity check to sweep them up.
    trimmed_variant_sequences = []
    untrimmed_ids = set()
    for variant_sequence in variant_sequences:
        trimmed = variant_sequence.trim_by_coverage(min_variant_sequence_coverage)
        if len(trimmed) > 0:
            if len(trimmed) == len(variant_sequence):
                untrimmed_ids.add(id(trimmed))
            trimmed_variant_sequences.append(trimmed)
    trimmed_variant_sequences = merge_identical_sequences(
        trimmed_variant_sequences)
    # Removing low-coverage overhangs can make previously conflicting reads
    # compatible. Recount only shortened or newly merged candidates against the
    # original pool without repeating work for untrimmed candidates.
    recount_positions = [
        idx
        for idx, seq in enumerate(trimmed_variant_sequences)
        if id(seq) not in untrimmed_ids
    ]
    if recount_positions:
        recounted = support_index.support_sequences(
            [trimmed_variant_sequences[idx] for idx in recount_positions],
            preserve_candidate_transcript_ids=preserve_candidate_transcript_ids,
        )
        for idx, seq in zip(recount_positions, recounted):
            trimmed_variant_sequences[idx] = seq
    n_after_trimming = len(trimmed_variant_sequences)
    logger.debug(
        "Kept %d/%d variant sequences after read coverage trimming to >=%dx",
        n_after_trimming,
        n_total,
        min_variant_sequence_coverage)
    return trimmed_variant_sequences
