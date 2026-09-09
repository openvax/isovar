"""Differential tests against the pre-#227 exhaustive support matcher."""

import random
from itertools import product
from unittest.mock import patch

import numpy as np
import pytest

from isovar.allele_read import AlleleRead
from isovar.allele_read_helpers import group_unique_sequences
from isovar.variant_sequence import VariantSequence
from isovar.protein_sequence_creator import ProteinSequenceCreator
from isovar import variant_sequence_helpers as helpers


def exhaustive_support(sequences, read_groups):
    """Frozen simple oracle: only the actual read can support a candidate."""
    result = []
    for sequence in sequences:
        reads = set()
        for (prefix, alt, suffix), group in read_groups.items():
            overlap = min(len(prefix), len(sequence.prefix)) + len(alt) + min(len(suffix), len(sequence.suffix))
            if (alt == sequence.alt and overlap > 0
                    and (prefix.endswith(sequence.prefix) or sequence.prefix.endswith(prefix))
                    and (suffix.startswith(sequence.suffix) or sequence.suffix.startswith(suffix))):
                reads.update(group)
        result.append(VariantSequence(sequence.prefix, sequence.alt, sequence.suffix, reads))
    return result


def signature(sequences):
    return [(s.prefix, s.alt, s.suffix, s.reads, s.coverage().tolist()) for s in sequences]


def random_case(seed):
    rng = random.Random(seed)
    reads = []
    for i in range(60):
        # Repeated names are NOT interchangeable with distinct allele objects.
        read = AlleleRead(
            "".join(rng.choices("AACCGT", k=rng.randrange(18))),
            rng.choice(["", "A", "TGC"]),
            "".join(rng.choices("AACCGT", k=rng.randrange(18))),
            name=str(i % 17), source_read_count=rng.choice([1, 2, 3]))
        reads.append(read)
        reads.append(AlleleRead(read.prefix, read.allele, read.suffix, read.name + "_copy"))
    sequences = [VariantSequence(r.prefix, r.allele, r.suffix, [r]) for r in reads]
    for read in reads[::3]:
        prefix = read.prefix[rng.randrange(len(read.prefix) + 1):]
        suffix = read.suffix[:rng.randrange(len(read.suffix) + 1)]
        sequences.append(VariantSequence(prefix, read.allele, suffix, [read]))
    rng.shuffle(sequences)
    return reads, sequences


@pytest.mark.parametrize("seed", range(60))
def test_shared_support_exactly_matches_exhaustive(seed):
    reads, sequences = random_case(seed)
    groups = group_unique_sequences(reads)
    assert signature(helpers._variant_sequences_with_shared_read_support(sequences, groups)) == signature(
        exhaustive_support(sequences, groups))


@pytest.mark.parametrize("seed", range(30))
@pytest.mark.parametrize("minimum", [0, 1, 2, 4])
def test_both_trimming_passes_match_exhaustive(seed, minimum):
    _, sequences = random_case(seed)
    observed = helpers.trim_variant_sequences(sequences, minimum)
    with patch.object(helpers, "_variant_sequences_with_shared_read_support", exhaustive_support):
        expected = helpers.trim_variant_sequences(sequences, minimum)
    assert signature(observed) == signature(expected)


@pytest.mark.parametrize("alt", ["", "C", "AACGT"])
def test_extreme_flanks_overlap_bounds_and_duplicate_names(alt):
    flanks = ["", "A", "AA", "AAA", "CAA", "N", "\x00", "\U0010ffff", "\U0010ffffA", "AC" * 3000]
    reads = [AlleleRead(left, alt, right, "same", source_read_count=n % 3 + 1)
             for n, (left, right) in enumerate((left, right) for left in flanks for right in flanks)]
    sequences = [VariantSequence(left, alt, right, []) for left in flanks[:9] for right in flanks[:9]]
    groups = group_unique_sequences(reads)
    assert signature(helpers._variant_sequences_with_shared_read_support(sequences, groups)) == signature(
        exhaustive_support(sequences, groups))


def test_empty_inputs_and_candidates_without_support():
    sequence = VariantSequence("A", "", "T", [])
    assert helpers._variant_sequences_with_shared_read_support([], {}) == []
    observed, = helpers._variant_sequences_with_shared_read_support([sequence], {})
    assert observed is sequence


@pytest.mark.parametrize("seed", range(30))
def test_interval_coverage_matches_original_slice_updates(seed):
    reads, sequences = random_case(seed)
    sequences += [VariantSequence("A" * left, alt, "T" * right, reads)
                  for left, right in [(0, 0), (0, 25), (25, 0), (3, 7), (30, 30)]
                  for alt in ("", "C", "TGC")]
    for sequence in sequences:
        expected = np.zeros(len(sequence), dtype="int32")
        start, end = sequence.variant_indices()
        for read in sequence.reads:
            expected[max(0, start - len(read.prefix)):min(len(sequence), end + len(read.suffix))] += 1
        assert sequence.coverage().dtype == expected.dtype
        np.testing.assert_array_equal(sequence.coverage(), expected)
        assert sequence.coverage() is sequence.coverage()


def test_disabled_coverage_logging_does_not_format_arrays():
    class UnformattableCoverage(np.ndarray):
        def __str__(self):
            raise AssertionError("eager formatting")

    sequence = VariantSequence("AA", "C", "GG", [AlleleRead("AA", "C", "GG", "x")])
    sequence._coverage_cache = sequence.coverage().view(UnformattableCoverage)
    with patch("isovar.variant_sequence.logger.isEnabledFor", return_value=False):
        assert sequence.trim_by_coverage(1).reads == sequence.reads


@pytest.mark.parametrize("minimum", [0, 1, 2])
def test_trimmed_coverage_cache_is_exact_and_does_not_alias_parent(minimum):
    reads = [AlleleRead("AAAA", "C", "GGGG", "long"), AlleleRead("AA", "C", "GG", "short")]
    sequence = VariantSequence("AAAA", "C", "GGGG", reads)
    trimmed = sequence.trim_by_coverage(minimum)
    recalculated = VariantSequence(trimmed.prefix, trimmed.alt, trimmed.suffix, trimmed.reads)
    np.testing.assert_array_equal(trimmed.coverage(), recalculated.coverage())
    before = sequence.coverage().copy()
    trimmed.coverage()[:] = -100
    np.testing.assert_array_equal(sequence.coverage(), before)


def test_group_coverage_uses_original_object_counts_not_names_or_multiplicity():
    reads = [AlleleRead("X" * i + "AA", "C", "GG" + "Y" * i, "same", source_read_count=i + 1)
             for i in range(1, 100)]
    sequence = VariantSequence("AA", "C", "GG", [])
    result, = helpers._variant_sequences_with_shared_read_support([sequence], group_unique_sequences(reads))
    assert result.reads == frozenset(reads)
    assert result.coverage().tolist() == [99] * 5


@pytest.mark.parametrize("upstream", [False, True])
def test_index_avoids_cartesian_scan_for_distinct_branches(upstream, monkeypatch):
    """An operation-count regression, independent of CPU/CI timing noise."""
    reads = [AlleleRead("".join(bases) if upstream else "A", "C",
                        "A" if upstream else "".join(bases), str(i))
             for i, bases in enumerate(product("ACGT", repeat=5))]
    sequences = [VariantSequence(r.prefix, r.allele, r.suffix, [r]) for r in reads]
    visited = []

    class CountingEntries(list):
        def __getitem__(self, index):
            visited.append(index)
            return super().__getitem__(index)

    original = helpers._FlankIndex.__init__

    def instrument(self, *args, **kwargs):
        original(self, *args, **kwargs)
        self.entries = CountingEntries(self.entries)

    monkeypatch.setattr(helpers._FlankIndex, "__init__", instrument)
    assert helpers._variant_sequences_with_shared_read_support(sequences, group_unique_sequences(reads)) == sequences
    assert len(visited) == len(reads)  # 1,024 comparisons, not 1,048,576.


@pytest.mark.parametrize("preference", ["support", "balanced", "context"])
@pytest.mark.parametrize("assembly", [False, True])
@pytest.mark.parametrize("alt", ["", "C", "GATTACA"])
@pytest.mark.parametrize("minimum", [0, 2])
def test_adaptive_contexts_preserve_all_candidates_and_support(preference, assembly, alt, minimum):
    rng = random.Random(227)
    reads = []
    for i in range(45):
        left, right = "ACGT" * 25, "TGCA" * 25
        if i % 3:
            offset = rng.randrange(100)
            left = left[:offset] + "G" + left[offset + 1:]
        if i % 5:
            offset = rng.randrange(100)
            right = right[:offset] + "C" + right[offset + 1:]
        reads.append(AlleleRead(left[rng.randrange(80):], alt, right[:rng.randrange(1, 101)],
                                str(i % 31), source_read_count=1 + i % 2))
    creator = ProteinSequenceCreator(protein_sequence_preference=preference, variant_sequence_assembly=assembly,
                                     min_variant_sequence_coverage=minimum)
    observed = creator.variant_sequences_from_reads(None, reads)
    with patch.object(helpers, "_variant_sequences_with_shared_read_support", exhaustive_support):
        expected = creator.variant_sequences_from_reads(None, reads)
    assert signature(observed) == signature(expected)
