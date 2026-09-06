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

from itertools import permutations
import random

import pytest

from isovar.allele_read import AlleleRead
from isovar.variant_sequence import VariantSequence
from isovar.variant_sequence_creator import VariantSequenceCreator
from isovar.variant_sequence_helpers import trim_variant_sequences


def _sequence(prefix, alt, suffix, name):
    return VariantSequence(prefix, alt, suffix, [
        AlleleRead(prefix, alt, suffix, name=name)])


@pytest.mark.parametrize("alt", ["C", "", "CTG"])
def test_support_is_shared_without_transferring_conflicting_reads(alt):
    left = _sequence("AAA", alt, "GGG", "left")
    right = _sequence("TAA", alt, "GGG", "right")
    core = _sequence("AA", alt, "GG", "core")
    suffix_conflict = _sequence("AAA", alt, "TGG", "suffix_conflict")
    allele_conflict = _sequence("AAA", "T", "GGG", "allele_conflict")
    # The core can carry reads with incompatible overhangs outside its own
    # bounds. Those reads must not leak into either longer branch by way of it.
    core = core.add_reads(left.reads | right.reads)
    expected = {
        ("AAA", alt, "GGG"): {"left", "core"},
        ("TAA", alt, "GGG"): {"right", "core"},
        ("AA", alt, "GG"): {"left", "right", "core"},
        ("AAA", alt, "TGG"): {"suffix_conflict"},
        ("AAA", "T", "GGG"): {"allele_conflict"},
    }
    for ordered in permutations([left, right, core, suffix_conflict, allele_conflict]):
        observed = trim_variant_sequences(list(ordered), min_variant_sequence_coverage=1)
        assert {(s.prefix, s.alt, s.suffix): s.read_names for s in observed} == expected


def test_identical_candidates_pool_support_before_coverage_filtering():
    first = _sequence("AA", "C", "GG", "first")
    second = _sequence("AA", "C", "GG", "second")
    result, = trim_variant_sequences([first, second, first], min_variant_sequence_coverage=2)
    assert result.read_names == {"first", "second"}
    assert result.coverage().tolist() == [2] * 5


def test_shared_support_respects_actual_read_bounds_and_invalidates_coverage_cache():
    left = _sequence("AAAA", "C", "GG", "left")
    right = _sequence("AA", "C", "GGGG", "right")
    assert left.min_coverage() == right.min_coverage() == 1

    # No new flanking sequence is assembled here: only their common 2x core
    # survives, even though both candidates have two compatible reads.
    result, = trim_variant_sequences([left, right], min_variant_sequence_coverage=2)
    assert (result.prefix, result.alt, result.suffix) == ("AA", "C", "GG")
    assert result.read_names == {"left", "right"}
    assert result.coverage().tolist() == [2] * 5
    assert left.min_coverage() == right.min_coverage() == 1


def test_newly_trimmed_core_recovers_reads_from_a_discarded_conflicting_candidate():
    left = _sequence("AAAA", "C", "GG", "left")
    right = _sequence("AA", "C", "GGGG", "right")
    # Conflicts with the left prefix and right suffix, but not with the 2x
    # core produced by trimming them. Its own candidate is discarded at 1x.
    outside_conflicts = _sequence("TAA", "C", "GGT", "outside_conflicts")
    for ordered in permutations([left, right, outside_conflicts]):
        result, = trim_variant_sequences(list(ordered), min_variant_sequence_coverage=2)
        assert (result.prefix, result.alt, result.suffix) == ("AA", "C", "GG")
        assert result.read_names == {"left", "right", "outside_conflicts"}
        assert result.coverage().tolist() == [3] * 5


def test_empty_candidates_neither_survive_nor_donate_support():
    real = _sequence("AA", "", "GG", "real")
    empty = VariantSequence("", "", "", [
        AlleleRead("AA", "", "GG", name="empty_%d" % i) for i in range(4)])
    assert trim_variant_sequences([empty], min_variant_sequence_coverage=1) == []
    assert trim_variant_sequences([real, empty], min_variant_sequence_coverage=2) == []


def test_identical_full_strings_at_different_variant_positions_do_not_share_support():
    first = _sequence("T", "C", "CG", "first")
    second = _sequence("TC", "C", "G", "second")
    assert first.sequence == second.sequence
    assert trim_variant_sequences([first, second], min_variant_sequence_coverage=2) == []


def test_deletion_reads_without_any_shared_bases_do_not_inflate_support():
    left = _sequence("AA", "", "", "left")
    right = _sequence("", "", "GG", "right")
    observed = trim_variant_sequences([left, right], min_variant_sequence_coverage=1)
    assert set(observed) == {left, right}


@pytest.mark.parametrize("assembly", [False, True])
def test_creator_tied_support_is_order_independent_without_requiring_a_merge(assembly):
    candidates = [_sequence("AAA", "C", suffix, "%s_%d" % (suffix, i))
                  for suffix in ("GGG", "TTT") for i in range(2)]
    reads = [next(iter(candidate.reads)) for candidate in candidates]
    creator = VariantSequenceCreator(variant_sequence_assembly=assembly)
    expected = None
    for ordered in permutations(reads):
        # This method uses the reads, not the variant metadata.
        observed = creator.reads_to_variant_sequences(None, iter(ordered))
        assert len(observed) == 2
        if expected is None:
            expected = observed
        assert observed == expected


def test_shared_support_matches_independent_base_coordinate_oracle():
    """Exercise anchored matching against a per-base oracle, not string tests."""
    rng = random.Random(209)

    def coordinates(prefix, alt, suffix):
        return dict(enumerate(prefix + alt + suffix, start=-len(prefix)))

    for _ in range(20):
        reads = [AlleleRead(
            prefix="".join(rng.choices("ACGT", k=rng.randrange(7))),
            allele=rng.choice(["C", "", "CTG"]),
            suffix="".join(rng.choices("ACGT", k=rng.randrange(7))),
            name=str(i)) for i in range(30)]
        reads = [read for read in reads if len(read)]
        candidates = [VariantSequence(r.prefix, r.allele, r.suffix, [r]) for r in reads]
        for read in reads:
            prefix = read.prefix[rng.randrange(len(read.prefix) + 1):]
            suffix = read.suffix[:rng.randrange(len(read.suffix) + 1)]
            if prefix or read.allele or suffix:
                candidates.append(VariantSequence(prefix, read.allele, suffix, [read]))

        expected = {}
        for candidate in candidates:
            bases = coordinates(candidate.prefix, candidate.alt, candidate.suffix)
            supporting_reads = set()
            coverage = dict.fromkeys(bases, 0)
            for read in reads:
                read_bases = coordinates(read.prefix, read.allele, read.suffix)
                common_positions = bases.keys() & read_bases.keys()
                if (read.allele == candidate.alt and common_positions
                        and all(bases[p] == read_bases[p] for p in common_positions)):
                    supporting_reads.add(read)
                    for position in common_positions:
                        coverage[position] += 1
            expected[(candidate.prefix, candidate.alt, candidate.suffix)] = (
                supporting_reads, list(coverage.values()))

        for ordered in (candidates, list(reversed(candidates))):
            results = trim_variant_sequences(ordered, min_variant_sequence_coverage=1)
            assert {(s.prefix, s.alt, s.suffix): (s.reads, s.coverage().tolist())
                    for s in results} == expected
