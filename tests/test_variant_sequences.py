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

from varcode import Variant

from isovar import (
    VariantSequenceCreator,
    VariantSequence
)
from isovar.allele_read import AlleleRead
from isovar.read_collector import ReadCollector
from isovar.variant_sequence_helpers import (
    filter_variant_sequences,
    filter_variant_sequences_by_length,
    trim_variant_sequences,
)

from .testing_helpers import load_bam
from .genomes_for_testing import grch38
from .common import eq_ 

def test_sequence_counts_snv():
    samfile = load_bam("data/cancer-wgs-primary.chr12.bam")
    chromosome = "chr12"
    base1_location = 65857041
    ref = "G"
    alt = "C"
    variant = Variant(chromosome, base1_location, ref, alt, "hg38")
    read_creator = ReadCollector()
    variant_reads = read_creator.allele_reads_supporting_variant(
        alignment_file=samfile,
        variant=variant)
    variant_sequence_creator = VariantSequenceCreator(preferred_sequence_length=61)
    variant_sequences = variant_sequence_creator.reads_to_variant_sequences(
        variant=variant,
        reads=variant_reads)
    assert len(variant_sequences) == 1
    for variant_sequence in variant_sequences:
        print(variant_sequence)
        eq_(variant_sequence.alt, alt)
        eq_(len(variant_sequence.prefix), 30)
        eq_(len(variant_sequence.suffix), 30)
        eq_(
            variant_sequence.prefix + variant_sequence.alt + variant_sequence.suffix,
            variant_sequence.sequence)


def test_variant_sequence_read_names():
    vs = VariantSequence(
        prefix="A",
        alt="C",
        suffix="T",
        reads=[
            AlleleRead(prefix="A", allele="C", suffix="T", name="1"),
            AlleleRead(prefix="A", allele="C", suffix="T", name="2")])
    eq_(vs.read_names, {"1", "2"})


def test_variant_sequence_contains():
    # AA|C|T
    vs_longer_prefix = VariantSequence(
        prefix="AA",
        alt="C",
        suffix="T",
        reads=[
            AlleleRead(
                prefix="AA", allele="C", suffix="T", name="longer_prefix")])
    # A|C|TT
    vs_longer_suffix = VariantSequence(
        prefix="A",
        alt="C",
        suffix="TT",
        reads=[
            AlleleRead(
                prefix="A", allele="C", suffix="TT", name="longer_suffix")])
    # A|C|T
    vs_short = VariantSequence(
        prefix="A",
        alt="C",
        suffix="T",
        reads=[
            AlleleRead(
                prefix="A", allele="C", suffix="T", name="short")])

    # two longer sequences contain the shorter subsequence
    assert vs_longer_prefix.contains(vs_short), \
        "Expected %s to contain %s" % (vs_longer_prefix, vs_short)
    assert vs_longer_suffix.contains(vs_short), \
        "Expected %s to contain %s" % (vs_longer_suffix, vs_short)
    # other pairs do not contain each other
    assert not vs_longer_prefix.contains(vs_longer_suffix), \
        "Expected %s to not contain %s" % (vs_longer_prefix, vs_longer_suffix)
    assert not vs_longer_suffix.contains(vs_longer_prefix), \
        "Expected %s to not contain %s" % (vs_longer_suffix, vs_longer_prefix)
    assert not vs_short.contains(vs_longer_prefix), \
        "Expected %s to not contain %s" % (vs_short, vs_longer_prefix)
    assert not vs_short.contains(vs_longer_suffix), \
        "Expected %s to not contain %s" % (vs_short, vs_longer_suffix)

    # Sequences above has 'C' allele whereas this one has 'G'
    # A|G|T
    vs_different_allele = VariantSequence(
        prefix="A",
        alt="G",
        suffix="T",
        reads=[
            AlleleRead(
                prefix="A", allele="G", suffix="T", name="short")])

    for vs in [vs_longer_suffix, vs_longer_prefix, vs_short]:
        assert not vs.contains(vs_different_allele), \
            "Expected %s to not contain %s" % (vs, vs_different_allele)
        assert not vs_different_allele.contains(vs), \
            "Expected %s to not contain %s" % (vs_different_allele, vs)


def test_variant_sequence_overlaps():
    # AAA|GG|TT
    vs_3A = VariantSequence(
        prefix="AAA",
        alt="GG",
        suffix="TT",
        reads=[
            AlleleRead(
                prefix="AAA", allele="GG", suffix="TT", name="1")])
    # AA|GG|TT
    vs_2A = VariantSequence(
        prefix="AA",
        alt="GG",
        suffix="TT",
        reads=[
            AlleleRead(
                prefix="AA", allele="GG", suffix="TT", name="1")])
    for min_overlap_size in [1, 2, 3, 4, 5, 6]:
        assert vs_3A.left_overlaps(vs_2A, min_overlap_size=min_overlap_size), \
            "Expected %s to overlap %s from left (min overlap size=%d)" % (
                vs_3A, vs_2A, min_overlap_size)

        assert not vs_2A.left_overlaps(vs_3A, min_overlap_size=min_overlap_size), \
            "Expected %s to not overlap %s from left (min overlap size=%d)" % (
                vs_2A, vs_3A, min_overlap_size)
    assert not vs_3A.left_overlaps(vs_2A, min_overlap_size=7), \
        "Unexpected overlap between %s and %s for min_overlap_size=7" % (
            vs_3A, vs_2A)


def test_variant_sequence_add_reads():
    vs = VariantSequence(prefix="A", alt="C", suffix="G", reads={"1"})
    # adding reads '2' and '3', sometimes multiple times
    vs_result = vs.add_reads("2").add_reads("1").add_reads("2").add_reads("3")
    expected = VariantSequence(prefix="A", alt="C", suffix="G", reads={"1", "2", "3"})
    eq_(vs_result, expected)


def test_variant_sequence_combine():
    vs1 = VariantSequence(prefix="A", alt="C", suffix="GG", reads={"1"})
    vs2 = VariantSequence(prefix="AA", alt="C", suffix="GG", reads={"2"})
    vs_result_1_to_2 = vs1.combine(vs2)
    expected = VariantSequence(prefix="AA", alt="C", suffix="GG", reads={"1", "2"})
    eq_(vs_result_1_to_2, expected)

    # shouldn't matter which sequence is first as an argument to the combine
    # function
    vs_result_2_to_1 = vs2.combine(vs1)
    eq_(vs_result_2_to_1, expected)


def test_variant_sequence_trim_by_coverage():
    reads = [
        AlleleRead(
            prefix="AA", allele="C", suffix="T", name="1"),
        AlleleRead(
            prefix="A", allele="C", suffix="T", name="2")
    ]
    vs = VariantSequence(
        prefix="AA",
        alt="C",
        suffix="T",
        reads=reads)
    # every nucleotide is spanned by one read
    eq_(vs.trim_by_coverage(1), vs)

    vs_expected_trim_by_2 = VariantSequence(
        prefix="A",
        alt="C",
        suffix="T",
        reads=reads)
    eq_(vs.trim_by_coverage(2), vs_expected_trim_by_2)


def test_variant_sequence_min_coverage():
    # 1: AA|C|TT
    # 2: AA|C|T
    # 3:  A|C|TT
    reads = [
        AlleleRead(
            prefix="AA", allele="C", suffix="TT", name="1"),
        AlleleRead(
            prefix="AA", allele="C", suffix="T", name="2"),
        AlleleRead(
            prefix="A", allele="C", suffix="TT", name="3")
    ]
    vs = VariantSequence(
        prefix="AA",
        alt="C",
        suffix="TT",
        reads=reads)
    eq_(vs.min_coverage(), 2)


def test_variant_sequence_mean_coverage():
    # 1: AA|C|TT
    # 2: AA|C|T
    # 3:  A|C|TT
    reads = [
        AlleleRead(
            prefix="AA", allele="C", suffix="TT", name="1"),
        AlleleRead(
            prefix="AA", allele="C", suffix="T", name="2"),
        AlleleRead(
            prefix="A", allele="C", suffix="TT", name="3")
    ]
    vs = VariantSequence(
        prefix="AA",
        alt="C",
        suffix="TT",
        reads=reads)
    # count the number of times a nucleotide in the sequences above
    # is contained in a read
    expected_mean_coverage = (2 + 3 + 3 + 3 + 2) / 5
    eq_(vs.mean_coverage(), expected_mean_coverage)


def test_variant_sequence_len():
    vs = VariantSequence(
        prefix="AA",
        alt="C",
        suffix="TT",
        reads=[])
    eq_(len(vs), 5)


def test_variant_sequence_coverage_is_cached():
    reads = [
        AlleleRead(prefix="AA", allele="C", suffix="TT", name="r1"),
        AlleleRead(prefix="A", allele="C", suffix="TT", name="r2"),
    ]
    vs = VariantSequence(prefix="AA", alt="C", suffix="TT", reads=reads)
    first_call = vs.coverage()
    second_call = vs.coverage()
    assert first_call is second_call, \
        "coverage() should return the same cached array on repeated calls"


def test_variant_sequence_coverage_cache_excluded_from_equality():
    reads = [AlleleRead(prefix="AA", allele="C", suffix="TT", name="r1")]
    vs1 = VariantSequence(prefix="AA", alt="C", suffix="TT", reads=reads)
    vs2 = VariantSequence(prefix="AA", alt="C", suffix="TT", reads=reads)
    vs1.coverage()
    eq_(vs1, vs2, "Cached coverage should not affect equality")
    eq_(hash(vs1), hash(vs2), "Cached coverage should not affect hash")


def _variant_sequence(prefix, alt, suffix, n_fragments, name_prefix):
    """
    Build a VariantSequence fully spanned by n_fragments distinct reads.
    """
    reads = [
        AlleleRead(
            prefix=prefix,
            allele=alt,
            suffix=suffix,
            name="%s_%d" % (name_prefix, i))
        for i in range(n_fragments)
    ]
    return VariantSequence(prefix=prefix, alt=alt, suffix=suffix, reads=reads)


def test_filter_by_length_keeps_shorter_sequence_with_better_support():
    # a long sequence assembled from few reads should not evict a shorter
    # sequence supported by many more fragments, since length at this stage
    # only reflects how far a chain of reads happened to extend
    well_supported = _variant_sequence("A" * 20, "C", "T" * 20, 7, "concordant")
    barely_supported = _variant_sequence("G" * 30, "C", "T" * 30, 2, "outlier")
    kept = filter_variant_sequences_by_length(
        [well_supported, barely_supported],
        preferred_sequence_length=61)
    eq_(set(kept), {well_supported, barely_supported})


def test_filter_by_length_defers_when_longer_sequence_has_better_support():
    # The longer sequence may fail reference matching, which cannot be known
    # at this stage, so even a less-supported shorter candidate must survive.
    short_sequence = _variant_sequence("A" * 20, "C", "T" * 20, 2, "short")
    long_sequence = _variant_sequence("G" * 30, "C", "T" * 30, 7, "long")
    kept = filter_variant_sequences_by_length(
        [short_sequence, long_sequence],
        preferred_sequence_length=61)
    eq_(kept, [short_sequence, long_sequence])


def test_filter_by_length_defers_when_sequences_have_equal_fragment_support():
    # Fragment equality is not enough to choose: raw read counts, reference
    # mismatches, and aggregate support are only available downstream.
    short_sequence = _variant_sequence("A" * 20, "C", "T" * 20, 5, "short")
    long_sequence = _variant_sequence("G" * 30, "C", "T" * 30, 5, "long")
    kept = filter_variant_sequences_by_length(
        [short_sequence, long_sequence],
        preferred_sequence_length=61)
    eq_(kept, [short_sequence, long_sequence])


def test_filter_by_length_defers_raw_read_count_tie_break():
    # Merged mates have one fragment name but preserve both raw reads through
    # source_read_count. The final ProteinSequence sort ranks this count before
    # length, so a fragment-only early filter would choose incorrectly.
    short_reads = [
        AlleleRead(
            prefix="A" * 20,
            allele="C",
            suffix="T" * 20,
            name="short_%d" % i,
            source_read_count=2)
        for i in range(3)
    ]
    short_sequence = VariantSequence(
        prefix="A" * 20,
        alt="C",
        suffix="T" * 20,
        reads=short_reads)
    long_sequence = _variant_sequence("G" * 30, "C", "T" * 30, 3, "long")

    kept = filter_variant_sequences_by_length(
        [short_sequence, long_sequence],
        preferred_sequence_length=61)

    eq_(kept, [short_sequence, long_sequence])
    eq_(sum(r.source_read_count for r in short_sequence.reads), 6)
    eq_(sum(r.source_read_count for r in long_sequence.reads), 3)


def test_filter_by_length_defers_aggregate_protein_support():
    # Distinct cDNA candidates can translate to the same amino-acid sequence.
    # Their disjoint reads are unioned only after translation, so neither can
    # be discarded based on its individual support here.
    short_1 = _variant_sequence("A" * 20, "C", "T" * 20, 2, "short_1")
    short_2 = _variant_sequence("G" * 20, "C", "T" * 20, 2, "short_2")
    long_sequence = _variant_sequence("C" * 30, "C", "T" * 30, 3, "long")

    kept = filter_variant_sequences_by_length(
        [short_1, short_2, long_sequence],
        preferred_sequence_length=61)

    eq_(kept, [short_1, short_2, long_sequence])


def test_filter_by_length_ignores_read_counts_of_empty_sequences():
    # trim_by_coverage returns an empty sequence which keeps all of the
    # original reads, so an empty sequence can claim more fragments than any
    # real one. It must never be rescued on those grounds.
    real_sequence = _variant_sequence("G" * 30, "C", "T" * 30, 2, "real")
    empty_sequence = VariantSequence(
        prefix="",
        alt="",
        suffix="",
        reads=[
            AlleleRead(prefix="A", allele="C", suffix="T", name="empty_%d" % i)
            for i in range(14)
        ])
    kept = filter_variant_sequences_by_length(
        [real_sequence, empty_sequence],
        preferred_sequence_length=61)
    eq_(kept, [real_sequence])


def test_filter_by_length_drops_empty_sequence_when_it_is_the_only_candidate():
    degenerate = VariantSequence(
        prefix="",
        alt="",
        suffix="",
        reads=[AlleleRead(prefix="A", allele="", suffix="T", name="empty")])

    kept = filter_variant_sequences_by_length(
        [degenerate],
        preferred_sequence_length=61)

    eq_(kept, [])


def test_filter_by_length_retains_nonempty_candidates_for_every_input_order():
    short = _variant_sequence("A" * 10, "C", "T" * 10, 2, "short")
    long = _variant_sequence("G" * 30, "C", "T" * 30, 5, "long")
    deletion = _variant_sequence("C" * 20, "", "A" * 20, 3, "deletion")
    degenerate = VariantSequence(
        prefix="",
        alt="",
        suffix="",
        reads=[AlleleRead(prefix="A", allele="", suffix="T", name="empty")])

    for preferred_sequence_length in (0, 21, 41, 61, 1000):
        for ordered_candidates in permutations(
                [short, long, deletion, degenerate]):
            kept = filter_variant_sequences_by_length(
                list(ordered_candidates),
                preferred_sequence_length=preferred_sequence_length)
            expected = [s for s in ordered_candidates if len(s) > 0]
            eq_(kept, expected)


def test_filter_by_length_keeps_deletion_sequences():
    # deletions legitimately have an empty alt allele, so degenerate
    # sequences must be identified by total length rather than by alt
    deletion = _variant_sequence("A" * 20, "", "T" * 20, 7, "deletion")
    longer_deletion = _variant_sequence("G" * 30, "", "T" * 30, 2, "long_deletion")
    kept = filter_variant_sequences_by_length(
        [deletion, longer_deletion],
        preferred_sequence_length=61)
    eq_(set(kept), {deletion, longer_deletion})


def test_trim_variant_sequences_drops_sequences_without_coverage():
    # one fragment can't meet a 2x coverage threshold, so trimming leaves an
    # empty sequence which should be discarded rather than passed along
    uncovered = _variant_sequence("A" * 20, "C", "T" * 20, 1, "uncovered")
    covered = _variant_sequence("G" * 20, "C", "T" * 20, 3, "covered")
    trimmed = trim_variant_sequences(
        [uncovered, covered],
        min_variant_sequence_coverage=2)
    eq_(trimmed, [covered])
    assert all(len(s) > 0 for s in trimmed), \
        "trim_variant_sequences should not return empty sequences"


# Deletions legitimately carry an empty alt allele (see read_collector.py,
# where is_deletion is defined as len(trimmed_alt) == 0). It is therefore
# always wrong to identify a useless "degenerate" sequence by an empty alt;
# the test has to be total sequence length. The tests below pin that down at
# every stage which discards sequences, since a guard written against alt
# would silently drop every deletion in a run while looking perfectly
# reasonable in review.

def test_trim_variant_sequences_keeps_deletion_sequences():
    deletion = _variant_sequence("A" * 20, "", "T" * 20, 3, "deletion")
    trimmed = trim_variant_sequences([deletion], min_variant_sequence_coverage=2)
    eq_(trimmed, [deletion])
    eq_(trimmed[0].alt, "")


def test_trim_variant_sequences_distinguishes_deletion_from_degenerate():
    # both have an empty alt; only the one whose bases lack coverage should go
    covered_deletion = _variant_sequence("A" * 20, "", "T" * 20, 3, "covered")
    uncovered_deletion = _variant_sequence("G" * 20, "", "C" * 20, 1, "uncovered")
    trimmed = trim_variant_sequences(
        [covered_deletion, uncovered_deletion],
        min_variant_sequence_coverage=2)
    eq_(trimmed, [covered_deletion])


def test_filter_by_length_distinguishes_deletion_from_degenerate_sequence():
    # a deletion and a degenerate sequence both have an empty alt, so only
    # total length can tell them apart. The deletion is shorter and better
    # supported, so it should be rescued; the empty one never should be,
    # despite claiming the most fragments of all.
    long_sequence = _variant_sequence("G" * 30, "", "T" * 30, 2, "long")
    short_deletion = _variant_sequence("A" * 20, "", "T" * 20, 7, "deletion")
    degenerate = VariantSequence(
        prefix="",
        alt="",
        suffix="",
        reads=[
            AlleleRead(prefix="A", allele="", suffix="T", name="degenerate_%d" % i)
            for i in range(14)
        ])
    kept = filter_variant_sequences_by_length(
        [long_sequence, short_deletion, degenerate],
        preferred_sequence_length=61)
    eq_(set(kept), {long_sequence, short_deletion})


def test_filter_variant_sequences_keeps_well_supported_deletion():
    # Coverage trimming and defensive validity filtering together, on a
    # deletion whose best supported sequence is not its longest.
    long_sequence = _variant_sequence("G" * 30, "", "T" * 30, 2, "long")
    short_deletion = _variant_sequence("A" * 20, "", "T" * 20, 7, "deletion")
    kept = filter_variant_sequences(
        [long_sequence, short_deletion],
        preferred_sequence_length=61,
        min_variant_sequence_coverage=2)
    eq_(set(kept), {long_sequence, short_deletion})
    assert all(s.alt == "" for s in kept)


def test_reads_to_variant_sequences_keeps_deletion():
    # end to end through the creator: a deletion supported by well covered
    # reads must survive assembly, coverage trimming, and validity filtering
    variant = Variant("chr12", 65857041, "G", "", grch38)
    prefix, suffix = "A" * 30, "T" * 30
    reads = [
        AlleleRead(prefix=prefix, allele="", suffix=suffix, name="deletion_%d" % i)
        for i in range(4)
    ]
    creator = VariantSequenceCreator(preferred_sequence_length=61)
    variant_sequences = creator.reads_to_variant_sequences(
        variant=variant, reads=reads)
    eq_(len(variant_sequences), 1)
    eq_(variant_sequences[0].alt, "")
    eq_(variant_sequences[0].sequence, prefix + suffix)
