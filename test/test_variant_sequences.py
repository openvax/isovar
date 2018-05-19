# Copyright (c) 2016-2018. Mount Sinai School of Medicine
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

from nose.tools import eq_
from varcode import Variant
from isovar.variant_sequences import (
    reads_to_variant_sequences,
    VariantSequence
)
from isovar.variant_reads import reads_supporting_variant
from isovar.allele_reads import AlleleRead

from testing_helpers import load_bam


def test_sequence_counts_snv():
    samfile = load_bam("data/cancer-wgs-primary.chr12.bam")
    chromosome = "chr12"
    base1_location = 65857041
    ref = "G"
    alt = "C"
    variant = Variant(chromosome, base1_location, ref, alt)

    variant_reads = reads_supporting_variant(
        samfile=samfile,
        chromosome=chromosome,
        variant=variant)

    variant_sequences = reads_to_variant_sequences(
        variant=variant,
        reads=variant_reads,
        preferred_sequence_length=61)
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
