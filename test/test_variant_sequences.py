# Copyright (c) 2016. Mount Sinai School of Medicine
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


"""
class VariantSequence(VariantSequenceBase):
    def __new__(cls, prefix, alt, suffix, reads):
        # construct sequence from prefix + alt + suffix
        return VariantSequenceBase.__new__(
            cls,
            prefix=prefix,
            alt=alt,
            suffix=suffix,
            sequence=prefix + alt + suffix,
            reads=frozenset(reads))

    def __len__(self):
        return len(self.sequence)

    @property
    def read_names(self):
        return {r.name for r in self.reads}

    def contains(self, other):
        return (self.alt == other.alt and
                self.prefix.endswith(other.prefix) and
                self.suffix.startswith(other.suffix))

    def left_overlaps(self, other, min_overlap_size=1):
        if self.alt != other.alt:
            # allele must match!
            return False

        if len(other.prefix) > len(self.prefix):
            # only consider strings that overlap like:
            #   variant_sequence1: ppppAssss
            #   variant_sequence2:   ppAsssssss
            # which excludes cases where variant_sequence2 has a longer
            # prefix
            return False
        elif len(other.suffix) < len(self.suffix):
            # similarly, we throw cases where variant_sequence2 is shorter
            # after the alt nucleotides than variant_sequence1
            return False

        # is the candidate sequence is a prefix of the accepted?
        # Example:
        # p1 a1 s1 = XXXXXXXX Y ZZZZZZ
        # p2 a2 s2 =       XX Y ZZZZZZZZZ
        # ...
        # then we can combine them into a longer sequence
        sequence_overlaps = (
            self.prefix.endswith(other.prefix) and
            other.suffix.startswith(self.suffix)
        )
        prefix_overlap_size = len(self.prefix) - len(other.prefix)
        suffix_overlap_size = len(other.suffix) - len(self.suffix)
        overlap_size = prefix_overlap_size + suffix_overlap_size + len(self.alt)
        return sequence_overlaps and overlap_size >= min_overlap_size

    def add_reads(self, reads):
        return VariantSequence(
            prefix=self.prefix,
            alt=self.alt,
            suffix=self.suffix,
            reads=self.reads.union(reads))

    def combine(self, other_sequence):
        if other_sequence.alt != self.alt:
            raise ValueError(
                "Cannot combine %s and %s with mismatching alt sequences" % (
                    self,
                    other_sequence))
        elif self.left_overlaps(other_sequence):
            # If sequences are like AABC and ABCC
            return VariantSequence(
                prefix=self.prefix,
                alt=self.alt,
                suffix=other_sequence.suffix,
                reads=self.reads.union(other_sequence.reads))
        elif other_sequence.left_overlaps(self):
            return other_sequence.combine(self)
        elif self.contains(other_sequence):
            return self.add_reads(other_sequence.reads)
        elif other_sequence.contains(self):
            return other_sequence.add_reads(self.reads)
        else:
            raise ValueError("%s does not overlap with %s" % (self, other_sequence))

    def variant_indices(self):
        variant_start_index = len(self.prefix)
        variant_len = len(self.alt)
        variant_end_index = variant_start_index + variant_len
        return variant_start_index, variant_end_index

    def coverage(self):
        variant_start_index, variant_end_index = self.variant_indices()
        n_nucleotides = len(self)
        coverage_array = np.zeros(n_nucleotides, dtype="int32")
        for read in self.reads:
            coverage_array[
                max(0, variant_start_index - len(read.prefix)):
                min(n_nucleotides, variant_end_index + len(read.suffix))] += 1
        return coverage_array

    def min_coverage(self):
        return np.min(self.coverage())

    def mean_coverage(self):
        return np.mean(self.coverage())

    def trim_by_coverage(self, min_reads):
        read_count_array = self.coverage()
        logger.info("Coverage: %s (len=%d)" % (
            read_count_array, len(read_count_array)))
        sufficient_coverage_mask = read_count_array >= min_reads
        sufficient_coverage_indices = np.argwhere(sufficient_coverage_mask)
        if len(sufficient_coverage_indices) == 0:
            logger.debug("No bases in %s have coverage >= %d" % (self, min_reads))
            return VariantSequence(prefix="", alt="", suffix="", reads=self.reads)
        variant_start_index, variant_end_index = self.variant_indices()
        # assuming that coverage drops off monotonically away from
        # variant nucleotides
        first_covered_index = sufficient_coverage_indices.min()
        last_covered_index = sufficient_coverage_indices.max()
        # adding 1 to last_covered_index since it's an inclusive index
        # whereas variant_end_index is the end of a half-open interval
        if (first_covered_index > variant_start_index or
                last_covered_index + 1 < variant_end_index):
            # Example:
            #   Nucleotide sequence:
            #       ACCCTTTT|AA|GGCGCGCC
            #   Coverage:
            #       12222333|44|33333211
            # Then the mask for bases covered >= 4x would be:
            #       ________|**|________
            # with indices:
            #       first_covered_index = 9
            #       last_covered_index = 10
            #       variant_start_index = 9
            #       variant_end_index = 11
            logger.debug("Some variant bases in %s don't have coverage >= %d" % (
                self, min_reads))
            return VariantSequence(prefix="", alt="", suffix="", reads=self.reads)
        return VariantSequence(
            prefix=self.prefix[first_covered_index:],
            alt=self.alt,
            suffix=self.suffix[:last_covered_index - variant_end_index + 1],
            reads=self.reads)

    def sort_key_decreasing_read_count(self):
        return -len(self.reads)
"""
