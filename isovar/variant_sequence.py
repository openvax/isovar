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

import numpy as np


from .value_object import ValueObject
from .logging import get_logger

logger = get_logger(__name__)


class VariantSequence(ValueObject):
    """
    Representation of a cDNA sequence containing a mutation
    """

    __slots__ = [
        # nucleotides before a variant
        "prefix",
        # nucleotide sequence of a variant
        "alt",
        # nucleotides after a variant
        "suffix",
        # since we often want to look at prefix+alt+suffix, let's cache it
        "sequence",
        # reads which were used to determine this sequences
        "reads"
    ]

    def __init__(self, prefix, alt, suffix, reads):
        self.prefix = prefix
        self.alt = alt
        self.suffix = suffix
        self.sequence = prefix + alt + suffix
        self.reads = frozenset(reads)

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return (
            "VariantSequence("
            "prefix='%s', alt='%s', "
            "suffix='%s', sequence='%s', "
            "reads=<%d AlleleRead object%s>)") % (
                self.prefix,
                self.alt,
                self.suffix,
                self.sequence,
                len(self.reads),
                "s" if len(self.reads) == 0 or len(self.reads) > 1 else "")

    def __repr__(self):
        return str(self)

    @property
    def read_names(self):
        """
        Unique read names from which this VariantSequence was constructed.

        Returns set of str
        """
        return {r.name for r in self.reads}

    def contains(self, other):
        """
        Is the other VariantSequence a subsequence of this one?

        The two sequences must agree on the alt nucleotides, the prefix of the
        longer must contain the prefix of the shorter, and the suffix of the
        longer must contain the suffix of the shorter.
        """
        return (self.alt == other.alt and
                self.prefix.endswith(other.prefix) and
                self.suffix.startswith(other.suffix))

    def left_overlaps(self, other, min_overlap_size=1):
        """
        Does this VariantSequence overlap another on the left side?
        """

        if self.alt != other.alt:
            # allele must match!
            return False

        if len(other.prefix) > len(self.prefix):
            # only consider strings that overlap like:
            #   self: ppppAssss
            #   other:   ppAsssssss
            # which excludes cases where the other sequence has a longer
            # prefix
            return False
        elif len(other.suffix) < len(self.suffix):
            # similarly, we throw away cases where the other sequence is shorter
            # after the alt nucleotides than this sequence
            return False

        # is the other sequence a prefix of this sequence?
        # Example:
        # p1 a1 s1 = XXXXXXXX Y ZZZZZZ
        # p2 a2 s2 =       XX Y ZZZZZZZZZ
        # ...
        # then we can combine them into a longer sequence
        sequence_overlaps = (
            self.prefix.endswith(other.prefix) and
            other.suffix.startswith(self.suffix)
        )
        prefix_overlap_size = min(len(self.prefix), len(other.prefix))
        suffix_overlap_size = min(len(other.suffix), len(self.suffix))
        overlap_size = (
            prefix_overlap_size + suffix_overlap_size + len(self.alt))

        return sequence_overlaps and overlap_size >= min_overlap_size

    def add_reads(self, reads):
        """
        Create another VariantSequence with more supporting reads.
        """
        if len(reads) == 0:
            return self
        new_reads = self.reads.union(reads)
        if len(new_reads) > len(self.reads):
            return VariantSequence(
                prefix=self.prefix,
                alt=self.alt,
                suffix=self.suffix,
                reads=new_reads)
        else:
            return self

    def combine(self, other_sequence, min_overlap_size=1):
        """
        If this sequence is the prefix of another sequence, combine
        them into a single VariantSequence object. If the other sequence
        is contained in this one, then add its reads to this VariantSequence.
        Also tries to flip the order (e.g. this sequence is a suffix or
        this sequence is a subsequence). If sequences can't be combined
        then returns None.
        """
        if other_sequence.alt != self.alt:
            logger.warn(
                "Cannot combine %s and %s with mismatching alt sequences",
                self,
                other_sequence)
            return None
        elif self.contains(other_sequence):
            if len(other_sequence) >= min_overlap_size:
                return self.add_reads(other_sequence.reads)
            else:
                return None
        elif other_sequence.contains(self):
            if len(self) >= min_overlap_size:
                return other_sequence.add_reads(self.reads)
            else:
                return None
        elif self.left_overlaps(other_sequence, min_overlap_size=min_overlap_size):
            # If sequences are like AABC and ABCC
            return VariantSequence(
                prefix=self.prefix,
                alt=self.alt,
                suffix=other_sequence.suffix,
                reads=self.reads.union(other_sequence.reads))
        elif other_sequence.left_overlaps(self, min_overlap_size=min_overlap_size):
            return VariantSequence(
                prefix=other_sequence.prefix,
                alt=self.alt,
                suffix=self.suffix,
                reads=self.reads.union(other_sequence.reads))
        else:
            # sequences don't overlap
            return None

    def variant_indices(self):
        """
        When we combine prefix + alt + suffix into a single string,
        what are is base-0 index interval which gets us back the alt
        sequence? First returned index is inclusive, the second is exclusive.
        """
        variant_start_index = len(self.prefix)
        variant_len = len(self.alt)
        variant_end_index = variant_start_index + variant_len
        return variant_start_index, variant_end_index

    def coverage(self):
        """
        Returns NumPy array indicating number of reads covering each
        nucleotides of this sequence.
        """
        variant_start_index, variant_end_index = self.variant_indices()
        n_nucleotides = len(self)
        coverage_array = np.zeros(n_nucleotides, dtype="int32")
        for read in self.reads:
            coverage_array[
                max(0, variant_start_index - len(read.prefix)):
                min(n_nucleotides, variant_end_index + len(read.suffix))] += 1
        return coverage_array

    def min_coverage(self):
        """
        Minimum number of reads covering any base in the cDNA sequence

        Returns int
        """
        return np.min(self.coverage())

    def mean_coverage(self):
        """
        Average number of reads covering each base in the cDNA sequence.

        Returns float
        """
        return np.mean(self.coverage())

    def trim_by_coverage(self, min_reads):
        """
        Given the min number of reads overlapping each nucleotide of
        a variant sequence, trim this sequence by getting rid of positions
        which are overlapped by fewer reads than specified.
        """
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

