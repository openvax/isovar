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

from isovar.nucleotide_counts import most_common_nucleotides
from isovar.variant_reads import gather_reads_for_single_variant
from pysam import AlignmentFile
from varcode import Variant
from nose.tools import eq_

def test_most_common_nucleotides_for_chr12_deletion():
    samfile = AlignmentFile("data/cancer-wgs-primary.chr12.bam")
    chromosome = "chr12"
    base1_location = 70091490
    ref = "TTGTAGATGCTGCCTCTCC"
    alt = ""
    variant = Variant(chromosome, base1_location, ref, alt)
    seq_parts = gather_reads_for_single_variant(
        samfile=samfile,
        chromosome=chromosome,
        variant=variant)
    consensus_sequence, chosen_counts, other_counts = most_common_nucleotides(
        seq_parts)
    eq_(len(chosen_counts), len(consensus_sequence))
    eq_(len(other_counts), len(consensus_sequence))
    eq_(other_counts.sum(), 0, "Didn't expect disagreement among reads")
    for variant_read in seq_parts:
        assert variant_read.prefix in consensus_sequence
        assert variant_read.suffix in consensus_sequence

if __name__ == "__main__":
    test_most_common_nucleotides_for_chr12_deletion()
