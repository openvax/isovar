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

from isovar.variant_reads import reads_supporting_variant
from pyensembl import ensembl_grch38

from varcode import Variant
from nose.tools import eq_

from testing_helpers import load_bam

def test_partition_variant_reads_snv():
    samfile = load_bam("data/cancer-wgs-primary.chr12.bam")
    chromosome = "chr12"
    base1_location = 65857041
    ref = "G"
    alt = "C"
    variant = Variant(
        contig=chromosome,
        start=base1_location,
        ref=ref,
        alt=alt,
        ensembl=ensembl_grch38)
    variant_reads = reads_supporting_variant(
        samfile=samfile,
        chromosome=chromosome,
        variant=variant)
    assert len(variant_reads) > 1
    for variant_read in variant_reads:
        eq_(variant_read.allele, alt)

def test_partition_variant_reads_deletion():
    samfile = load_bam("data/cancer-wgs-primary.chr12.bam")
    chromosome = "chr12"
    base1_location = 70091490
    ref = "TTGTAGATGCTGCCTCTCC"
    alt = ""
    variant = Variant(
        contig=chromosome,
        start=base1_location,
        ref=ref,
        alt=alt,
        ensembl=ensembl_grch38)
    variant_reads = reads_supporting_variant(
        samfile=samfile,
        chromosome=chromosome,
        variant=variant)
    assert len(variant_reads) > 1
    for variant_read in variant_reads:
        eq_(variant_read.allele, alt)

if __name__ == "__main__":
    test_partition_variant_reads_snv()
    test_partition_variant_reads_deletion()
