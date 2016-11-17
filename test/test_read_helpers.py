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

from isovar.read_helpers import group_unique_sequences
from isovar.variant_reads import reads_supporting_variant
from varcode import Variant
from pyensembl import ensembl_grch38

from testing_helpers import load_bam

def test_group_unique_sequences():
    samfile = load_bam("data/cancer-wgs-primary.chr12.bam")
    chromosome = "chr12"
    base1_location = 65857041
    ref = "G"
    alt = "C"
    variant = Variant(
        contig=chromosome,
        start=base1_location,
        ref=ref, alt=alt,
        ensembl=ensembl_grch38)
    variant_reads = reads_supporting_variant(
        samfile=samfile,
        chromosome=chromosome,
        variant=variant)
    print("%d variant reads: %s" % (
        len(variant_reads), variant_reads))
    groups = group_unique_sequences(
        variant_reads,
        max_prefix_size=30,
        max_suffix_size=30)
    print("%d unique sequences: %s" % (
        len(groups), groups))
    # there are some redundant reads, so we expect that the number of
    # unique entries should be less than the total read partitions
    assert len(variant_reads) > len(groups)
