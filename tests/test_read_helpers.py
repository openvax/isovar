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

from isovar import ReadCollector
from isovar.allele_read_helpers import  group_unique_sequences

from varcode import Variant

from .testing_helpers import load_bam



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
        ensembl="hg38")
    read_collector = ReadCollector()
    variant_reads = read_collector.allele_reads_supporting_variant(
        alignment_file=samfile,
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
