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

from pyensembl import ensembl_grch38

from varcode import Variant


from isovar import ReadCollector

from .common import eq_
from .testing_helpers import load_bam

def test_partition_variant_reads_snv():
    alignment_file = load_bam("data/cancer-wgs-primary.chr12.bam")
    chromosome = "12"
    base1_location = 65857041
    ref = "G"
    alt = "C"
    variant = Variant(
        contig=chromosome,
        start=base1_location,
        ref=ref,
        alt=alt,
        ensembl=ensembl_grch38)
    read_collector = ReadCollector()
    read_evidence = read_collector.read_evidence_for_variant(
        alignment_file=alignment_file,
        variant=variant)
    alt_reads = read_evidence.alt_reads
    assert len(alt_reads) > 1
    for variant_read in alt_reads:
        eq_(variant_read.allele, alt)


def test_partition_variant_reads_deletion():
    alignment_file = load_bam("data/cancer-wgs-primary.chr12.bam")
    chromosome = "12"
    base1_location = 70091490
    ref = "TTGTAGATGCTGCCTCTCC"
    alt = ""
    variant = Variant(
        contig=chromosome,
        start=base1_location,
        ref=ref,
        alt=alt,
        ensembl=ensembl_grch38)
    read_collector = ReadCollector()
    read_evidence = read_collector.read_evidence_for_variant(
        alignment_file=alignment_file,
        variant=variant)
    assert len(read_evidence.alt_reads) > 1
    for variant_read in read_evidence.alt_reads:
        eq_(variant_read.allele, alt)

if __name__ == "__main__":
    test_partition_variant_reads_snv()
    test_partition_variant_reads_deletion()
