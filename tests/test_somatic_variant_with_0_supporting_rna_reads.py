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

from varcode import Variant


from isovar.read_collector import ReadCollector

from .genomes_for_testing import grch38
from .testing_helpers import load_bam
from .common import eq_ 


def test_somatic_variant_with_0_supporting_rna_reads():
    variant = Variant("6", 90411765, "G", "A", grch38)
    base_dir = "data/somatic-variant-with-0-supporting-rna-reads/"
    normal_reads = load_bam(base_dir + "normal.6.90411765.G.A.sorted.bam")
    tumor_reads = load_bam(base_dir + "tumor.6.90411765.G.A.sorted.bam")
    rna_reads = load_bam(base_dir + "rna.6.90411765.G.A.sorted.bam")
    read_creator = ReadCollector()
    normal_sample_variant_reads = read_creator.allele_reads_supporting_variant(
        variant=variant,
        alignment_file=normal_reads)
    eq_(len(normal_sample_variant_reads), 0)
    print(normal_sample_variant_reads)

    tumor_sample_variant_reads = read_creator.allele_reads_supporting_variant(
        variant=variant,
        alignment_file=tumor_reads)
    print(tumor_sample_variant_reads)
    eq_(len(tumor_sample_variant_reads), 5)

    rna_sample_variant_reads = read_creator.allele_reads_supporting_variant(
        variant=variant,
        alignment_file=rna_reads)
    print(rna_sample_variant_reads)
    eq_(len(rna_sample_variant_reads), 0)

if __name__ == "__main__":
    test_somatic_variant_with_0_supporting_rna_reads()
