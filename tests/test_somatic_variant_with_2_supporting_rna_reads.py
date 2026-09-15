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
import pytest

from isovar.read_collector import ReadCollector

from .genomes_for_testing import grch38
from .testing_helpers import load_bam
from .common import eq_ 


@pytest.mark.parametrize("merge,expected_alt_reads", [(False, 2), (True, 0)])
def test_somatic_variant_with_2_supporting_rna_reads(merge, expected_alt_reads):
    variant = Variant("14", 105849746, "G", "A", grch38)
    base_dir = "data/somatic-variant-with-2-supporting-rna-reads/"
    normal_reads = load_bam(base_dir + "normal.14.105849746.G.A.no-alt.sorted.bam")
    tumor_reads = load_bam(base_dir + "tumor.14.105849746.G.A.many-alt.sorted.bam")
    rna_reads = load_bam(base_dir + "rna.14.105849746.G.A.2-alt.sorted.bam")
    read_creator = ReadCollector(merge_overlapping_fragments=merge)
    normal_sample_variant_reads = read_creator.allele_reads_supporting_variant(
        variant=variant,
        alignment_file=normal_reads)
    eq_(len(normal_sample_variant_reads), 0)
    print(normal_sample_variant_reads)

    tumor_sample_variant_reads = read_creator.allele_reads_supporting_variant(
        variant=variant,
        alignment_file=tumor_reads)
    print(tumor_sample_variant_reads)
    eq_(sum(r.source_read_count for r in tumor_sample_variant_reads), 8)

    rna_sample_variant_reads = read_creator.allele_reads_supporting_variant(
        variant=variant,
        alignment_file=rna_reads)
    print(rna_sample_variant_reads)
    # Both raw A calls have quality 12; their overlapping mates call G at
    # quality 41. The existing consensus rule selects G when merging.
    eq_(len(rna_sample_variant_reads), expected_alt_reads)
    # Arun went through the hassle of pulling out the exact read names
    # in IGV
    expected_variant_rna_read_names = {
        "K00193:50:H5NKVBBXX:5:2202:6421:24964",
        "K00193:50:H5NKVBBXX:5:2119:30908:1138",
    }
    for variant_read in rna_sample_variant_reads:
        assert variant_read.name in expected_variant_rna_read_names
