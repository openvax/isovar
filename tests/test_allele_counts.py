# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from isovar.dataframe_helpers import allele_counts_dataframe
from isovar.allele_read import AlleleRead
from isovar.read_evidence import ReadEvidence
from varcode import Variant
from .common import eq_


def test_allele_count_dataframe():
    variant = Variant("test_contig", 50, "C", "G")
    read_evidence = ReadEvidence(
            trimmed_base1_start=50,
            trimmed_ref="C",
            trimmed_alt="G",
            ref_reads=[
                AlleleRead(prefix="AAA", allele="C", suffix="TTT", name="C1"),
                AlleleRead(prefix="AAC", allele="C", suffix="TTA", name="C2"),
            ],
            alt_reads=[
                AlleleRead(prefix="AAA", allele="G", suffix="TTT", name="G1")
            ],
            other_reads=[])
    df = allele_counts_dataframe([(variant, read_evidence)])
    assert len(df) == 1, "Wrong number of rows in DataFrame: %s" % (df,)
    row = df.iloc[0]
    eq_(row.num_ref_reads, 2)
    eq_(row.num_alt_reads, 1)
    eq_(row.num_other_reads, 0)


if __name__ == "__main__":
    test_allele_count_dataframe()
