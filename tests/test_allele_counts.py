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
from isovar.read_collector import ReadCollector
from varcode import Variant

from .mock_objects import MockAlignmentFile, make_pysam_read
from .common import eq_


def test_allele_count_dataframe():
    variant = Variant("test_contig", 50, "C", "G")
    read_evidence = ReadEvidence(
            trimmed_base1_start=50,
            trimmed_ref="C",
            trimmed_alt="G",
            ref_reads=[
                AlleleRead(prefix="AAA", allele="C", suffix="TTT", name="C1"),
                AlleleRead(prefix="AAC", allele="C", suffix="TTA", name="C1"),
            ],
            alt_reads=[
                AlleleRead(prefix="AAA", allele="G", suffix="TTT", name="G1"),
                AlleleRead(prefix="AAT", allele="G", suffix="TTC", name="G1"),
            ],
            other_reads=[
                AlleleRead(prefix="CCA", allele="T", suffix="GGA", name="T1"),
                AlleleRead(prefix="CCG", allele="T", suffix="GGT", name="T2"),
                AlleleRead(prefix="CCT", allele="T", suffix="GGC", name="T2"),
            ])
    df = allele_counts_dataframe([(variant, read_evidence)])
    assert len(df) == 1, "Wrong number of rows in DataFrame: %s" % (df,)
    row = df.iloc[0]
    eq_(row.num_ref_reads, 2)
    eq_(row.num_alt_reads, 2)
    eq_(row.num_other_reads, 3)
    eq_(row.num_ref_fragments, 1)
    eq_(row.num_alt_fragments, 1)
    eq_(row.num_other_fragments, 2)


def test_allele_count_dataframe_preserves_raw_read_count_after_fragment_merge():
    """
    Overlapping mates from one fragment should be merged for assembly while the
    read-count columns still report both raw reads.

    Regression test for GitHub issue #59.
    """
    variant = Variant("1", 4, "T", "G", normalize_contig_names=False)
    left_mate = make_pysam_read(
        seq="ACCGTG",
        cigar="6M",
        name="fragment-1",
        reference_start=0,
    )
    right_mate = make_pysam_read(
        seq="CGTGAA",
        cigar="6M",
        name="fragment-1",
        reference_start=2,
    )
    read_evidence = ReadCollector(
        merge_overlapping_fragments=True
    ).read_evidence_for_variant(
        variant=variant,
        alignment_file=MockAlignmentFile(
            references=("1",),
            reads=[left_mate, right_mate],
        ),
    )
    eq_(len(read_evidence.alt_reads), 1)
    eq_(read_evidence.alt_reads[0].sequence, "ACCGTGAA")

    df = allele_counts_dataframe([(variant, read_evidence)])
    row = df.iloc[0]
    eq_(row.num_alt_reads, 2)
    eq_(row.num_alt_fragments, 1)


if __name__ == "__main__":
    test_allele_count_dataframe()
