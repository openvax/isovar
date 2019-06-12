from isovar.dataframe_helpers import allele_counts_dataframe
from isovar.allele_read import AlleleRead
from isovar.read_evidence import ReadEvidence
from varcode import Variant
from nose.tools import eq_


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
