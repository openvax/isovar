from isovar.allele_counts import allele_counts_dataframe
from isovar.allele_reads import AlleleRead
from varcode import Variant
from nose.tools import eq_

def test_allele_count_dataframe():
    variant = Variant("test_contig", 50, "C", "G")
    reads = [
        AlleleRead(prefix="AAA", allele="C", suffix="TTT", name="C1"),
        AlleleRead(prefix="AAC", allele="C", suffix="TTA", name="C2"),
        AlleleRead(prefix="AAA", allele="G", suffix="TTT", name="G1"),
    ]
    df = allele_counts_dataframe([(variant, reads)])
    assert len(df) == 1, "Wrong number of rows in DataFrame: %s" % (df,)
    row = df.irow(0)
    eq_(row.n_ref, 2)
    eq_(row.n_alt, 1)
    eq_(row.n_other, 0)

if __name__ == "__main__":
    test_allele_count_dataframe()
