from isovar.common import reverse_complement_dna
from nose.tools import eq_

def test_reverse_complement_dna():
    eq_(
        "ATGCAATTGGCC",
        reverse_complement_dna("GGCCAATTGCAT"))
