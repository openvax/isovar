from isovar.dna import reverse_complement_dna, complement_dna
from nose.tools import eq_


def test_reverse_complement_dna():
    eq_("ATGCAATTGGCC", reverse_complement_dna("GGCCAATTGCAT"))


def test_complement_dna():
    eq_("ATGC", complement_dna("TACG"))
