from isovar.alignment_score import alignment_score
from nose.tools import eq_

def test_alignment_score_same_strings():
    eq_(alignment_score("", ""), 0)
    eq_(alignment_score("a", "a"), 0)
    eq_(alignment_score("WUZZLE", "WUZZLE"), 0)

def test_alignment_score_strings_of_different_length():
    eq_(alignment_score("ab", "a"), 1)
    eq_(alignment_score("a", "ab"), 1)
    eq_(alignment_score("WUZZLE", "WUZZLE?"), 1)
    eq_(alignment_score("SNUZZLE", "UZZLE"), 2)


def test_alignment_score_totally_different_strings():
    eq_(alignment_score("", "a"), 1)
    eq_(alignment_score("", "ab"), 2)
    eq_(alignment_score("WUZZLE", "HEAVY"), 6)
    eq_(alignment_score("DOG", "CATCATCAT"), 9)
    eq_(alignment_score("a", ""), 1)
    eq_(alignment_score("ab", ""), 2)
    eq_(alignment_score("WUZZLE", "HEAVY"), 6)
    eq_(alignment_score("CATCATCAT", "DOG"), 9)


def test_alignment_min_subsequence_length():
    # if matching subsequence isn't long enough then misalignment score
    # will be sum of two sequence lengths
    eq_(alignment_score("aaa", "aaa", min_subsequence_length=10), 6)
