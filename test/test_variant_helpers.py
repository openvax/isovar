from isovar.variant_helpers import (
    trim_variant_fields,
    trim_variant,
    base0_interval_for_variant
)
from nose.tools import eq_
from varcode import Variant


def test_trim_variant_substitution():
    loc, ref, alt = trim_variant(Variant("chr1", 10, "C", "G"))
    eq_(loc, 10)
    eq_(ref, "C")
    eq_(alt, "G")

def test_trim_variant_fields_substitution():
    loc, ref, alt = trim_variant_fields(10, "C", "G")
    eq_(loc, 10)
    eq_(ref, "C")
    eq_(alt, "G")

def test_trim_variant_insertion():
    loc, ref, alt = trim_variant(Variant("chr1", 10, "C", "CG"))
    eq_(loc, 10)
    eq_(ref, "")
    eq_(alt, "G")


def test_trim_variant_fields_insertion():
    loc, ref, alt = trim_variant_fields(10, "C", "CG")
    eq_(loc, 10)
    eq_(ref, "")
    eq_(alt, "G")


def test_trim_variant_deletion():
    loc, ref, alt = trim_variant(Variant("chr1", 10, "CG", "C"))
    eq_(loc, 11)
    eq_(ref, "G")
    eq_(alt, "")

def test_trim_variant_fields_deletion():
    loc, ref, alt = trim_variant_fields(10, "CG", "C")
    eq_(loc, 11)
    eq_(ref, "G")
    eq_(alt, "")

def test_base0_interval_for_variant_substitution():
    (start, end) = base0_interval_for_variant(Variant("chr1", 10, "C", "G"))
    eq_(start, 9)
    eq_(end, 10)

def test_base0_interval_for_variant_insertion():
    (start, end) = base0_interval_for_variant(Variant("chr1", 10, "C", "CG"))
    eq_(start, 9)
    eq_(end, 9)

def test_base0_interval_for_variant_deletion():
    (start, end) = base0_interval_for_variant(Variant("chr1", 10, "CG", "C"))
    eq_(start, 10)
    eq_(end, 11)
