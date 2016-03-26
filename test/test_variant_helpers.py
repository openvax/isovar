from isovar.variant_helpers import trim_variant_fields

def test_trim_variant_fields_substitution():
    loc, ref, alt = trim_variant_fields(10, "C", "G")
    assert loc == 10
    assert ref == "C"
    assert alt == "G"


def test_trim_variant_fields_insertion():
    loc, ref, alt = trim_variant_fields(10, "C", "CG")
    assert loc == 11
    assert ref == ""
    assert alt == "G"

def test_trim_variant_fields_insertion_before_ref():
    loc, ref, alt = trim_variant_fields(10, "C", "TC")
    # ideally we'd want this variant to get normalized to the previous reference
    # base but since this info isn't available we settle for treating
    # this variant as a substitution.
    assert loc == 10
    assert ref == "C"
    assert alt == "TC"

def test_trim_variant_fields_deletion():
    loc, ref, alt = trim_variant_fields(10, "CG", "C")
    assert loc == 11
    assert ref == "G"
    assert alt == ""

def test_trim_variant_fields_deletion_before_ref():
    loc, ref, alt = trim_variant_fields(10, "CG", "G")
    assert loc == 10
    assert ref == "C"
    assert alt == ""
