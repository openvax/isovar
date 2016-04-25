import isovar
import isovar.variant_read
from isovar.variant_read import (
    variant_read_from_single_read_at_locus,
    VariantRead,
)
from isovar.read_at_locus import ReadAtLocus
from nose.tools import eq_

def make_read_at_locus(prefix, alt, suffix, base_quality=30, name="dummy"):
    dummy_sequence = prefix + alt + suffix
    return ReadAtLocus(
        name="dummy",
        sequence=dummy_sequence,
        reference_positions=list(range(1, len(dummy_sequence) + 1)),
        quality_scores=[base_quality] * len(dummy_sequence),
        base0_read_position_before_variant=len(prefix) - 1,
        base0_read_position_after_variant=len(prefix) + len(alt),
    )

def test_variant_read_from_single_read_at_locus_trim_N_nucleotides():
    read_at_locus =make_read_at_locus(prefix="NCCN", alt="A", suffix="TNNA")
    variant_read = variant_read_from_single_read_at_locus(
        read_at_locus, ref="T", alt="A")
    print(variant_read)
    expected = VariantRead(prefix="", alt="A", suffix="T", name="dummy")
    eq_(variant_read, expected)

if __name__ == "__main__":
    test_variant_read_from_single_read_at_locus_trim_N_nucleotides()
