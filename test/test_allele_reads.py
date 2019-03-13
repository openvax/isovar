
from isovar.allele_read import AlleleRead
from isovar.locus_read import LocusRead
from nose.tools import eq_


def make_read_at_locus(prefix, alt, suffix, base_quality=30, name="dummy"):
    dummy_sequence = prefix + alt + suffix
    return LocusRead(
        name="dummy",
        sequence=dummy_sequence,
        reference_positions=list(range(1, len(dummy_sequence) + 1)),
        quality_scores=[base_quality] * len(dummy_sequence),
        read_base0_start_inclusive=len(prefix),
        read_base0_end_exclusive=len(prefix) + len(alt),
        reference_base0_start_inclusive=len(prefix),
        reference_base0_end_exclusive=len(prefix) + len(alt))


def test_allele_read_from_single_read_at_locus_trim_N_nucleotides():
    read_at_locus = make_read_at_locus(prefix="NCCN", alt="A", suffix="TNNA")
    allele_read = AlleleRead.from_locus_read(read_at_locus)
    print(allele_read)
    expected = AlleleRead(prefix="", allele="A", suffix="T", name="dummy")
    eq_(allele_read, expected)

if __name__ == "__main__":
    test_allele_read_from_single_read_at_locus_trim_N_nucleotides()
