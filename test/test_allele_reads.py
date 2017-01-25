
from isovar.allele_reads import AlleleRead
from isovar.locus_reads import LocusRead
from nose.tools import eq_

def make_read_at_locus(prefix, alt, suffix, base_quality=30, name="dummy"):
    dummy_sequence = prefix + alt + suffix
    return LocusRead(
        name="dummy",
        sequence=dummy_sequence,
        reference_positions=list(range(1, len(dummy_sequence) + 1)),
        quality_scores=[base_quality] * len(dummy_sequence))

def test_allele_read_from_single_read_at_locus_trim_N_nucleotides():
    prefix = "NCCN"
    alt = "A"
    suffix = "TNNA"
    read_at_locus = make_read_at_locus(prefix, alt, suffix)
    base0_locus_start = len(prefix)
    base0_locus_end = len(prefix) + len(alt)
    allele_read = AlleleRead.from_locus_read(
        base0_locus_start,
        base0_locus_end,
        read_at_locus)
    print(allele_read)
    expected = AlleleRead(prefix="", allele="A", suffix="T", name="dummy")
    eq_(allele_read, expected)

if __name__ == "__main__":
    test_allele_read_from_single_read_at_locus_trim_N_nucleotides()
