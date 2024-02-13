from isovar.allele_read import AlleleRead
from isovar.locus_read import LocusRead
from isovar.read_collector import ReadCollector
from varcode import Variant

from mock_objects import MockAlignmentFile, make_pysam_read
from testing_helpers import assert_equal_fields

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
        reference_base0_end_exclusive=len(prefix) + len(alt),
    )


def test_allele_read_from_single_read_at_locus_trim_N_nucleotides():
    read_at_locus = make_read_at_locus(prefix="NCCN", alt="A", suffix="TNNA")
    allele_read = AlleleRead.from_locus_read(read_at_locus)
    print(allele_read)
    expected = AlleleRead(prefix="", allele="A", suffix="T", name="dummy")
    eq_(allele_read, expected)


def test_allele_read_insertion():
    """
    Test that read gets assigned to the correct allele when there's an insertion on a reference "chr1" which is assumed to be "ACCTTG"
    and the variant sequence is "ACCTGTG"
    """
    pysam_read_with_insertion = make_pysam_read(
        seq="ACCTGTG", cigar="4M1I2M", mdtag="6"
    )
    pysam_read_without_insertion = make_pysam_read(seq="ACCTTG", cigar="6M", mdtag="6")
    variant = Variant("1", 4, ref="T", alt="TG")
    samfile_with_insertion = MockAlignmentFile(
        references={"chromosome"}, reads=[pysam_read_with_insertion]
    )
    read_creator = ReadCollector()
    reads_with_insertion = read_creator.get_locus_reads(
        samfile_with_insertion, "chromosome", variant.start, variant.start
    )

    assert (
        len(reads_with_insertion) == 1
    ), "Expected to get back one read but instead got %d" % (len(reads_with_insertion),)
    locus_read_with_insertion = reads_with_insertion[0]
    allele_read_with_insertion = AlleleRead.from_locus_read(locus_read_with_insertion)

    samfile_without_insertion = MockAlignmentFile(
        references={"chromosome"}, reads=[pysam_read_without_insertion]
    )
    reads_without_insertion = read_creator.get_locus_reads(
        samfile_without_insertion, "chromosome", variant.start, variant.start
    )

    assert (
        len(reads_without_insertion) == 1
    ), "Expected to get back one read but instead got %d" % (
        len(reads_without_insertion),
    )
    locus_read_without_insertion = reads_without_insertion[0]
    allele_read_without_insertion = AlleleRead.from_locus_read(
        locus_read_without_insertion
    )
    expected_allele_read_with_insertion = AlleleRead(
        prefix="ACCT",
        allele="G",
        suffix="TG",
        name="dummy",
    )
    expected_allele_read_without_insertion = AlleleRead(
        prefix="ACCT",
        allele="",
        suffix="TG",
        name="dummy",
    )

    assert_equal_fields(allele_read_with_insertion, expected_allele_read_with_insertion)
    assert_equal_fields(
        allele_read_without_insertion, expected_allele_read_without_insertion
    )


if __name__ == "__main__":
    test_allele_read_from_single_read_at_locus_trim_N_nucleotides()
    test_allele_read_insertion()
