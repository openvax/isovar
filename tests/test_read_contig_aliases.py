"""Contig naming is not assembly liftover; ambiguous aliases must fail."""

from itertools import product

import pytest

from isovar.read_collector import ReadCollector


MITO_NAMES = ("M", "MT", "chrM", "chrMT", "m", "mt", "chrm", "chrmt", "CHRM", "CHRMT")


@pytest.mark.parametrize("variant,bam", product(MITO_NAMES, repeat=2))
def test_unique_mitochondrial_alias(variant, bam):
    assert ReadCollector._infer_chromosome_name(variant, {"1", bam}) == bam


@pytest.mark.parametrize("name", MITO_NAMES + ("1", "chr1", "X", "chrX"))
def test_exact_contig_preferred_to_aliases(name):
    assert ReadCollector._infer_chromosome_name(name, set(MITO_NAMES) | {"1", "chr1", "X", "chrX"}) == name


@pytest.mark.parametrize("name,names", [("MT", {"M", "chrM"}), ("MT", {"chrM", "chrMT"}),
                                       ("x", {"X", "chrX"}), ("chr1", {"1", "CHR1"})])
def test_ambiguous_alias_is_not_arbitrary_or_zero_reads(name, names):
    with pytest.raises(ValueError, match="Ambiguous alignment contig aliases"):
        ReadCollector._infer_chromosome_name(name, names)


@pytest.mark.parametrize("name,names,expected", [("1", {"chr1"}, "chr1"), ("chr1", {"1"}, "1"),
                                                ("X", {"chrx"}, "chrx"), ("MT", {"M_random", "chr1"}, None),
                                                ("2", {"chr1"}, None), ("GL0001.1", {"GL0001.1"}, "GL0001.1")])
def test_nuclear_and_missing_contigs(name, names, expected):
    assert ReadCollector._infer_chromosome_name(name, names) == expected
