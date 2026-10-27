"""Contig naming is not assembly liftover; ambiguous aliases must fail."""

from itertools import product

import pytest
import pysam
from pyensembl import EnsemblRelease
from varcode import Variant

from isovar.read_collector import ReadCollector


MITO_NAMES = ("M", "MT", "chrM", "chrMT", "m", "mt", "chrm", "chrmt", "CHRM", "CHRMT")


@pytest.mark.parametrize("variant,bam", list(product(MITO_NAMES, repeat=2)))
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


def test_original_accession_alias_survives_annotation_normalization(tmp_path):
    path = tmp_path / 'accession.bam'
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
              'SQ': [{'SN': 'chrUn_KI270442v1', 'LN': 1000}]}
    with pysam.AlignmentFile(path, 'wb', header=header) as bam:
        read = pysam.AlignedSegment(bam.header)
        read.query_name = 'original-witness'
        read.query_sequence = 'A' * 9 + 'G' + 'A' * 20
        read.query_qualities = pysam.qualitystring_to_array('I' * 30)
        read.reference_id = 0
        read.reference_start = 90
        read.mapping_quality = 60
        read.cigarstring = '30M'
        bam.write(read)
    pysam.index(str(path))
    variant = Variant('chrUn_KI270442v1', 100, 'A', 'G', EnsemblRelease(87),
                      convert_ucsc_contig_names=True)
    assert variant.contig == 'KI270442.1'
    with pysam.AlignmentFile(path) as bam:
        evidence = ReadCollector().read_evidence_for_variant(variant, bam)
    assert len(evidence.alt_reads) == 1
    assert evidence.alt_reads[0].name == 'original-witness'
