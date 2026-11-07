"""Contig naming is not assembly liftover; ambiguous aliases must fail."""

from itertools import product
from types import SimpleNamespace

import pytest
import pysam
from pyensembl import EnsemblRelease
from varcode import Variant

from isovar.read_collector import ReadCollector
from isovar import UnresolvedContigError


def resolve(name, names, *, reference_name=None):
    return ReadCollector().resolve_contig(SimpleNamespace(references=names), name,
                                         reference_name=reference_name)


MITO_NAMES = ("M", "MT", "chrM", "chrMT", "m", "mt", "chrm", "chrmt", "CHRM", "CHRMT")


@pytest.mark.parametrize("variant,bam", list(product(MITO_NAMES, repeat=2)))
def test_unique_mitochondrial_alias(variant, bam):
    assert resolve(variant, {"1", bam}).resolved_name == bam


@pytest.mark.parametrize("name", MITO_NAMES + ("1", "chr1", "X", "chrX"))
def test_exact_contig_preferred_to_aliases(name):
    assert resolve(name, set(MITO_NAMES) | {"1", "chr1", "X", "chrX"}).resolved_name == name


@pytest.mark.parametrize("name,names", [("MT", {"M", "chrM"}), ("MT", {"chrM", "chrMT"}),
                                       ("x", {"X", "chrX"}), ("chr1", {"1", "CHR1"})])
def test_ambiguous_alias_is_not_arbitrary_or_zero_reads(name, names):
    result = resolve(name, names)
    assert result.status == "ambiguous"
    assert result.resolved_name is None
    assert set(result.candidates) == names


@pytest.mark.parametrize("name,names,expected", [("1", {"chr1"}, "chr1"), ("chr1", {"1"}, "1"),
                                                ("X", {"chrx"}, "chrx"), ("MT", {"M_random", "chr1"}, None),
                                                ("2", {"chr1"}, None), ("GL0001.1", {"GL0001.1"}, "GL0001.1")])
def test_nuclear_and_missing_contigs(name, names, expected):
    assert resolve(name, names).resolved_name == expected


@pytest.mark.parametrize("names,status", [(None, "reference_unavailable"), ([], "missing"),
                                         (["M", "chrM"], "ambiguous"), (["1"], "missing")])
def test_unresolved_collection_is_not_zero_coverage(names, status):
    variant = Variant("MT", 10, "A", "G")
    alignment = SimpleNamespace(references=names, filename=b"rna.bam")
    with pytest.raises(UnresolvedContigError, match="unresolved_chromosome") as caught:
        ReadCollector().read_evidence_for_variant(variant, alignment)
    assert caught.value.resolution.status == status
    assert caught.value.to_dict()["source"] == "rna.bam"


def test_build_context_and_noncanonical_names_are_preserved():
    assert resolve("GL000216.1", ["chrUn_gl000216"], reference_name="GRCh37").status == "alias"
    assert resolve("GL000216.1", ["chrUn_gl000216"], reference_name="GRCh38").status == "missing"
    names = ["0", "HPV16", "NC_007605.1", "hs37d5", "HLA-A*01:01:01:01", "chr6_apd_hap1"]
    for name in names:
        assert resolve(name, names).status == "exact"
    for name in ["EBV", "HPV-16", "NC_007605.2", "6"]:
        assert resolve(name, names).status == "missing"


def test_cache_is_bounded_and_scoped_to_header_and_assembly(monkeypatch):
    import isovar.contig_resolution as module
    original = module.ContigResolver
    constructions = []

    def construct(names, **kwargs):
        constructions.append((names, kwargs))
        return original(names, **kwargs)

    monkeypatch.setattr(module, "ContigResolver", construct)
    collector = ReadCollector()
    first = SimpleNamespace(references=("chrUn_gl000216",))
    second = SimpleNamespace(references=("GL000216.1",))
    for _ in range(20):
        assert collector.resolve_contig(first, "GL000216.1", reference_name="GRCh37").status == "alias"
    assert len(constructions) == 1
    assert collector.resolve_contig(first, "GL000216.1", reference_name="GRCh38").status == "missing"
    assert collector.resolve_contig(second, "GL000216.1", reference_name="GRCh37").status == "exact"
    assert len(constructions) == 3
    for i in range(12):
        collector.resolve_contig(SimpleNamespace(references=(str(i + 1),)), str(i + 1))
    assert len(collector._contig_cache._indexes) == 8
    assert collector.resolve_contig(first, "GL000216.1", reference_name="GRCh37").status == "alias"
    assert len(constructions) == 16


def test_actual_header_cache_does_not_reuse_closed_or_reopened_handle(tmp_path):
    path = tmp_path / "header.bam"
    with pysam.AlignmentFile(path, "wb", header={"SQ": [{"SN": "chr1", "LN": 1000}]}):
        pass
    collector = ReadCollector()
    with pysam.AlignmentFile(path) as bam:
        assert bam.header is bam.header
        assert collector.resolve_contig(bam, "1").resolved_name == "chr1"
        assert collector.resolve_contig(bam, "chr1").status == "exact"
    with pytest.raises(ValueError, match="closed alignment"):
        collector.resolve_contig(bam, "1")
    with pysam.AlignmentFile(path) as reopened:
        assert collector.resolve_contig(reopened, "1").resolved_name == "chr1"
    assert len(collector._contig_cache._indexes) == 2


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


def test_original_exact_name_wins_over_normalized_exact_name(tmp_path):
    path = tmp_path / "two-spellings.bam"
    with pysam.AlignmentFile(path, "wb", header={"SQ": [
            {"SN": "1", "LN": 1000}, {"SN": "chr1", "LN": 1000}]}) as bam:
        read = pysam.AlignedSegment(bam.header)
        read.query_name = "source-exact"
        read.query_sequence = "A" * 9 + "G" + "A" * 20
        read.query_qualities = pysam.qualitystring_to_array("I" * 30)
        read.reference_id, read.reference_start = 1, 90
        read.mapping_quality, read.cigarstring = 60, "30M"
        bam.write(read)
    pysam.index(str(path))
    variant = Variant("chr1", 100, "A", "G", EnsemblRelease(87), convert_ucsc_contig_names=True)
    assert variant.contig == "1"
    with pysam.AlignmentFile(path) as bam:
        evidence = ReadCollector().read_evidence_for_variant(variant, bam)
    assert [read.name for read in evidence.alt_reads] == ["source-exact"]


@pytest.mark.parametrize("assembly", ["GRCh37", "hg19", "b37", "NCBI37", " grch37 "])
@pytest.mark.parametrize("query,target", [
    ("MT", "chrM"), ("chrM", "MT"), ("M", "chrMT"), ("chrMT", "M"),
    ("mt", "CHRM"), ("chrm", "CHRMT"),
])
def test_grch37_distinct_mitochondrial_sequences_are_not_aliases(assembly, query, target):
    result = resolve(query, [target], reference_name=assembly)
    assert result.status == "missing"
    assert result.resolved_name is None and result.candidates == ()
    assert result.original_name == query and result.reference_name == assembly
    assert result.normalized_name == ("MT" if query.upper() in ("MT", "CHRMT") else "chrM")


@pytest.mark.parametrize("query,names,status,selected,candidates", [
    ("M", ["chrM", "MT", "chrMT"], "alias", "chrM", ("chrM",)),
    ("mt", ["M", "chrM", "chrMT"], "alias", "chrMT", ("chrMT",)),
    ("M", ["chrM", "CHRM", "MT"], "ambiguous", None, ("CHRM", "chrM")),
    ("MT", ["chrMT", "mt", "chrM"], "ambiguous", None, ("chrMT", "mt")),
    ("chrM", ["chrM", "M", "MT"], "exact", "chrM", ("chrM",)),
    ("MT", ["MT", "M", "chrMT"], "exact", "MT", ("MT",)),
    ("MT", None, "reference_unavailable", None, ()),
])
def test_grch37_resolution_retains_only_same_sequence_candidates(query, names, status, selected, candidates):
    result = resolve(query, names, reference_name="GRCh37")
    assert (result.status, result.resolved_name, result.candidates) == (status, selected, candidates)


def test_cached_mitochondrial_lookup_remains_assembly_scoped():
    collector = ReadCollector()
    bam = SimpleNamespace(references=("chrM",))
    for assembly, status in [("GRCh38", "alias"), ("GRCh37", "missing"), ("GRCh38", "alias")]:
        assert collector.resolve_contig(bam, "MT", reference_name=assembly).status == status


def mitochondrial_bam(path, names):
    header = {"SQ": [{"SN": name, "LN": 16571 if name == "chrM" else 16569} for name in names]}
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for index, name in enumerate(names):
            read = pysam.AlignedSegment(bam.header)
            read.query_name = "witness-" + name
            read.query_sequence = "A" * 9 + "G" + "A" * 20
            read.query_qualities = pysam.qualitystring_to_array("I" * 30)
            read.reference_id, read.reference_start = index, 90
            read.mapping_quality, read.cigarstring = 60, "30M"
            bam.write(read)
    pysam.index(str(path))
    return path


@pytest.mark.parametrize("query,target", [("MT", "chrM"), ("chrM", "MT")])
def test_grch37_collection_does_not_read_another_mitochondrial_sequence(tmp_path, query, target):
    path = mitochondrial_bam(tmp_path / "mitochondria.bam", [target])
    variant = Variant(query, 100, "A", "G", EnsemblRelease(75), normalize_contig_names=False)
    with pysam.AlignmentFile(path) as bam:
        with pytest.raises(UnresolvedContigError) as caught:
            ReadCollector().read_evidence_for_variant(variant, bam)
    assert caught.value.resolution.status == "missing"
    assert caught.value.resolution.original_name == query


@pytest.mark.parametrize("query,expected", [("M", "chrM"), ("MT", "chrMT"), ("chrM", "chrM")])
def test_grch37_collection_selects_the_same_reference_from_mixed_header(tmp_path, query, expected):
    path = mitochondrial_bam(tmp_path / "mixed.bam", ["chrM", "chrMT"])
    variant = Variant(query, 100, "A", "G", EnsemblRelease(75), convert_ucsc_contig_names=True)
    assert variant.contig == ("chrM" if expected == "chrM" else "MT")
    with pysam.AlignmentFile(path) as bam:
        evidence = ReadCollector().read_evidence_for_variant(variant, bam)
    assert [read.name for read in evidence.alt_reads] == ["witness-" + expected]


def test_grch37_original_library_query_keeps_mitochondrial_references_separate(tmp_path):
    from isovar.rna_library import prepare_rna_library
    path = mitochondrial_bam(tmp_path / "library.bam", ["chrM"])
    with prepare_rna_library(path, tmp_path / "index", reference_name="GRCh37").open() as library:
        with pytest.raises(ValueError, match="absent"):
            list(library.templates([("MT", 90, 120)]))
        matches = list(library.templates([("M", 90, 120)]))
        assert [[read.query_name for read in template] for template in matches] == [["witness-chrM"]]
