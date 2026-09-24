"""Platform-neutral collection, native tags, and #330/#295 regressions."""

from array import array
import json

import pysam
import pytest

from isovar.read_collector import ReadCollector
from isovar.read_metadata import record_evidence
from tests.mock_objects import MockAlignmentFile
from tests.test_sv_rna import Scenario, write_bam


def record(platform=None):
    group = dict(ID="run")
    if platform:
        group["PL"] = platform
    header = pysam.AlignmentHeader.from_dict(dict(SQ=[dict(SN="1", LN=1000)], RG=[group]))
    read = pysam.AlignedSegment(header)
    read.query_name = "read"
    read.query_sequence = "AACCGGTT"
    read.cigarstring = "8M"
    read.reference_start = 100
    read.reference_id = 0
    read.mapping_quality = 60
    read.set_tag("RG", "run")
    return read


@pytest.mark.parametrize("platform", ["ILLUMINA", "ONT", "PACBIO", "DNBSEQ", "IONTORRENT", None])
@pytest.mark.parametrize("compact", [False, True])
def test_platform_independent_quality_and_qc_policy(platform, compact):
    read = record(platform)
    bam = MockAlignmentFile(["1"], [read])

    def collect(collector):
        return collector.collect_locus_reads_with_optional_compaction(bam, "1", 103, 104, compact=compact)

    assert collect(ReadCollector())[0].quality_scores == [None] * 8
    assert collect(ReadCollector(use_reads_without_base_qualities=False)) == []
    # High read-level confidence must never fabricate missing QUAL.
    read.set_tag("rq", 0.99999)
    read.set_tag("qs", 40.0)
    read.set_tag("OQ", "I" * 8)
    assert collect(ReadCollector(use_reads_without_base_qualities=False)) == []
    read.query_qualities = [0, 10, 20, 30, 40, 30, 20, 10]
    assert list(collect(ReadCollector())[0].quality_scores) == list(read.query_qualities)
    read.flag |= 512
    assert collect(ReadCollector()) == []
    assert ReadCollector().alignment_filter_reason(read) == "qc_fail"


def test_native_evidence_keeps_ont_pacbio_and_sam_fields_separate():
    read = record("ONT")
    read.mapping_quality = 255
    read.set_tag("qs", 21.5)
    read.set_tag("dx", -1)  # A simplex parent of a duplex, not a duplex read.
    read.set_tag("pi", "parent-read")
    read.set_tag("sp", 0)
    read.set_tag("NH", 1)
    read.set_tag("AS", 7)
    read.set_tag("RX", "ACGT")
    read.set_tag("QX", "!!!!")
    read.set_tag("MI", "molecule")
    read.set_tag("OQ", "I" * 8)
    read.set_tag("rq", -1.0)  # CCS unpolished sentinel, not an accuracy.
    read.set_tag("ff", 0)
    read.set_tag("mv", array("b", [5] * 10000))
    evidence = record_evidence(read)
    assert evidence["mapping_quality"] is None and evidence["mapping_quality_raw"] == 255
    assert not evidence["base_qualities_available"] and evidence["original_base_qualities_tag_present"]
    tags = evidence["tags"]
    assert {key: tags[key] for key in ("qs", "dx", "pi", "sp", "NH", "AS", "RX", "QX", "MI", "rq", "ff")} == dict(
        qs=21.5, dx=-1, pi="parent-read", sp=0, NH=1, AS=7, RX="ACGT", QX="!!!!", MI="molecule", rq=-1.0, ff=0)
    assert tags["np"] is None and "mv" not in tags and "OQ" not in tags
    json.dumps(evidence, allow_nan=False)
    assert read.query_qualities is None  # No mutation or implicit OQ recovery.


def test_tag_extraction_never_decodes_unrequested_arrays():
    class NoBulkTags:
        def __getattr__(self, name):
            if name in ("get_tags", "tags"):
                raise AssertionError("Do not materialize every auxiliary tag")
            return getattr(read, name)

    read = record("PACBIO")
    read.set_tag("ec", 12.5)
    read.set_tag("ip", array("B", [1] * 10000))
    assert record_evidence(NoBulkTags())["tags"]["ec"] == 12.5


def test_star_mapq_255_is_retained_but_not_claimed_as_phred_255():
    read = record("ILLUMINA")
    read.mapping_quality = 255
    read.set_tag("NH", 1)
    collector = ReadCollector(min_mapping_quality=30)
    assert collector.locus_read_from_pysam_aligned_segment(read, 103, 104) is not None
    assert record_evidence(read)["mapping_quality"] is None


def test_custom_filter_applies_before_small_variant_and_sv_reconstruction(tmp_path):
    s = Scenario()
    reads = s.reads()
    for r in reads:
        r.set_tag("qs", 5.0)
    visited = []

    def keep(read):
        visited.append(read.query_name)
        return not read.has_tag("qs") or read.get_tag("qs") >= 10

    collector = ReadCollector(read_filter=keep)
    for r in reads:
        assert collector.locus_read_from_pysam_aligned_segment(
            r, r.reference_start + 1, r.reference_start + 2) is None
    assert len(visited) == len(reads)
    result = s.run(write_bam(tmp_path / "filtered.bam", reads), read_collector=collector)
    assert result["paths"] == []
    assert result["excluded_records"]["read_filter"] > 0
    assert result["parameters"]["read_collection"]["custom_read_filter"]
    assert result["alignment_metadata"]["read_groups"] == [dict(ID="a"), dict(ID="b")]


def test_custom_filter_is_not_called_on_ineligible_records_and_errors_propagate():
    def fail(read):
        raise RuntimeError("invalid filter configuration")

    read = record()
    collector = ReadCollector(read_filter=fail)
    read.flag = 512
    assert collector.alignment_filter_reason(read) == "qc_fail"
    read.flag = 0
    with pytest.raises(RuntimeError, match="invalid filter configuration"):
        collector.locus_read_from_pysam_aligned_segment(read, 103, 104)
    with pytest.raises(TypeError, match="read_filter"):
        ReadCollector(read_filter=17)


def test_allele_only_collection_never_loads_annotation_for_logging():
    class VariantWithoutAnnotation:
        contig, start, ref, alt = "1", 104, "C", "T"

        @property
        def gene_names(self):
            raise AssertionError("Allele collection must not load a gene database")

        def __str__(self):
            raise AssertionError("Formatting a Variant may access annotation")

    bam = MockAlignmentFile(["1"], [record()])
    evidence = ReadCollector().read_evidence_for_variant(VariantWithoutAnnotation(), bam)
    assert len(evidence.ref_reads) == 1
