"""Original CC0 tumor RNA reads at osteosarc.com somatic mutation loci.

The independent ledger describes CIGAR evidence, not somatic truth or
immunogenicity. No network, transcript annotation or reference downloads.
"""

from collections import Counter
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import traceback
from types import SimpleNamespace

import pysam
import pytest

from isovar.allele_read import AlleleRead
from isovar.read_collector import ReadCollector
from isovar.variant_helpers import base0_interval_for_variant
from isovar.variant_sequence_creator import VariantSequenceCreator

from .real_rna_helpers import cigar_observation, record_digest


DATA = Path(__file__).parent / "data" / "osteosarc"
MANIFEST = json.loads((DATA / "manifest.json").read_text())
DATASETS = MANIFEST["datasets"]
CASES = [pytest.param(name, v, id=f"{name}-{v['gene']}")
         for name, data in DATASETS.items() for v in data["variants"]]
SNPS = [pytest.param(name, v, id=f"{name}-{v['gene']}")
        for name, data in DATASETS.items() for v in data["variants"]
        if len(v["ref"]) == len(v["alt"]) == 1]


class MissingSpliceReferenceBase(AssertionError):
    """Only the documented normalization crash may count as an expected failure."""


SPLICE_NORMALIZATION_BUG = pytest.mark.xfail(
    strict=True, raises=MissingSpliceReferenceBase,
    reason="https://github.com/openvax/isovar/issues/217: CIGAR N has no reference base to normalize")


def variant_object(record):
    return SimpleNamespace(contig=record["chrom"].removeprefix("chr"),
                           start=record["pos"], ref=record["ref"], alt=record["alt"],
                           gene_names=(), short_description=record["variant_id"])


def primary_observations(record):
    return [o for o in record["observations"] if not o["flag"] & (256 | 1024 | 2048)]


@pytest.fixture(scope="module")
def alignments(tmp_path_factory):
    directory = tmp_path_factory.mktemp("osteosarc")
    result = {}
    for name, data in DATASETS.items():
        sam = directory / (name + ".sam")
        sam.write_bytes(gzip.decompress((DATA / data["file"]).read_bytes()))
        bam = directory / (name + ".bam")
        pysam.sort("--no-PG", "-o", str(bam), str(sam))
        pysam.index(str(bam))
        result[name] = bam
    return result


@pytest.mark.parametrize("name", DATASETS)
def test_osteosarc_original_records_and_selection_are_preserved(name, alignments):
    data = DATASETS[name]
    raw = gzip.decompress((DATA / data["file"]).read_bytes())
    assert sha256(raw).hexdigest() == data["sam_sha256"]
    with pysam.AlignmentFile(alignments[name]) as bam:
        reads = list(bam)
    assert len(reads) == data["fixture_records"]
    # Do not turn original duplicate SAM records into unique observations.
    assert Counter(map(record_digest, reads)) == Counter(data["records"])
    assert {r.query_name for r in reads} == (
        set(data["baseline_template_names"]) | set(data["quality_stress_template_names"]))
    for variant in data["variants"]:
        observed = [o for read in reads if (o := cigar_observation(read, variant)) is not None]
        assert sorted(observed, key=lambda o: o["record"]) == sorted(
            variant["observations"], key=lambda o: o["record"])


@pytest.mark.parametrize("name,record", CASES)
@pytest.mark.parametrize("secondary", [False, True])
def test_osteosarc_direct_cigar_alleles_and_original_quality_indices(name, record, secondary, alignments):
    start, end = base0_interval_for_variant(variant_object(record))
    expected = {o["record"]: o for o in record["observations"]}
    collector = ReadCollector(min_mapping_quality=0, use_secondary_alignments=secondary,
                              use_soft_clipped_bases=True)
    tested = 0
    with pysam.AlignmentFile(alignments[name]) as bam:
        for read in bam:
            observation = expected.get(record_digest(read))
            if observation is None or read.is_duplicate or (read.is_secondary and not secondary):
                continue
            # This tests direct CIGAR extraction, not repeat normalization.
            locus = collector.locus_read_from_pysam_aligned_segment(read, start, end)
            assert locus is not None, read.query_name
            allele = AlleleRead.from_locus_read(locus)
            assert allele is not None
            assert allele.allele == observation["allele"], read.query_name
            assert [locus.read_base0_start_inclusive, locus.read_base0_end_exclusive] == observation["query_interval"]
            assert [locus.quality_scores[i] for i in observation["quality_indices"]] == observation["qualities"]
            tested += 1
    assert tested > 0


@pytest.mark.parametrize("name,record", SNPS)
@pytest.mark.parametrize("min_mapq", [0, 1, 4, 20, 60])
def test_osteosarc_snp_partition_matches_independent_evidence(name, record, min_mapq, alignments):
    expected = Counter((o["name"], o["allele"]) for o in record["observations"]
                       if not o["flag"] & (256 | 1024) and o["mapq"] >= min_mapq)
    with pysam.AlignmentFile(alignments[name]) as bam:
        evidence = ReadCollector(min_mapping_quality=min_mapq, use_secondary_alignments=False).read_evidence_for_variant(
            variant_object(record), bam)
    actual = Counter((r.name, r.allele) for r in evidence.ref_reads + evidence.alt_reads + evidence.other_reads)
    assert actual == expected


@pytest.mark.parametrize("name,gene", [
    ("bulk_star_t0", "H1-2"), ("bulk_star_t0", "GTF3C5"),
    pytest.param("ont_t1", "H1-2", marks=SPLICE_NORMALIZATION_BUG),
    pytest.param("ont_t1", "GTF3C5", marks=SPLICE_NORMALIZATION_BUG),
    pytest.param("ont_t1", "PIP5K1A", marks=SPLICE_NORMALIZATION_BUG),
])
def test_osteosarc_cancer_deletion_support_reaches_sequence_creation(name, gene, alignments):
    record = next(v for v in DATASETS[name]["variants"] if v["gene"] == gene)
    direct_alt_names = {o["name"] for o in primary_observations(record) if o["allele"] == ""}
    assert len(direct_alt_names) >= 2
    variant = variant_object(record)
    with pysam.AlignmentFile(alignments[name]) as bam:
        try:
            evidence = ReadCollector(use_secondary_alignments=False).read_evidence_for_variant(variant, bam)
        except AttributeError as error:
            frames = [frame.f_code.co_name for frame, _ in traceback.walk_tb(error.__traceback__)]
            if (str(error) == "'NoneType' object has no attribute 'upper'"
                    and "_left_aligned_indel_interval_for_variant" in frames):
                raise MissingSpliceReferenceBase(gene) from error
            raise
    assert direct_alt_names <= {r.name for r in evidence.alt_reads}
    sequences = VariantSequenceCreator(min_variant_sequence_coverage=1,
                                       variant_sequence_assembly=False).reads_to_variant_sequences(variant, evidence.alt_reads)
    assert sequences
    assert all(s.alt == "" and len(s) > 0 for s in sequences)
    assert all(s.read_names <= {r.name for r in evidence.alt_reads} for s in sequences)


@pytest.mark.parametrize("name", DATASETS)
def test_map2_dna_report_does_not_create_rna_deletion_support(name, alignments):
    record = next(v for v in DATASETS[name]["variants"] if v["gene"] == "MAP2")
    with (DATA / "source_variant_counts.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    dna = [r for r in rows if r["variant_id"] == record["variant_id"] and r["bam_file"] == "tumor_bg22.recal.bam"]
    assert len(dna) == 1
    assert int(dna[0]["alt_reads"]) == 135  # Source report, not our RNA oracle.
    observations = primary_observations(record)
    assert observations and all(o["allele"] != "" for o in observations)
    with pysam.AlignmentFile(alignments[name]) as bam:
        evidence = ReadCollector(use_secondary_alignments=False).read_evidence_for_variant(variant_object(record), bam)
    assert not evidence.alt_reads
    # This sparse/soft-clipped RNA alignment is NOT proof of absent expression.


def test_osteosarc_star_missing_md_and_mapping_conventions(alignments):
    with pysam.AlignmentFile(alignments["bulk_star_t0"]) as bam:
        reads = list(bam)
        assert any(p.get("PN") == "STAR" and p.get("VN") == "2.7.11b"
                   for p in bam.header.to_dict()["PG"])
    assert all(not r.has_tag("MD") for r in reads)  # Despite source filename .md.bam.
    assert all((r.mapping_quality == 255) == (r.get_tag("NH") == 1) for r in reads)
    assert any(r.is_secondary for r in reads)
    assert any(n > 1 for n in Counter(r.query_name for r in reads).values())


def test_osteosarc_ont_cell_umi_tags_are_preserved_not_reinterpreted(alignments):
    with pysam.AlignmentFile(alignments["ont_t1"]) as bam:
        reads = list(bam)
        assert any(p.get("PN") == "minimap2" and p.get("VN") == "2.24-r1122"
                   for p in bam.header.to_dict()["PG"])
    assert all(r.has_tag("MD") and r.has_tag("CB") and r.has_tag("UB") for r in reads)
    assert any("N" in r.cigarstring for r in reads)
    assert any(r.is_reverse for r in reads)
    assert not any(r.is_secondary for r in reads)
    # This is already tagged/deduplicated single-cell long-read RNA, not an
    # independent direct-RNA replicate of the older NA12878 collection.


@pytest.mark.parametrize("name,record", SNPS)
def test_osteosarc_stress_subset_keeps_low_quality_cancer_alleles(name, record):
    alternate = [o for o in primary_observations(record) if o["allele"] == record["alt"]]
    assert alternate
    assert any(min(o["qualities"]) < 20 for o in alternate)
    assert any(min(o["qualities"]) >= 30 for o in alternate)
    # Raw MAPQ/BQ must remain separate; this does not endorse a cutoff or a
    # calibrated probabilistic combination of the reported scores.
    assert all(o["mapq"] >= 60 for o in alternate)


def test_osteosarc_cigar_oracle_requires_the_same_chromosome(alignments):
    record = next(v for v in DATASETS["ont_t1"]["variants"] if v["gene"] == "EXOC4")
    with pysam.AlignmentFile(alignments["ont_t1"]) as bam:
        read = next(r for r in bam if cigar_observation(r, record) is not None)
    assert cigar_observation(read, dict(record, chrom="chr8")) is None


def test_osteosarc_variant_definitions_match_source_records():
    selection = json.loads((DATA / "selection.json").read_text())
    with (DATA / "source_variant_counts.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for record in selection["variants"]:
        reports = [r for r in rows if r["variant_id"] == record["variant_id"]]
        assert reports
        assert all(all(r[k] == v for k, v in record.items()) for r in reports)
        for data in DATASETS.values():
            matches = [v for v in data["variants"] if v["variant_id"] == record["variant_id"]]
            assert len(matches) == 1
            assert all(str(matches[0][k]) == v for k, v in record.items())
