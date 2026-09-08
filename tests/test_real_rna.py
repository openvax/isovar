"""Offline regressions using original public RNA alignments, not simulated reads."""

from collections import Counter
import gzip
from hashlib import sha256
import json
from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest

from isovar.allele_read import AlleleRead
from isovar.read_collector import ReadCollector
from isovar.variant_helpers import base0_interval_for_variant, trim_variant
from isovar.variant_sequence_creator import VariantSequenceCreator

from .real_rna_helpers import cigar_observation, record_digest


DATA = Path(__file__).parent / "data" / "real_rna"
MANIFEST = json.loads((DATA / "manifest.json").read_text())
DATASETS = MANIFEST["datasets"]
CASES = [(name, v) for name, data in DATASETS.items() for v in data["variants"]]


@pytest.fixture(scope="module")
def alignments(tmp_path_factory):
    """Generate disposable BAM indexes locally; no network/genome downloads."""
    directory = tmp_path_factory.mktemp("real_rna")
    result = {}
    for name, data in DATASETS.items():
        sam = directory / (name + ".sam")
        sam.write_bytes(gzip.decompress((DATA / data["file"]).read_bytes()))
        bam = directory / (name + ".bam")
        pysam.sort("--no-PG", "-o", str(bam), str(sam))
        pysam.index(str(bam))
        result[name] = bam
    return result


def variant_from_record(record):
    # Read extraction uses these fields only. No gene annotation/reference
    # downloads are needed to validate a CIGAR against a known VCF allele.
    return SimpleNamespace(contig=record["chrom"].removeprefix("chr"),
                           start=record["pos"], ref=record["ref"], alt=record["alt"],
                           gene_names=(), short_description=str(record["pos"]))


@pytest.mark.parametrize("name", DATASETS)
def test_public_rna_integrity_and_frozen_cigar_ledger(name, alignments):
    data = DATASETS[name]
    raw = gzip.decompress((DATA / data["file"]).read_bytes())
    assert sha256(raw).hexdigest() == data["sam_sha256"]
    with pysam.AlignmentFile(alignments[name]) as bam:
        reads = list(bam)
    assert len(reads) == data["fixture_records"]
    assert Counter(record_digest(r) for r in reads) == Counter(data["records"])
    for variant in data["variants"]:
        observations = [obs for r in reads if (obs := cigar_observation(r, variant)) is not None]
        assert sorted(observations, key=lambda o: o["record"]) == sorted(
            variant["observations"], key=lambda o: o["record"])


@pytest.mark.parametrize("name,record", CASES, ids=[f"{n}-{v['pos']}" for n, v in CASES])
@pytest.mark.parametrize("secondary", [False, True])
def test_public_rna_direct_observation_preserves_allele_and_quality(name, record, secondary, alignments):
    variant = variant_from_record(record)
    collector = ReadCollector(min_mapping_quality=0, use_secondary_alignments=secondary,
                              use_soft_clipped_bases=True)
    start, end = base0_interval_for_variant(variant)
    trimmed_start, ref, alt = trim_variant(variant)
    expected = {obs["record"]: obs for obs in record["observations"]}
    tested = 0
    with pysam.AlignmentFile(alignments[name]) as bam:
        for read in bam:
            obs = expected.get(record_digest(read))
            if obs is None or read.is_duplicate or (read.is_secondary and not secondary):
                continue
            # Direct extraction and normalization are separate contracts.
            # The STAR insertion has real right-shifted equivalents; the
            # unnormalized CIGAR ledger must not relabel these as reference
            # support in the variant-aware normalization test below.
            locus = collector.locus_read_from_pysam_aligned_segment(read, start, end)
            assert locus is not None, (record["pos"], read.query_name)
            allele = AlleleRead.from_locus_read(locus)
            assert allele is not None
            assert allele.allele == obs["allele"], (record["pos"], read.query_name)
            assert [locus.read_base0_start_inclusive, locus.read_base0_end_exclusive] == obs["query_interval"]
            assert [locus.quality_scores[i] for i in obs["quality_indices"]] == obs["qualities"]
            tested += 1
    assert tested > 0


@pytest.mark.parametrize("name,record", [(n, v) for n, v in CASES if len(v["ref"]) == len(v["alt"]) == 1])
@pytest.mark.parametrize("min_mapq", [0, 1, 4, 20, 60])
def test_public_rna_snp_partition_matches_independent_cigar_counts(name, record, min_mapq, alignments):
    expected = Counter((obs["name"], obs["allele"]) for obs in record["observations"]
                       if not obs["flag"] & (256 | 1024) and obs["mapq"] >= min_mapq)
    with pysam.AlignmentFile(alignments[name]) as bam:
        evidence = ReadCollector(min_mapping_quality=min_mapq, use_secondary_alignments=False).read_evidence_for_variant(
            variant_from_record(record), bam)
    actual = Counter((r.name, r.allele) for r in evidence.ref_reads + evidence.alt_reads + evidence.other_reads)
    assert actual == expected


def test_star_mapq_255_and_nh_are_retained_as_reported(alignments):
    with pysam.AlignmentFile(alignments["illumina_star"]) as bam:
        reads = list(bam)
    assert {r.mapping_quality for r in reads} == {0, 3, 255}
    assert all((r.mapping_quality == 255) == (r.get_tag("NH") == 1) for r in reads)
    assert any(r.is_secondary for r in reads)
    assert any("N" in r.cigarstring for r in reads)
    assert all(r.has_tag("MD") for r in reads)


def test_pacbio_fixture_contains_constant_qualities_and_split_alignments(alignments):
    with pysam.AlignmentFile(alignments["pacbio_minimap2"]) as bam:
        reads = list(bam)
    assert {tuple(sorted(set(r.query_qualities))) for r in reads} == {(20,), (40,)}
    assert all(not r.has_tag("NH") for r in reads)
    assert any(r.is_supplementary and r.has_tag("SA") for r in reads)
    assert any(r.is_reverse for r in reads)
    assert all("N" in r.cigarstring for r in reads)


def test_real_nanopore_snp_demonstrates_q20_support_loss():
    record = next(v for v in DATASETS["ont_minimap2"]["variants"] if v["pos"] == 43785588)
    alt = [o for o in record["observations"] if o["base"] == "A" and not o["flag"] & (256 | 1024 | 2048)]
    # This is a sensitivity observation, not a recommendation or calibration.
    assert len(alt) == 106
    assert sum(o["base_quality"] >= 10 for o in alt) == 63
    assert sum(o["base_quality"] >= 20 for o in alt) == 4
    assert sum(o["base_quality"] >= 30 for o in alt) == 0


def test_real_star_snp_contains_low_quality_alternate_bases():
    record = DATASETS["illumina_star"]["variants"][0]
    alt = [o for o in record["observations"] if o["allele"] == "C" and not o["flag"] & (256 | 1024 | 2048)]
    assert len(alt) == 8
    assert sum(min(o["qualities"]) < 20 for o in alt) == 2
    assert all(o["mapq"] in (3, 255) for o in alt)


def test_all_selected_variants_are_in_independent_dna_benchmark_regions():
    regions = [line.split("\t") for line in (DATA / "HG001.selected.bed").read_text().splitlines()]
    records = {line for line in (DATA / "HG001.selected.vcf").read_text().splitlines() if not line.startswith("#")}
    for _, variant in CASES:
        assert variant["vcf_record"] in records
        assert any(chrom == variant["chrom"] and int(start) <= variant["pos"] - 1
                   and int(end) >= variant["pos"] - 1 + len(variant["ref"])
                   for chrom, start, end in regions)


def test_real_star_repeat_shifted_insertion_keeps_canonical_sequence_and_original_quality(alignments):
    with pysam.AlignmentFile(alignments["illumina_star"]) as bam:
        read = next(r for r in bam if r.query_name == "SRR1258218.15092602" and r.flag == 163)
    assert read.reference_start == 337853
    assert read.cigarstring == "52M1I22M"
    assert read.get_tag("MD") == "70T3"
    assert read.query_sequence[51:53] == "CC"
    assert list(read.query_qualities[51:53]) == [28, 35]
    # The original inserted C is q52 at reference boundary 337905. Since
    # the preceding reference base is also C, the equivalent left shift is
    # q51 at boundary 337904. These TWO quality scores must not be confused.
    locus = ReadCollector(use_soft_clipped_bases=True).locus_read_from_pysam_aligned_segment(
        read, 337904, 337904, trimmed_base1_start=337904, trimmed_ref="", trimmed_alt="C")
    assert locus is not None
    assert locus.read_base0_start_inclusive == 51
    assert locus.read_base0_end_exclusive == 52
    assert list(locus.quality_scores[51:53]) == [28, 35]
    assert AlleleRead.from_locus_read(locus).allele == "C"


@pytest.mark.parametrize("name,record", [(n, v) for n, v in CASES if len(v["ref"]) != len(v["alt"])])
def test_real_indel_support_reaches_variant_sequence_creation(name, record, alignments):
    variant = variant_from_record(record)
    with pysam.AlignmentFile(alignments[name]) as bam:
        evidence = ReadCollector(use_secondary_alignments=False).read_evidence_for_variant(variant, bam)
    assert len(evidence.alt_reads) >= 2
    sequences = VariantSequenceCreator(min_variant_sequence_coverage=1,
                                       variant_sequence_assembly=False).reads_to_variant_sequences(
                                           variant, evidence.alt_reads)
    assert sequences
    expected_alt = record["alt"][1:]
    assert all(s.alt == expected_alt and len(s) > 0 for s in sequences)
    assert all(s.read_names <= {r.name for r in evidence.alt_reads} for s in sequences)


def test_real_overlapping_mates_preserve_fragment_names_and_source_counts(alignments):
    # A genuine 26-base mate overlap outside the benchmark loci. This is a
    # fragment-assembly probe, not an assertion of a mutation at this position.
    with pysam.AlignmentFile(alignments["illumina_star"]) as bam:
        raw = ReadCollector(use_secondary_alignments=False).get_locus_reads(bam, "chr4", 15560, 15561)
        merged = ReadCollector(use_secondary_alignments=False, merge_overlapping_fragments=True).get_locus_reads(bam, "chr4", 15560, 15561)
    assert {r.name for r in raw} == {"SRR1258218.22918678"}
    assert {r.name for r in raw} == {r.name for r in merged}
    assert sum(r.source_read_count for r in merged) == len(raw)
    assert len(merged) == 1
    assert len(raw) == 2


@pytest.mark.parametrize("name", DATASETS)
def test_real_splice_skip_is_not_deletion_evidence(name, alignments):
    # Negative probe at a real splice junction, NOT a claimed DNA variant.
    with pysam.AlignmentFile(alignments[name]) as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            pos = read.reference_start
            for op, length in read.cigartuples:
                if op == 3:
                    # The independent oracle distinguishes N from D; both
                    # adjacent exon bases must exist for a meaningful probe.
                    anchors = read.get_reference_positions()
                    if pos - 1 in anchors and pos + length in anchors:
                        result = ReadCollector().locus_read_from_pysam_aligned_segment(read, pos, pos + length)
                        assert result is None, (read.query_name, pos, pos + length)
                        return
                if op in (0, 2, 3, 7, 8):
                    pos += length
    pytest.fail("Fixture did not contain a suitable real splice junction")
