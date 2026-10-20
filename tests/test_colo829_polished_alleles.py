"""Sequence-level SNV calls must not treat absent coverage as reference."""

import csv
import gzip
from pathlib import Path

import pysam
import pytest

from reports.hcc1395_colo829.build_polished_alleles import (
    observations,
    read_counts,
    read_variants,
    write_outputs,
)


def test_polished_isoform_alleles_require_a_callable_base(tmp_path):
    fasta_path = tmp_path / "reference.fa"
    fasta_path.write_text(">chr1\n" + "A" * 30 + "\n")
    pysam.faidx(str(fasta_path))
    vcf_path = tmp_path / "snvs.vcf.gz"
    with gzip.open(vcf_path, "wt") as target:
        target.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        target.write("chr1\t3\t.\tA\tC\t.\tPASS\t.\n")
        target.write("chr1\t8\t.\tA\tG\t.\tPASS\t.\n")
    counts_path = tmp_path / "counts.csv"
    counts_path.write_text("id,COLO829T,COLO829BL\n"
                           "PB.alt,4,0\nPB.ref,2,3\nPB.other,1,0\n"
                           "PB.distant,8,0\nPB.low_mapq,9,0\nPB.secondary,9,0\n")

    bam_path = tmp_path / "polished.bam"
    with pysam.AlignmentFile(bam_path, "wb", header={
        "HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 30}],
    }) as target:
        for name, sequence, start, mapq, flag in (
            ("PB.alt", "AACAA", 0, 60, 0),
            ("PB.ref", "AAAAA", 0, 60, 0),
            ("PB.other", "AAGAA", 0, 60, 0),
            ("PB.low_mapq", "AACAA", 0, 0, 0),
            ("PB.secondary", "AACAA", 0, 60, 256),
            ("PB.distant", "AAAAA", 12, 60, 0),
        ):
            read = pysam.AlignedSegment()
            read.query_name = name + "|source-header"
            read.query_sequence = sequence
            read.flag = flag
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = mapq
            read.cigarstring = "5M"
            target.write(read)

    with pysam.FastaFile(fasta_path) as reference:
        variants = read_variants(vcf_path, reference)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        rows = list(observations(bam, variants, read_counts(counts_path)))
    assert {(row["isoform_id"], row["allele"]) for row in rows} == {
        ("PB.alt", "1"), ("PB.ref", "0"), ("PB.other", "?"),
    }
    assert all(row["variant_id"] == "chr1:3:A>C" for row in rows)
    assert write_outputs(rows, variants, tmp_path) == (3, 2)
    with gzip.open(tmp_path / "colo829_polished_snv_summary.tsv.gz", "rt") as source:
        summary = list(csv.DictReader(source, delimiter="\t"))
    assert summary[0]["tumor_alt_isoforms"] == "1"
    assert summary[0]["tumor_ref_isoforms"] == "1"
    assert summary[0]["tumor_other_isoforms"] == "1"
    assert summary[0]["tumor_alt_fl_count"] == "4"
    assert summary[0]["normal_ref_fl_count"] == "3"
    assert summary[1]["tumor_ref_isoforms"] == "0"
    assert summary[1]["tumor_callability"] == "not_callable"


def test_rejects_a_variant_with_the_wrong_reference_build(tmp_path):
    fasta_path = tmp_path / "reference.fa"
    fasta_path.write_text(">chr1\nAAAAA\n")
    pysam.faidx(str(fasta_path))
    vcf_path = tmp_path / "snvs.vcf.gz"
    with gzip.open(vcf_path, "wt") as target:
        target.write("chr1\t3\t.\tG\tC\t.\tPASS\t.\n")
    with pysam.FastaFile(fasta_path) as reference:
        with pytest.raises(ValueError, match="Reference mismatch"):
            read_variants(vcf_path, reference)


def test_real_colo829_representatives_have_reference_and_alt_bases():
    """Three published TAGET representatives aligned at one SMaHT SNV."""
    fixture = Path(__file__).parent / "data" / "colo829_polished" / "representatives.bam"
    site = 948136 - 1
    variants = {"chr1": ([site], {site: ("G", "A")})}
    counts = {"PB.7.4": (68, 99), "PB.7.5": (2, 0), "PB.7.8": (7, 0)}
    with pysam.AlignmentFile(fixture, "rb") as bam:
        rows = list(observations(bam, variants, counts))
    assert {(row["isoform_id"], row["allele"], row["observed_base"]) for row in rows} == {
        ("PB.7.4", "0", "G"), ("PB.7.5", "1", "A"),
        ("PB.7.8", "1", "A"),
    }
    assert all(row["mapping_quality"] == 60 for row in rows)
