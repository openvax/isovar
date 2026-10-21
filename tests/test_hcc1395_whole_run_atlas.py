"""Offline checks of HCC1395 source and original direct-RNA evidence."""

import csv
import gzip
from pathlib import Path

import pysam
from pyensembl import cached_release

from isovar import ProteinSequenceCreator, ReadCollector
from reports.hcc1395_colo829.build_atlas_regions import merged_regions
from reports.hcc1395_colo829.build_colo829_atlas import evaluate_call, source_calls
from reports.hcc1395_colo829.build_hcc1395_regions import regions
from reports.hcc1395_colo829.build_hcc1395_sv_screen import screen_rows


ROOT = Path(__file__).parents[1]
REPORT = ROOT / "reports/hcc1395_colo829"
DATA = ROOT / "tests/data/hcc1395_colo829_original"
FULL_DATA = ROOT / "tests/data/hcc1395_whole_run"


def test_region_builder_covers_all_hcc_source_records():
    selected = list(regions(REPORT / "somatic_call_inventory.tsv.gz"))
    small_calls = list(source_calls(REPORT / "somatic_call_inventory.tsv.gz",
                                    sample="HCC1395", dna_source="SEQC2_v1.2.1"))
    assert len(small_calls) == 41072
    assert len(selected) == 41072 + 2 * 1788
    assert len(list(merged_regions(selected))) < len(selected)
    prdx5 = next(row for row in small_calls if row["call_id"] == "COSM32759")
    assert (prdx5["chrom1"], int(prdx5["pos1"]) - 1,
            int(prdx5["pos1"])) in selected


def test_original_direct_rna_alt_is_observed_without_asserting_a_protein():
    source = next(row for row in source_calls(
        REPORT / "somatic_call_inventory.tsv.gz", sample="HCC1395",
        dna_source="SEQC2_v1.2.1") if row["call_id"] == "COSM32759")
    with pysam.AlignmentFile(str(DATA / "hcc1395_srr23456709.bam"), "rb") as alignment:
        call, proteins, witnesses = evaluate_call(
            source, alignment, cached_release(87), ReadCollector(),
            ProteinSequenceCreator(max_protein_sequences_per_variant=0),
            rna_run="SRR23456709")
    assert call["rna_run"] == "SRR23456709"
    assert (call["ref_reads"], call["alt_reads"], call["other_reads"]) == (40, 18, 14)
    assert call["rna_status"] == "alternate_observed"
    assert call["protein_status"] == "no_translation"
    assert len(witnesses) == 18
    assert not proteins


def test_full_run_fixture_retains_local_protein_and_allele_statuses():
    chosen = {
        "chr1:877772:G>C", "chr2:46384706:C>G", "COSM32736", "COSM32759",
        "chr8:123220221:A>AGGGCTAC",
    }
    sources = source_calls(REPORT / "somatic_call_inventory.tsv.gz", chosen,
                           "HCC1395", "SEQC2_v1.2.1")
    observed = {}
    with pysam.AlignmentFile(str(FULL_DATA / "SRR23456709.fixture.grch38.bam"), "rb") as bam:
        for source in sources:
            observed[source["call_id"]] = evaluate_call(
                source, bam, cached_release(87), ReadCollector(),
                ProteinSequenceCreator(max_protein_sequences_per_variant=0),
                rna_run="SRR23456709")
    assert observed["chr1:877772:G>C"][0]["rna_status"] == "no_alignment_overlap"
    assert observed["chr2:46384706:C>G"][0]["rna_status"] == "reference_only"
    assert observed["chr2:46384706:C>G"][0]["ref_reads"] == 28
    ostc, ostc_proteins, ostc_reads = observed["COSM32736"]
    assert (ostc["alt_reads"], len(ostc_reads)) == (60, 60)
    assert ostc_proteins[0]["protein_sequence"] == "METLYRVPLLVLECPNLKLKKPPWL"
    assert ostc_proteins[0]["mutation_start_idx"] == 8
    assert ostc_proteins[0]["supporting_fragments"] == 5
    prdx5, prdx5_proteins, _ = observed["COSM32759"]
    assert prdx5["alt_reads"] == 193
    assert prdx5_proteins[0]["protein_sequence"] == "VRLLADPTGALGKETDLLLD"
    assert prdx5_proteins[0]["supporting_fragments"] == 6
    indel, indel_proteins, _ = observed["chr8:123220221:A>AGGGCTAC"]
    assert indel["alt_reads"] == 3
    assert indel_proteins[0]["mutation_end_idx"] == 25
    assert indel_proteins[0]["supporting_fragments"] == 2


def test_hcc_sv_split_candidate_does_not_claim_an_rna_path():
    with pysam.AlignmentFile(str(FULL_DATA / "SRR23456709.fixture.grch38.bam"), "rb") as bam:
        rows = list(screen_rows(REPORT / "somatic_call_inventory.tsv.gz", bam))
    row = next(row for row in rows
               if row["call_id"] == "chr9_83969153_chr9_84015113_-+")
    assert row["observed_split_alignment_candidates"] == 7
    assert row["rna_status"] == "split_alignment_candidate_unresolved"
    assert row["path_status"] == "unresolved"
    assert row["protein_status"] == "not_called_unresolved_path"
    fusion = next(row for row in rows if row["published_rna_fusion"] == "POU2F1-NME7")
    assert fusion["observed_split_alignment_candidates"] == 1
    assert fusion["candidate_read_ids"] == "SRR23456709.320596"
    assert fusion["protein_status"] == "not_called_unresolved_path"


def test_full_run_tables_cover_every_source_call_and_alt_observation():
    expected = {row["call_id"] for row in source_calls(
        REPORT / "somatic_call_inventory.tsv.gz", sample="HCC1395",
        dna_source="SEQC2_v1.2.1")}
    with gzip.open(REPORT / "hcc1395_direct_rna_atlas_calls.tsv.gz", "rt") as source:
        calls = list(csv.DictReader(source, delimiter="\t"))
    with gzip.open(REPORT / "hcc1395_direct_rna_atlas_alt_witnesses.tsv.gz", "rt") as source:
        witnesses = list(csv.DictReader(source, delimiter="\t"))
    with gzip.open(REPORT / "hcc1395_direct_rna_atlas_proteins.tsv.gz", "rt") as source:
        proteins = list(csv.DictReader(source, delimiter="\t"))
    assert len(calls) == len(expected) == 41072
    assert {row["call_id"] for row in calls} == expected
    assert sum(int(row["alt_reads"]) for row in calls) == len(witnesses) == 2681
    assert sum(int(row["protein_count"]) for row in calls) == len(proteins) == 119
    assert all(row["call_id"] in expected for row in witnesses + proteins)
