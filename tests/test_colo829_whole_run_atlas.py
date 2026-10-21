"""Offline checks of whole-run atlas logic using exact public RNA read fixtures."""

import csv
import gzip
from pathlib import Path

import pysam
from pyensembl import cached_release

from isovar import ProteinSequenceCreator, ReadCollector
from reports.hcc1395_colo829.build_atlas_regions import breakpoints, merged_regions, regions
from reports.hcc1395_colo829.build_colo829_atlas import evaluate_call, source_calls
from reports.hcc1395_colo829.build_colo829_sv_screen import screen_rows


ROOT = Path(__file__).parents[1]
REPORT = ROOT / "reports/hcc1395_colo829"
DATA = ROOT / "tests/data/hcc1395_colo829_original"
FULL_DATA = ROOT / "tests/data/colo829_whole_run"


def test_source_allele_and_local_protein_are_linked_to_original_reads():
    source = next(row for row in source_calls(REPORT / "somatic_call_inventory.tsv.gz")
                  if row["call_id"] == "chr13:110634738:C>T")
    with pysam.AlignmentFile(str(DATA / "colo829_srr8615617.bam"), "rb") as alignment:
        call, proteins, witnesses = evaluate_call(
            source, alignment, cached_release(87), ReadCollector(),
            ProteinSequenceCreator(max_protein_sequences_per_variant=0))
    assert (call["ref_reads"], call["alt_reads"], call["other_reads"]) == (104, 81, 0)
    assert (call["ref_fragments"], call["alt_fragments"]) == (100, 78)
    assert call["alt_source_reads"] == 94
    assert call["rna_status"] == "alternate_observed"
    assert call["assembly_status"] == "default"
    assert len(witnesses) == 81
    assert proteins[0]["protein_sequence"] == "LWLVAQQPALIHGYRKAVLTSNHVEFSRLYDAV"
    assert proteins[0]["supporting_fragments"] == 59
    assert proteins[0]["sequence_scope"] == "local_context"


def test_no_coverage_does_not_mean_reference_or_translation():
    source = next(row for row in source_calls(REPORT / "somatic_call_inventory.tsv.gz")
                  if row["call_id"] == "chr1:183180:G>A")
    with pysam.AlignmentFile(str(DATA / "colo829_srr8615617.bam"), "rb") as alignment:
        call, proteins, witnesses = evaluate_call(
            source, alignment, cached_release(87), ReadCollector(),
            ProteinSequenceCreator(max_protein_sequences_per_variant=0))
    assert call["rna_status"] == "no_alignment_overlap"
    assert call["raw_alignment_overlap"] == 0
    assert call["protein_status"] == "not_attempted_no_alt"
    assert not proteins and not witnesses


def test_target_regions_include_indel_bases_and_lifted_sv_partner():
    lifted = ["14", "104093751", "truthset_50_1", "G", "G[14:104093806[",
              ".", "PASS", "SVTYPE=DEL"]
    assert breakpoints(lifted) == (("14", 104093751), ("14", 104093806))
    source_indel = next(row for row in source_calls(REPORT / "somatic_call_inventory.tsv.gz")
                        if row["class"] == "indel")
    selected = list(regions(REPORT / "somatic_call_inventory.tsv.gz", ()))
    start = int(source_indel["pos1"]) - 1
    assert (source_indel["chrom1"], start,
            start + max(1, len(source_indel["ref"]))) in selected
    assert list(merged_regions([("chr1", 2, 4), ("chr1", 4, 5)])) == [("chr1", 2, 5)]


def test_sv_screen_keeps_sources_separate_and_withholds_unresolved_protein(tmp_path):
    v2 = tmp_path / "v2.vcf"
    lifted = tmp_path / "lifted.vcf"
    v2.write_text("chr1\t207807886\tCOLO829_SOMATIC_SV1_DEL\tN\t<DEL>\t.\tPASS\t"
                  "SVTYPE=DEL;END=207841474;CHR2=chr1\n")
    lifted.write_text("14\t104093751\ttruthset_50_1\tG\tG[14:104093806[\t.\tPASS\t"
                      "SVTYPE=DEL\n")
    with pysam.AlignmentFile(str(DATA / "colo829_srr8615617.bam"), "rb") as alignment:
        rows = list(screen_rows((v2, lifted), alignment))
    assert len(rows) == 2
    assert [row["dna_source"] for row in rows] == [
        "PaulinPaz_2024_v2", "ValleInclan_2022_lifted"]
    assert rows[1]["chrom2"] == "14" and rows[1]["pos2"] == 104093806
    assert rows[1]["vcf_alt"] == "G[14:104093806["
    assert all(row["rna_status"] == "no_breakend_window_alignment" for row in rows)
    assert all(row["protein_status"] == "not_called_unresolved_path" for row in rows)


def test_full_run_sv_mate_context_is_not_a_junction_or_protein_call():
    vcfs = (FULL_DATA / "colo829-sv-grch38.vcf",
            FULL_DATA / "colo829-sv-2022-hg38lifted.vcf")
    with pysam.AlignmentFile(str(FULL_DATA / "SRR8615617.all_sites.grch38.bam"), "rb") as alignment:
        row = next(row for row in screen_rows(vcfs, alignment)
                   if row["call_id"] == "COLO829_SOMATIC_SV3_DEL")
    assert (row["end1_fragments"], row["end2_fragments"]) == (427, 956)
    assert row["cross_end_mate_pair_candidates"] == 6
    assert row["observed_split_alignment_candidates"] == 0
    assert row["rna_status"] == "cross_end_pair_context_unresolved"
    assert row["protein_status"] == "not_called_unresolved_path"


def test_deep_indel_uses_all_alt_reads_without_overlap_assembly():
    source = next(row for row in source_calls(REPORT / "somatic_call_inventory.tsv.gz")
                  if row["call_id"] == "chr9:21971154:CCG>C")
    with pysam.AlignmentFile(str(FULL_DATA / "SRR8615617.all_sites.grch38.bam"), "rb") as alignment:
        call, proteins, witnesses = evaluate_call(
            source, alignment, cached_release(87), ReadCollector(),
            ProteinSequenceCreator(max_protein_sequences_per_variant=0))
    assert (call["ref_reads"], call["alt_reads"]) == (552, 6984)
    assert len(witnesses) == 6984
    assert call["assembly_status"] == "direct_read_contexts_high_depth"
    assert call["protein_status"] == "multiple_local_contexts"
    assert call["protein_count"] == len(proteins) > 0
    assert proteins[0]["protein_sequence"] == "MMMGSARVAELLLLHGGAQLRRPR"
    assert proteins[0]["supporting_reads"] == 5338


def test_atlas_tables_cover_every_source_call_and_alt_observation():
    expected_ids = {row["call_id"] for row in source_calls(
        REPORT / "somatic_call_inventory.tsv.gz")}
    with gzip.open(REPORT / "colo829_atlas_calls.tsv.gz", "rt") as source:
        calls = list(csv.DictReader(source, delimiter="\t"))
    with gzip.open(REPORT / "colo829_atlas_alt_witnesses.tsv.gz", "rt") as source:
        witnesses = list(csv.DictReader(source, delimiter="\t"))
    with gzip.open(REPORT / "colo829_atlas_proteins.tsv.gz", "rt") as source:
        proteins = list(csv.DictReader(source, delimiter="\t"))
    assert len(calls) == len(expected_ids) == 46064
    assert {row["call_id"] for row in calls} == expected_ids
    assert sum(int(row["alt_reads"]) for row in calls) == len(witnesses)
    assert sum(int(row["protein_count"]) for row in calls) == len(proteins)
    assert all(row["protein_status"] != "not_attempted_read_limit" for row in calls)
    assert all(row["call_id"] in expected_ids for row in witnesses + proteins)
