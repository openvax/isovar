"""Replay public original RNA reads at matched-normal somatic sites."""

from pathlib import Path

import pysam
from varcode import load_vcf

from isovar import ProteinSequenceCreator, ReadCollector

from .genomes_for_testing import grch38


DATA = Path(__file__).parent / "data/hcc1395_colo829_original"


def test_colo829_naxd_original_rna_mutant_protein():
    # SMaHT COLO829/BL SNV, with original CCLE RNA pairs from SRR8615617.
    variants = load_vcf(str(DATA / "colo829_selected.vcf"), genome=grch38)
    variant = next(v for v in variants if (v.contig, v.start) == ("13", 110634738))
    assert (variant.ref, variant.alt) == ("C", "T")
    with pysam.AlignmentFile(str(DATA / "colo829_srr8615617.bam"), "rb") as alignment:
        evidence = ReadCollector().read_evidence_for_variant(variant, alignment)
    assert (len(evidence.ref_reads), len(evidence.alt_reads), len(evidence.other_reads)) == (104, 81, 0)
    assert sum(read.source_read_count for read in evidence.alt_reads) == 94

    proteins = ProteinSequenceCreator().sorted_protein_sequences_for_variant(variant, evidence)
    assert len(proteins) == 1
    protein = proteins[0]
    assert protein.amino_acids == "LWLVAQQPALIHGYRKAVLTSNHVEFSRLYDAV"
    assert (protein.mutation_start_idx, protein.mutation_end_idx) == (20, 21)
    assert protein.amino_acids[20:21] == "S"
    assert protein.num_supporting_fragments == 59


def test_hcc1395_prdx5_direct_rna_support_without_default_translation():
    # SEQC2 HCC1395/BL SNV, with original direct RNA from SRR23456709.
    variants = load_vcf(str(DATA / "hcc1395_selected.vcf"), genome=grch38,
                        only_passing=False)
    variant = next(v for v in variants if (v.contig, v.start) == ("11", 64320895))
    assert (variant.ref, variant.alt) == ("T", "C")
    with pysam.AlignmentFile(str(DATA / "hcc1395_srr23456709.bam"), "rb") as alignment:
        evidence = ReadCollector().read_evidence_for_variant(variant, alignment)
    assert (len(evidence.ref_reads), len(evidence.alt_reads), len(evidence.other_reads)) == (40, 18, 14)
    assert sum(read.source_read_count for read in evidence.alt_reads) == 18
    assert ProteinSequenceCreator().sorted_protein_sequences_for_variant(variant, evidence) == []


def test_hcc1395_prdx5_cdna_mutant_protein():
    # The independent cDNA assay translates this same SEQC2 truth-set allele.
    variants = load_vcf(str(DATA / "hcc1395_selected.vcf"), genome=grch38,
                        only_passing=False)
    variant = next(v for v in variants if (v.contig, v.start) == ("11", 64320895))
    with pysam.AlignmentFile(str(DATA / "hcc1395_srr23456710.bam"), "rb") as alignment:
        evidence = ReadCollector().read_evidence_for_variant(variant, alignment)
    assert (len(evidence.ref_reads), len(evidence.alt_reads), len(evidence.other_reads)) == (41, 24, 2)

    proteins = ProteinSequenceCreator().sorted_protein_sequences_for_variant(variant, evidence)
    assert len(proteins) == 1
    protein = proteins[0]
    assert protein.amino_acids == "GKVRLLADPTGALGKETDLLLDDSL"
    assert (protein.mutation_start_idx, protein.mutation_end_idx) == (12, 13)
    assert protein.amino_acids[12:13] == "L"
    assert protein.num_supporting_fragments == 6
