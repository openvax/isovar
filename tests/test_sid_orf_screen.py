"""Sequence novelty and contributing reads must not imply full-ORF witnesses."""

from types import SimpleNamespace
from hashlib import sha256

import pytest
import pysam

from isovar import AlleleRead
from examples.screen_sid_non_snv_orfs import full_window_support
from examples.screen_sid_orf_candidates import proteome_screen
from examples.check_sid_coding_windows import check_window


def test_proteome_screen_distinguishes_exact_il_and_recombined_sequence(tmp_path):
    path = tmp_path / "proteins.fa"
    path.write_text(">p1\nMALKGHTQV\n>p2\nABCDEFGH\n>p3\nBCDEFGHK\n")
    result = proteome_screen(["MAIKGHTQV", "ABCDEFGHK", "MZZZZZZZZ"], path)
    assert result["proteins_scanned"] == 3
    variant = result["sequences"]["MAIKGHTQV"]
    assert variant["novel_exact_peptides"] and variant["exact_full_sequence_hits"] == []
    assert variant["novel_il_collapsed_peptides"] == []
    assert variant["il_collapsed_full_sequence_hits"] == ["p1"]
    recombined = result["sequences"]["ABCDEFGHK"]
    assert recombined["novel_exact_peptides"] == []
    assert recombined["exact_full_sequence_hits"] == []
    assert result["sequences"]["MZZZZZZZZ"]["novel_il_collapsed_peptides"]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_full_window_support_checks_allele_placement_quality_and_fragment_identity(strand):
    def read(name, prefix, allele, suffix, quality):
        return AlleleRead(prefix, allele, suffix, name, quality_scores=quality,
                         source_alignments=((("rg", name, 0), (0, 100, "9M", strand == "-")),))
    parts = ("ATG", "GCT", "TAA") if strand == "+" else ("TTA", "AGC", "CAT")
    complete = read("full", *parts, (30,) * 9)
    duplicate = read("full", *parts, (30,) * 9)
    missing = read("unknown", *parts, None)
    # Same nt occurs later, but not at the translation's focal allele offset.
    unrelated = read("unrelated", "CCC", "AAA", "ATGGCTTAA", (30,) * 15)
    partial = read("partial", "ATG", "GCT", "", (30,) * 6)
    translation = SimpleNamespace(
        variant_orf=SimpleNamespace(offset_to_first_complete_codon=0,
                                    cdna_sequence="ATGGCTTAA", variant_cdna_interval_start=3),
        reference_context=SimpleNamespace(strand=strand, mitochondrial=False),
        reads=[complete, duplicate, missing, unrelated, partial])
    protein = SimpleNamespace(translations=[translation], amino_acids="MA", ends_with_stop_codon=True)
    support = full_window_support(protein)
    assert support["full_window_fragments"] == 2
    assert support["q20_full_window_fragments"] == support["fragments_with_unavailable_qualities"] == 1
    protein.amino_acids = "MV"
    with pytest.raises(ValueError, match="translation"):
        full_window_support(protein)


def test_native_window_check_requires_original_nt_and_preserves_missing_quality(tmp_path):
    path = tmp_path / "original.bam"
    header = dict(SQ=[dict(SN="chr1", LN=1000)], RG=[dict(ID="rg")])
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for name, seq, qualities in [("full", "ATGGCTTAA", [30] * 9),
                                     ("unknown", "TTAAGCCAT", None)]:
            read = pysam.AlignedSegment(bam.header)
            read.query_name, read.query_sequence = name, seq
            read.reference_id, read.reference_start = 0, 100
            read.cigarstring, read.mapping_quality = "9M", 60
            read.query_qualities = qualities
            read.set_tag("RG", "rg")
            bam.write(read)
    nt_hash = sha256(b"ATGGCTTAA").hexdigest()
    protein = dict(amino_acids="MA", ends_with_stop_codon=True, mutation_interval=[1, 2],
        full_window_fragments=2, q20_full_window_fragments=1, witnesses=[dict(
            fragment_id=sha256(repr(("rg", name)).encode()).hexdigest(), nucleotide_sha256=nt_hash)
            for name in ("full", "unknown")])
    check = check_window(protein, path)
    assert check["full_window_fragments"] == 2 and check["native_q20_fragments"] == 1
    assert len(check["original_record_hashes"]) == 2
    protein["witnesses"][0]["nucleotide_sha256"] = sha256(b"ATGGTGTAA").hexdigest()
    with pytest.raises(ValueError, match="absent"):
        check_window(protein, path)


def test_full_window_witness_retains_the_nt_that_earned_quality_support():
    translations = []
    for suffix, quality in [("GCTTAA", 30), ("GCCTAA", 5)]:
        read = AlleleRead("", "ATG", suffix, "same_fragment", quality_scores=(quality,) * 9,
                         source_alignments=((("rg", "same_fragment", 0), (0, 100, "9M", False)),))
        translations.append(SimpleNamespace(
            variant_orf=SimpleNamespace(offset_to_first_complete_codon=0,
                                        cdna_sequence=read.sequence, variant_cdna_interval_start=0),
            reference_context=SimpleNamespace(strand="+", mitochondrial=False), reads=[read]))
    protein = SimpleNamespace(translations=translations, amino_acids="MA", ends_with_stop_codon=True)
    support = full_window_support(protein)
    assert support["full_window_fragments"] == support["q20_full_window_fragments"] == 1
    assert support["witnesses"][0]["nucleotide_sha256"] == sha256(b"ATGGCTTAA").hexdigest()
