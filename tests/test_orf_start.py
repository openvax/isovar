"""#350 synthetic builders: explicit transcript sequence, exon map and RNA path.

No clinical sequence is invented here. Original-read regressions live in
test_sv_rna_orfs and test_sv_rna_three_fusions with their pinned manifests.
"""

from dataclasses import replace
import csv
import json

import pysam
import pytest

from isovar import (
    FusionBreakpoint, FusionReference, annotate_orf_start, export_sv_rna_orfs,
    reconstruct_sv_rna, summarize_orf_start_evidence, write_sv_rna_orfs,
)
from tests.test_sv_rna import aligned, write_bam


def reference(leader="CCCATGAAATAACCC", strand="+"):
    """Two exons, one intron, a complete main CDS and a separate leader ATG."""
    coding = "ATG" + "GCC" * 8 + "TAA"
    sequence = leader + coding + "CCC"
    split = len(leader)
    exons = ((100, 100 + split), (300, 300 + len(sequence) - split))
    if strand == "-":
        exons = ((100, 100 + len(sequence) - split), (300, 300 + split))
    return FusionReference("tx1", "synthetic", "explicit-fixture-v1", "1", strand, exons,
                           sequence, len(leader), len(leader) + len(coding))


def transcript_positions(ref):
    positions = [(ref.contig, p, ref.strand) for a, b in ref.exons for p in range(a, b)]
    return positions if ref.strand == "+" else positions[::-1]


def annotation(ref, offset=3):
    return annotate_orf_start(ref.sequence, transcript_positions(ref), offset, [ref])


@pytest.mark.parametrize("strand", ["+", "-"])
def test_known_reference_start_utr_and_intronic_origins_are_distinct(strand):
    ref = reference(strand=strand)
    utr = annotation(ref)
    coding = annotation(ref, ref.cds_start)
    intronic_positions = [("1", p, strand) for p in ([200, 201, 202] if strand == "+" else [202, 201, 200])]
    intronic = annotate_orf_start("ATG", intronic_positions, 0, [ref])
    assert [coding["priority"], utr["priority"], intronic["priority"]] == [1, 2, 4]
    assert utr["assessments"][0]["reference_orf_type"] == "upstream_ORF"
    row, = intronic["assessments"]
    assert row["region"] == "intronic" and row["transcript_inclusion"] == "not_established"
    assert row["splice_inference"] == "not_assessed"
    assert all(not a["initiation_observed"] and not a["translation_observed"] for a in (coding, utr, intronic))


@pytest.mark.parametrize("leader,expected", [
    ("CCCATGAAATAACCC", "upstream_ORF"),
    ("CCCATGCCC", "in_frame_extension"),
    ("CCCATGCCCC", "stop_unresolved"),
    ("CCCATGCC", "stop_unresolved"),
])
def test_reference_utr_reading_is_not_automatically_a_uorf(leader, expected):
    evidence = annotation(reference(leader))
    assert evidence["assessments"][0]["reference_orf_type"] == expected
    assert evidence["tier"] == "five_prime_UTR_ATG"


def test_reference_overlapping_orf_requires_an_observed_reference_stop():
    ref = reference("CCCATGCCCC")
    ref = replace(ref, sequence=ref.sequence[:-3] + "CCTAA",
                  exons=(ref.exons[0], (ref.exons[1][0], ref.exons[1][1] + 2)))
    assert annotation(ref)["assessments"][0]["reference_orf_type"] == "overlapping_ORF"


def test_overlapping_isoforms_do_not_inherit_the_strongest_tier():
    utr = reference()
    coding = replace(utr, transcript_id="tx2", cds_start=3, cds_end=12)
    evidence = annotate_orf_start(utr.sequence, transcript_positions(utr), 3, [utr, coding])
    assert evidence["status"] == "ambiguous" and evidence["priority"] is None
    assert [r["priority"] for r in evidence["assessments"]] == [2, 1]
    assert evidence == annotate_orf_start(utr.sequence, transcript_positions(utr), 3, [coding, utr])
    assert summarize_orf_start_evidence([evidence, None])["status"] == "unavailable"


@pytest.mark.parametrize("positions,region", [
    ([None] * 3, "unplaced"),
    ([("1", 104, "+"), None, None], "partially_unplaced"),
    ([("1", 105, "-"), ("1", 104, "-"), ("1", 103, "-")], "antisense"),
    ([("1", 103, "+"), ("1", 105, "+"), ("1", 106, "+")], "junction_created"),
    ([("1", 103, "+"), ("2", 200, "+"), ("2", 201, "+")], "junction_created"),
    ([("3", p, "+") for p in range(3)], "outside_supplied_transcripts"),
])
def test_sequence_only_subtypes_remain_explicit(positions, region):
    evidence = annotate_orf_start("ATG", positions, 0, [reference()])
    assert evidence["tier"] == "sequence_only" and evidence["priority"] == 4
    assert evidence["assessments"][0]["region"] == region
    assert not evidence["reference_comparisons"]


def test_start_codon_must_match_all_three_spliced_positions_and_reference_bases():
    ref = reference()
    positions = transcript_positions(ref)
    # A reference transcript containing an ATG at the same offset is not enough
    # if even one placed base belongs to another locus.
    positions[4] = ("2", 500, "+")
    assert annotate_orf_start(ref.sequence, positions, 3, [ref])["priority"] == 4
    # A new UTR ATG still has a UTR location but is not a reference uORF.
    altered = replace(ref, sequence=ref.sequence[:3] + "ACG" + ref.sequence[6:])
    evidence = annotate_orf_start(ref.sequence, transcript_positions(ref), 3, [altered])
    assert evidence["priority"] == 2 and evidence["reference_comparisons"] == []
    assert evidence["assessments"][0]["reference_start_codon"] == "ACG"
    assert evidence["assessments"][0]["reference_orf_type"] is None


def test_start_can_cross_an_annotated_exon_boundary():
    ref = reference()
    # Split inside the leader ATG, while preserving the transcript sequence.
    ref = replace(ref, exons=((100, 104), (300, 300 + len(ref.sequence) - 4)))
    evidence = annotation(ref)
    assert evidence["priority"] == 2 and evidence["reference_comparisons"][0]["transcript_offset"] == 3


@pytest.mark.parametrize("sequence,positions,start", [
    ("CCC", [None] * 3, 0), ("ATG", [], 0), ("ATG", [None] * 3, -1),
    ("ATG", [None] * 3, True), ("ATG", [None] * 3, 1),
])
def test_invalid_atg_input_is_rejected(sequence, positions, start):
    with pytest.raises(ValueError, match="ATG offset"):
        annotate_orf_start(sequence, positions, start, [])


def reconstructed_start(tmp_path, *, intronic, assemble):
    """Equal RNA support and peptide; move the start from leader into intron."""
    donor_sequence = "CCCATG" + "GCC" * 6
    sequence = donor_sequence + "GAA" * 3 + "TAACCCCCCCCC"
    ref = reference(donor_sequence + "TAACCCCCCCCC")
    genomic_start = 200 if intronic else 100
    positions = [("1", genomic_start + q, "+") for q in range(len(donor_sequence))]
    positions += [("2", 500 + q, "+") for q in range(len(sequence) - len(donor_sequence))]
    records = [r for i in range(3) for r in aligned("fragment%d" % i, sequence, positions)]
    bam_path = write_bam(tmp_path / ("intronic.bam" if intronic else "utr.bam"), records)
    with pysam.AlignmentFile(str(bam_path)) as bam:
        return reconstruct_sv_rna(
            bam, event_id="synthetic-350", reference_name="synthetic", references=[ref],
            donor=FusionBreakpoint("1", genomic_start + len(donor_sequence), "+"),
            acceptor=FusionBreakpoint("2", 500, "+"), sample_id="sample", source="synthetic",
            regions=(), event_provenance=dict(fixture="test_orf_start.reconstructed_start"),
            min_orf_amino_acids=1, assemble=assemble)


@pytest.mark.parametrize("assemble", [True, False])
def test_reconstruction_through_export_ranks_equal_support_utr_above_intronic(tmp_path, assemble):
    exported = [export_sv_rna_orfs(reconstructed_start(tmp_path, intronic=i, assemble=assemble))
                for i in (True, False)]
    candidates = [e["candidates"][0] for e in exported]
    assert all(len(e["candidates"]) == 1 for e in exported)
    assert candidates[0]["amino_acids"] == candidates[1]["amino_acids"] == "MAAAAAAEEE"
    assert [c["rna_support"]["fragments"] for c in candidates] == [3, 3]
    assert [c["start_evidence_summary"]["priority"] for c in candidates] == [4, 2]
    ranked = sorted(candidates, key=lambda c: c["start_evidence_summary"]["priority"])
    assert ranked[0]["start_evidence_summary"]["tier"] == "five_prime_UTR_ATG"
    # Annotation and reconstruction mode do not alter sequence identity.
    assert candidates[0]["candidate_id"] == candidates[1]["candidate_id"]
    paths = write_sv_rna_orfs(exported[1], tmp_path / "ranked")
    with paths["tsv"].open() as handle:
        row, = csv.DictReader(handle, delimiter="\t")
    assert row["start_priority"] == "2" and row["start_tier"] == "five_prime_UTR_ATG"
    metadata, = json.loads(row["start_evidence_json"])
    assert metadata["start_evidence"]["assessments"][0]["annotation"] == "explicit-fixture-v1"
