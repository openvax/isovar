"""Original-read regression for the unresolved COLO829 SV3 protein hypotheses."""

import json

from Bio import SeqIO

from reports.hcc1395_colo829.build_colo829_sv3_protein_audit import write_audit


def test_sv3_junction_orfs_remain_distinct_and_do_not_become_a_proteoform(tmp_path):
    prefix = tmp_path / "sv3_orfs"
    assessment_path = write_audit(prefix)
    assessment = json.loads(assessment_path.read_text())
    assert assessment["event_id"] == "COLO829_SOMATIC_SV3_DEL"
    assert assessment["proteoform_call"] == "not_established"
    assert assessment["path_count"] == 1
    assert assessment["flanking_voting_fragments"] == 533
    assert assessment["unplaced_junction_bases"] == "CTATAGTAGAACTTCTTTT"
    assert assessment["original_junction_read_names"] == [
        "SRR8615617.101713327", "SRR8615617.7313798", "SRR8615617.97504216"]
    assert assessment["junction_base_quality_min_by_read"] == {
        "SRR8615617.101713327": 34,
        "SRR8615617.7313798": 38,
        "SRR8615617.97504216": 23,
    }
    assert assessment["locus"]["genes"] == ["CNIH3"]
    assert assessment["locus"]["overlapping_ensembl87_exons"] == []
    assert assessment["locus"]["overlapping_bundled_orfdata_orfs"] == []

    candidates = assessment["candidates"]
    assert [row["amino_acid_length"] for row in candidates] == [37, 162]
    assert {row["genomic_start_zero_based"] for row in candidates} == {224595019, 224595032}
    assert all(row["shared_junction_fragments"] == 3 for row in candidates)
    assert all(row["full_orf_junction_fragments"] == 0 for row in candidates)
    assert all(not row["annotated_start"] and not row["initiation_observed"] for row in candidates)
    assert all(not row["translation_observed"] for row in candidates)
    assert all(row["protein_status"] == "exploratory_orf_mature_transcript_unresolved"
               for row in candidates)
    assert [row["shared_reference_genomic_prefix_amino_acids"] for row in candidates] == [25, 20]
    assert [row["junction_ninemers"] for row in candidates] == [
        ["QAYSRTSFG", "AYSRTSFGV"],
        ["VRPIVELLL", "RPIVELLLE", "PIVELLLEL"],
    ]
    assert all(not matches for row in candidates for matches in
               row["ensembl87_protein_matches_by_junction_ninemer"].values())

    fasta = list(SeqIO.parse(str(prefix) + ".protein.fasta", "fasta"))
    assert {record.id: str(record.seq) for record in fasta} == {
        row["candidate_id"]: row["amino_acids"] for row in candidates}
    assert b"\r" not in (tmp_path / "sv3_orfs.tsv").read_bytes()
    report = (tmp_path / "sv3_orfs.assessment.md").read_text()
    assert "**Proteoform call: not established.**" in report
    assert "no single read is required" in report
