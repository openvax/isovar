"""Orientation warnings preserve the public reconstruction evidence contract."""

from copy import deepcopy
import csv
import json

import pysam
import pytest

from isovar import FusionBreakpoint, export_sv_rna_orfs, reconstruct_sv_rna, write_sv_rna_orfs
from tests.test_sv_rna import aligned, antisense, write_bam


@pytest.mark.parametrize("mode,counts,flag", [
    ("original", (2, 0, 0), None),
    ("reverse", (0, 2, 0), "reverse_complement_support_only"),
    ("both", (1, 1, 0), "mixed_read_orientations"),
    ("within_fragment", (0, 0, 1), "mixed_read_orientations"),
    ("absent", (0, 0, 0), "no_full_fragment_witness"),
])
def test_reconstruction_export_json_tsv_orientation_contract(tmp_path, mode, counts, flag):
    sequence = "CCC" + "ATG" + "GCT" * 16 + "GAA" * 5 + "TAA" + "CCC"
    positions = [("1", 1000 + q, "+") for q in range(30)] + [
        ("2", 2000 + q, "+") for q in range(len(sequence) - 30)]
    records = []
    for i in range(2):
        seq, pos = sequence, positions
        if mode == "absent":
            lo, hi = (0, 54) if i == 0 else (18, len(sequence))
            seq, pos = seq[lo:hi], pos[lo:hi]
        elif mode == "reverse" or mode in ("both", "within_fragment") and i == 1:
            seq, pos = antisense(seq, pos)
        name = "paired" if mode == "within_fragment" else "fragment%d" % i
        bits = (64 if i == 0 else 128) if mode == "within_fragment" else 0
        records.extend(aligned(name, seq, pos, flag=bits))
    # Duplicate observations do not create new templates or orientation votes.
    bam_path = write_bam(tmp_path / "input.bam", records + records)
    with pysam.AlignmentFile(bam_path) as bam:
        result = reconstruct_sv_rna(
            bam, event_id="synthetic-orientations", reference_name="synthetic", references=[],
            donor=FusionBreakpoint("1", 1030, "+"), acceptor=FusionBreakpoint("2", 2000, "+"),
            regions=(), sample_id="sample", source="library", event_provenance={"caller": "synthetic"},
            min_overlap=6, min_alternative_fragments=1, min_orf_amino_acids=1)
    export = export_sv_rna_orfs(result)
    candidate, = [c for c in export["candidates"] if c["nucleotide_sequence"] == sequence[3:-3]]
    support = candidate["rna_support"]
    assert support["fragment_query_orientations"] == dict(zip(
        ("original_query", "reverse_complement", "mixed"), counts))
    assert support["fragments"] == sum(counts)
    flags = set(candidate["uncertainty_flags"])
    assert "rna_strand_unresolved" in flags
    assert ("reverse_complement_support_only" in flags) == (mode == "reverse")
    assert ("mixed_read_orientations" in flags) == (mode in ("both", "within_fragment"))
    if flag:
        assert flag in flags
    paths = write_sv_rna_orfs(export, tmp_path / "output")
    assert json.loads(paths["json"].read_text()) == json.loads(json.dumps(export))
    with paths["tsv"].open() as handle:
        row, = [r for r in csv.DictReader(handle, delimiter="\t")
                if r["candidate_id"] == candidate["candidate_id"]]
    assert row["schema"] == "isovar.sv_rna_orfs.v4"
    assert set(row["uncertainty_flags"].split(";")) == flags
    assert int(row["fragments"]) == sum(counts)
    assert tuple(int(row[k]) for k in ("original_query_fragments", "reverse_complement_fragments",
                                      "mixed_orientation_fragments")) == counts
