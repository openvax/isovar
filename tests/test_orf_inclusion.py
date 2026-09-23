"""Same-sequence splice-linkage contrasts through public reconstruction/export."""

from copy import deepcopy
import csv
import json

import pysam
import pytest

from isovar import FusionBreakpoint, export_sv_rna_orfs, reconstruct_sv_rna, write_sv_rna_orfs
from tests.test_orf_start import reference
from tests.test_sv_rna import aligned, write_bam


def inclusion_result(tmp_path, mode, assemble):
    ref = reference()
    intronic = "CCCATG" + "GCC" * 6
    exon = ref.sequence[ref.cds_start:ref.cds_start + 24]
    sequence = intronic + exon + "GAAGAATAACCC"
    start = 200
    positions = [("1", start + q, "+") for q in range(len(intronic))]
    exon_start = start + len(intronic) if mode == "unspliced" else 300
    positions += [("1", exon_start + q, "+") for q in range(len(exon))]
    positions += [("2", 1000 + q, "+") for q in range(len(sequence) - len(intronic) - len(exon))]
    if mode == "wrong_strand":
        from dataclasses import replace
        ref = replace(ref, strand="-")
    reads = [r for i in range(2) for r in aligned("read%d" % i, sequence, positions,
                                               qualities=mode != "missing_quality")]
    for read in reads:
        if mode == "low_quality":
            read.query_qualities = [5] * len(read.query_sequence)
        if mode == "unknown_mapq":
            read.mapping_quality = 255
        if mode == "low_mapq":
            read.mapping_quality = 1
        if mode == "deletion":
            read.cigarstring = read.cigarstring.replace("N", "D")
    # Update SA after editing CIGAR so records remain compatible paths.
    if mode == "deletion":
        from tests.test_sv_rna import link
        for i in range(2):
            link([r for r in reads if r.query_name == "read%d" % i])
    path = write_bam(tmp_path / (mode + ".bam"), reads)
    from isovar import ReadCollector
    with pysam.AlignmentFile(path) as bam:
        result = reconstruct_sv_rna(
            bam, event_id="splice-linked-ATG", reference_name="synthetic", references=[ref],
            donor=FusionBreakpoint("1", exon_start + len(exon), "+"),
            acceptor=FusionBreakpoint("2", 1000, "+"), regions=(),
            sample_id="sample", source="synthetic-library", event_provenance={"fixture": "inclusion_result"},
            read_collector=ReadCollector(min_mapping_quality=0),
            min_orf_amino_acids=1, assemble=assemble)
    return result, sequence[3:-3]


@pytest.mark.parametrize("assemble", [False, True])
@pytest.mark.parametrize("mode,priority", [
    ("linked", 3), ("unspliced", 4), ("wrong_strand", 4),
    ("low_quality", 4), ("missing_quality", 4), ("unknown_mapq", 4),
    ("low_mapq", 4), ("deletion", 4),
])
def test_only_qualified_same_read_splice_linkage_promotes_intronic_atg(tmp_path, mode, priority, assemble):
    result, sequence = inclusion_result(tmp_path, mode, assemble)
    export = export_sv_rna_orfs(result)
    candidate, = [c for c in export["candidates"] if c["nucleotide_sequence"] == sequence]
    assert candidate["rna_support"]["fragments"] == 2
    assert candidate["start_evidence_summary"]["priority"] == priority
    evidence = candidate["occurrences"][0]["start_evidence"]["assessments"][0]["splice_inclusion"]
    assert evidence["qualified_fragments"] == (2 if priority == 3 else 0)
    assert "segment_ids" not in evidence["cell_umi_support"]
    assert "segment_ids" not in evidence["read_lineage"]
    assert not candidate["initiation_observed"] and not candidate["translation_observed"]
    if priority == 3:
        assert evidence["mechanisms"] == ["cryptic_donor"]
        assert all(w["minimum_base_quality"] == 30 for w in evidence["witnesses"])
        assert all(w["alignment_gap_kind"] == "N" for w in evidence["witnesses"])
    paths = write_sv_rna_orfs(export, tmp_path / "output")
    with paths["tsv"].open() as handle:
        row, = [r for r in csv.DictReader(handle, delimiter="\t") if r["candidate_id"] == candidate["candidate_id"]]
    assert int(row["start_priority"]) == priority
    assert json.loads(row["start_evidence_json"])[0]["start_evidence"]["priority"] == priority


def test_unrelated_splice_or_partial_witness_cannot_promote_start(tmp_path):
    from isovar import annotate_orf_inclusion, annotate_orf_start
    from isovar.sv_rna import RnaObservation
    ref = reference()
    sequence = "ATG" + "GCC" * 20
    positions = [("1", q, "+") for q in range(180, 204)] + [("1", q, "+") for q in range(220, 240)]
    positions += [("1", q, "+") for q in range(300, 300 + len(sequence) - len(positions))]
    obs = RnaObservation("read", ("rg", "read", 0), (), sequence, tuple(positions),
                         ((23, 24, "N"), (43, 44, "N")), False, (0, len(sequence)), False, False,
                         (30,) * len(sequence), 60)
    annotation = annotate_orf_start(sequence, positions, 0, [ref])
    before = deepcopy(annotation)
    evidence = annotate_orf_inclusion(annotation, sequence, positions, len(sequence), [ref], {"read": obs},
                                     [dict(observation="read", observation_interval=[0, len(sequence)])])
    assert evidence["priority"] == 4  # Exon-linked splice is on another contiguous block.
    assert annotation == before
    evidence = annotate_orf_inclusion(annotation, sequence, positions, len(sequence), [ref], {"read": obs}, [])
    assert evidence["priority"] == 4


@pytest.mark.parametrize("competitor,context,boundaries,expected", [
    (True, True, True, 3), (False, True, True, 4),
    (True, False, True, 4), (True, True, False, 4),
])
def test_retention_requires_both_boundaries_processed_context_and_competing_splice(
        competitor, context, boundaries, expected):
    from isovar import FusionReference, annotate_orf_inclusion, annotate_orf_start
    from isovar.sv_rna import RnaObservation
    ref = FusionReference("retained", "synthetic", "test", "1", "+",
                          ((100, 120), (200, 220), (300, 320)), "C" * 60)
    first = 100 if boundaries else 125
    positions = [("1", p, "+") for p in range(first, 220)] + [("1", p, "+") for p in range(300, 320)]
    start = 130 - first
    sequence = "C" * start + "ATG" + "C" * (len(positions) - start - 3)
    edge = 220 - first
    obs = RnaObservation("read", ("rg", "read", 0), (), sequence, tuple(positions),
                         ((edge - 1, edge, "N" if context else "D"),), False,
                         (0, len(sequence)), False, False, (30,) * len(sequence), 60)
    annotation = annotate_orf_start(sequence, positions, start, [ref])
    competitors = [dict(left=["1", 119, "+"], right=["1", 200, "+"], fragments=2)] if competitor else []
    evidence = annotate_orf_inclusion(annotation, sequence, positions, len(sequence), [ref], {"read": obs},
                                     [dict(observation="read", observation_interval=[start, len(sequence)])],
                                     competing_splices=competitors)
    assert evidence["priority"] == expected
    if expected == 3:
        assert evidence["assessments"][0]["splice_inclusion"]["mechanisms"] == ["intron_retention"]


@pytest.mark.parametrize("context,expected", [(True, 3), (False, 4)])
def test_rearranged_intronic_segment_needs_linked_partner_splicing(context, expected):
    from isovar import FusionReference, annotate_orf_inclusion, annotate_orf_start
    from isovar.sv_rna import RnaObservation
    donor = reference()
    partner = FusionReference("partner", "synthetic", "test", "2", "+",
                              ((500, 524), (700, 724)), "C" * 48)
    positions = [("1", p, "+") for p in range(200, 224)]
    positions += [("2", p, "+") for p in range(500, 524)] + [("2", p, "+") for p in range(700, 724)]
    sequence = "CCC" + "ATG" + "C" * (len(positions) - 6)
    obs = RnaObservation("read", ("rg", "read", 0), (), sequence, tuple(positions),
                         ((23, 24, "split"), (47, 48, "N" if context else "D")), False,
                         (0, len(sequence)), False, False, (30,) * len(sequence), 60)
    models = [donor, partner]
    evidence = annotate_orf_inclusion(annotate_orf_start(sequence, positions, 3, models),
                                     sequence, positions, len(sequence), models, {"read": obs},
                                     [dict(observation="read", observation_interval=[3, len(sequence)])])
    assert evidence["priority"] == expected
    if expected == 3:
        assert evidence["assessments"][0]["splice_inclusion"]["mechanisms"] == ["fusion_rearranged_exon"]
