"""Same-sequence splice-linkage contrasts through public reconstruction/export."""

from copy import deepcopy
import csv
import json

import pysam
import pytest

from isovar import FusionBreakpoint, export_sv_rna_orfs, reconstruct_sv_rna, write_sv_rna_orfs
from isovar.rna_evidence import CELL_UMI_FIELDS
from tests.test_orf_start import reference
from tests.test_sv_rna import aligned, write_bam
from tests.testing_helpers import LINEAGE_SUMMARY


def inclusion_result(tmp_path, mode, assemble, **options):
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
            min_orf_amino_acids=1, assemble=assemble, **options)
    return result, sequence[3:-3]


@pytest.mark.parametrize("assemble", [False, True])
@pytest.mark.parametrize("mode,priority", [
    ("linked", 3), ("unspliced", 4), ("wrong_strand", 4),
    ("low_quality", 4), ("missing_quality", 4), ("unknown_mapq", 3),
    ("low_mapq", 4), ("deletion", 4),
])
def test_only_qualified_same_read_splice_linkage_promotes_intronic_atg(tmp_path, mode, priority, assemble):
    result, sequence = inclusion_result(tmp_path, mode, assemble)
    export = export_sv_rna_orfs(result)
    candidate, = [c for c in export["candidates"] if c["nucleotide_sequence"] == sequence]
    assert candidate["rna_support"]["fragments"] == 2
    assert candidate["start_evidence_summary"]["priority"] == priority
    assessment, = candidate["occurrences"][0]["start_evidence"]["assessments"]
    assert not candidate["initiation_observed"] and not candidate["translation_observed"]
    if mode == "wrong_strand":
        # An antisense start is not intronic, so no inclusion is assessed.
        assert assessment["region"] == "antisense" and "splice_inclusion" not in assessment
    else:
        evidence = assessment["splice_inclusion"]
        assert evidence["qualified_fragments"] == (2 if priority == 3 else 0)
        # Counts only: no read identities leave the reconstruction.
        assert set(evidence["support"]) == {"reads", "fragments", *CELL_UMI_FIELDS}
        assert set(evidence["read_lineage"]) == LINEAGE_SUMMARY
    if priority == 3:
        assert evidence["mechanisms"] == ["cryptic_donor"]
        assert all(w["minimum_base_quality"] == 30 for w in evidence["witnesses"])
        assert all(w["alignment_gap_kind"] == "N" for w in evidence["witnesses"])
    paths = write_sv_rna_orfs(export, tmp_path / "output")
    with paths["tsv"].open() as handle:
        row, = [r for r in csv.DictReader(handle, delimiter="\t") if r["candidate_id"] == candidate["candidate_id"]]
    assert int(row["start_priority"]) == priority
    assert json.loads(row["start_evidence_json"])[0]["start_evidence"]["priority"] == priority


@pytest.mark.parametrize("assemble", [False, True])
def test_mapq_255_is_unique_by_default_and_unavailable_on_request(tmp_path, assemble):
    result, sequence = inclusion_result(tmp_path, "unknown_mapq", assemble)
    assert "mapq_255_excluded_from_inclusion" not in result["limitations"]
    assert result["parameters"]["inclusion"] == dict(
        min_splice_anchor_bases=8, min_base_quality=20, min_mapping_quality=20, mapq_255_is_unique=True)
    candidate, = [c for c in export_sv_rna_orfs(result)["candidates"] if c["nucleotide_sequence"] == sequence]
    assert candidate["start_evidence_summary"]["priority"] == 3

    (tmp_path / "unavailable").mkdir()
    result, sequence = inclusion_result(tmp_path / "unavailable", "unknown_mapq", assemble,
                                        inclusion_mapq_255_is_unique=False)
    assert "mapq_255_excluded_from_inclusion" in result["limitations"]
    assert not result["parameters"]["inclusion"]["mapq_255_is_unique"]
    candidate, = [c for c in export_sv_rna_orfs(result)["candidates"] if c["nucleotide_sequence"] == sequence]
    evidence = candidate["occurrences"][0]["start_evidence"]["assessments"][0]["splice_inclusion"]
    assert candidate["start_evidence_summary"]["priority"] == 4
    assert {w["status"] for w in evidence["witnesses"]} == {"mapping_quality_unavailable"}
    assert all(w["mapping_quality_255"] and w["minimum_mapping_quality"] is None for w in evidence["witnesses"])


def test_inclusion_mapping_quality_threshold_is_configurable(tmp_path):
    result, sequence = inclusion_result(tmp_path, "linked", False, inclusion_min_mapping_quality=61)
    candidate, = [c for c in export_sv_rna_orfs(result)["candidates"] if c["nucleotide_sequence"] == sequence]
    evidence = candidate["occurrences"][0]["start_evidence"]["assessments"][0]["splice_inclusion"]
    assert candidate["start_evidence_summary"]["priority"] == 4
    assert evidence["thresholds"]["min_mapping_quality"] == 61
    assert {w["status"] for w in evidence["witnesses"]} == {"low_mapping_quality"}


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
@pytest.mark.parametrize("donor_end", [224, 324])
def test_rearranged_intronic_segment_needs_linked_partner_splicing(context, expected, donor_end):
    from isovar import FusionReference, annotate_orf_inclusion, annotate_orf_start
    from isovar.sv_rna import RnaObservation
    donor = reference()
    partner = FusionReference("partner", "synthetic", "test", "2", "+",
                              ((500, 524), (700, 724)), "C" * 48)
    positions = [("1", p, "+") for p in range(200, donor_end)]
    edge = len(positions)
    positions += [("2", p, "+") for p in range(500, 524)] + [("2", p, "+") for p in range(700, 724)]
    sequence = "CCC" + "ATG" + "C" * (len(positions) - 6)
    obs = RnaObservation("read", ("rg", "read", 0), (), sequence, tuple(positions),
                         ((edge - 1, edge, "split"), (edge + 23, edge + 24, "N" if context else "D")), False,
                         (0, len(sequence)), False, False, (30,) * len(sequence), 60)
    models = [donor, partner]
    evidence = annotate_orf_inclusion(annotate_orf_start(sequence, positions, 3, models),
                                     sequence, positions, len(sequence), models, {"read": obs},
                                     [dict(observation="read", observation_interval=[3, len(sequence)])])
    # A rearrangement after another exon does not establish inclusion of the
    # upstream retained intron, even in the same contiguous aligned run.
    expected = expected if donor_end == 224 else 4
    assert evidence["priority"] == expected
    if expected == 3:
        assert evidence["assessments"][0]["splice_inclusion"]["mechanisms"] == ["fusion_rearranged_exon"]


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("start_position,priority", [(180, 3), (230, 4)])
def test_cryptic_splice_must_include_the_start_intron(strand, start_position, priority):
    from isovar import FusionReference, annotate_orf_inclusion, annotate_orf_start
    from isovar.sv_rna import RnaObservation
    ref = FusionReference("introns", "synthetic", "test", "1", strand,
                          ((100, 120), (200, 220), (300, 320)), "C" * 60)
    coordinates = list(range(100, 120)) + list(range(160, 320))
    if strand == "-":
        coordinates.reverse()
    positions = [("1", p, strand) for p in coordinates]
    start = coordinates.index(start_position)
    sequence = "C" * start + "ATG" + "C" * (len(positions) - start - 3)
    edge = 20 if strand == "+" else 160
    obs = RnaObservation("read", ("rg", "read", 0), (), sequence, tuple(positions),
                         ((edge - 1, edge, "N"),), False, (0, len(sequence)), False, False,
                         (30,) * len(sequence), 60)
    result = annotate_orf_inclusion(annotate_orf_start(sequence, positions, start, [ref]),
                                    sequence, positions, len(sequence), [ref], {"read": obs},
                                    [dict(observation="read", observation_interval=[start, len(sequence)])])
    # The same contiguous alignment also spans a different retained intron.
    # Processing the first intron cannot establish inclusion of that second one.
    assert result["priority"] == priority


@pytest.mark.parametrize("options,expected", [
    ([], dict(inclusion_min_splice_anchor_bases=8, inclusion_min_base_quality=20,
              inclusion_min_mapping_quality=20, inclusion_mapq_255_is_unique=True)),
    (["--inclusion-min-splice-anchor-bases", "10", "--inclusion-min-base-quality", "25",
      "--inclusion-min-mapping-quality", "30", "--inclusion-mapq-255-unavailable"],
     dict(inclusion_min_splice_anchor_bases=10, inclusion_min_base_quality=25,
          inclusion_min_mapping_quality=30, inclusion_mapq_255_is_unique=False)),
])
def test_cli_passes_inclusion_gates_to_reconstruction(tmp_path, monkeypatch, options, expected):
    from contextlib import nullcontext
    from isovar.cli import isovar_sv_rna
    seen = {}
    event = tmp_path / "input.json"
    event.write_text("{}")
    monkeypatch.setattr(isovar_sv_rna, "sv_rna_input_from_dict", lambda value: {})
    monkeypatch.setattr(isovar_sv_rna, "alignment_file_from_args", lambda options: nullcontext(None))
    monkeypatch.setattr(isovar_sv_rna, "reconstruct_sv_rna", lambda *args, **kwargs: seen.update(kwargs) or {})
    isovar_sv_rna.run(["--input", str(event), "--bam", "rna.bam", "--output", str(tmp_path / "out.json")] + options)
    assert {key: seen[key] for key in expected} == expected
