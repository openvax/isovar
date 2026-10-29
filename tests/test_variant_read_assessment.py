"""RNA survives missing annotation without weakening coding/path safeguards."""

import json
import pickle

import pytest
from varcode import Variant

from isovar import AlleleRead, IsovarResult, ProteinSequenceCreator, ReadEvidence, export_protein_hypotheses
from isovar.main import _assess_and_select
from isovar.rna_evidence import EvidenceSets
from isovar.variant_read_assessment import assemble_variant_reads

from .test_translation_regressions import _context, _read
from .test_transcript_compatibility import SHORT, LONG


def no_cds(monkeypatch):
    monkeypatch.setattr("isovar.protein_sequence_creator.reference_contexts_for_variant", lambda *a, **k: [])
    return Variant("1", 100, "G", "C", "GRCh38")


@pytest.mark.parametrize("mode", ["configured", "read_contexts"])
def test_no_coding_annotation_retains_rna_and_original_support(monkeypatch, mode):
    variant = no_cds(monkeypatch)
    reads = [_read("ACG" * 10, "C", "GTA" * 10, f"read{i}", "+") for i in range(3)]
    creator = ProteinSequenceCreator()
    assessment = creator.assess_variant_reads(variant, iter(reads), assembly_mode=mode)
    assert assessment.assembly.state == "assembled"
    assert assessment.frame_state == "no_coding_reference_context"
    assert not assessment.translations and not creator.protein_sequences_from_assessment(assessment)
    assert creator.translate_variant_reads(variant, reads) == []  # list API remains compatible
    exported = assessment.to_dict(EvidenceSets(["sample", "library"]))
    assert exported["assembly_mode"] == mode
    assert exported["reference_incompatible_support"]["fragments"] == 3
    assert exported["rna_candidates"]
    for row in exported["rna_candidates"]:
        left, right = row["variant_interval"]
        assert row["sequence"][left:right] == "C"
        assert row["frame_state"] == "no_coding_reference_context"
        assert row["rna_support"]["fragments"] == 3
    assert json.loads(json.dumps(exported)) == exported


def test_absent_reads_coverage_failure_and_wrong_allele_are_distinct(monkeypatch):
    variant = no_cds(monkeypatch)
    creator = ProteinSequenceCreator(min_variant_sequence_coverage=2)
    empty = creator.assess_variant_reads(variant, [])
    assert empty.frame_state == empty.assembly.state == "no_alt_reads"
    one = creator.assess_variant_reads(variant, [_read("ACG" * 10, "C", "GTA" * 10, "one", "+")])
    assert one.assembly.state == "no_retained_rna"
    assert one.frame_state == "no_coding_reference_context"
    with pytest.raises(ValueError, match="nominated ALT"):
        creator.assess_variant_reads(variant, [_read("ACG", "G", "GTA", "reference", "+")])
    with pytest.raises(ValueError, match="assembly mode"):
        creator.assess_variant_reads(variant, [], assembly_mode="unknown")


@pytest.mark.parametrize("mode", ["configured", "read_contexts"])
def test_incompatible_retention_is_kept_separately_from_spliced_path(mode):
    variant = Variant("1", 119, "G", "C", "GRCh38")
    reads = [AlleleRead("AAAA", "C", "TTTT", name=name, reference_blocks=blocks,
                        splice_junctions=junctions)
             for name, blocks, junctions in [
                 ("retained", ((0, 9, 114, 123),), ()),
                 ("spliced", ((0, 6, 114, 120), (6, 9, 200, 203)), ((120, 200),))]]
    creator = ProteinSequenceCreator(min_variant_sequence_coverage=1)
    assembly = assemble_variant_reads(variant, reads, [SHORT], creator=creator, mode=mode)
    assert {r.name for r in assembly.compatible_reads} == {"spliced"}
    assert {r.name for r in assembly.unmatched_reads} == {"retained"}
    assert assembly.compatible_sequences and assembly.unmatched_sequences
    assert all(s.compatible_transcript_ids == {"short"} for s in assembly.compatible_sequences)
    assert all(s.compatible_transcript_ids is None for s in assembly.unmatched_sequences)
    # Adding a credible retained-extension model makes the same observation compatible.
    competing = assemble_variant_reads(variant, reads, [SHORT, LONG], creator=creator, mode=mode)
    assert not competing.unmatched_reads
    assert {r.name: r.compatible_transcript_ids for r in competing.compatible_reads} == {
        "retained": {"long"}, "spliced": {"short"}}


@pytest.mark.parametrize("mode", ["configured", "read_contexts"])
def test_alternative_placements_cannot_satisfy_coverage_floor(monkeypatch, mode):
    variant = no_cds(monkeypatch)
    reads = [AlleleRead("ACG" * 10, "C", "GTA" * 10, name="same", source_alignments=(
        (("rg", "same", 0), (0, position, "61M", False)),)) for position in (69, 70)]
    creator = ProteinSequenceCreator(min_variant_sequence_coverage=2)
    assessment = creator.assess_variant_reads(variant, reads, assembly_mode=mode)
    assert assessment.assembly.state == "no_retained_rna"
    assert assessment.to_dict(EvidenceSets(["S", "R"]))["reference_incompatible_support"]["fragments"] == 1
    creator = ProteinSequenceCreator(min_variant_sequence_coverage=1, protein_sequence_preference="support")
    rows = creator.assess_variant_reads(variant, reads, assembly_mode=mode).to_dict(
        EvidenceSets(["S", "R"]))["rna_candidates"]
    assert len(rows) == len({row["rna_candidate_id"] for row in rows}) == 2


def test_reference_matches_reuse_one_assembly_and_keep_per_path_failures(monkeypatch):
    variant = Variant("1", 100, "G", "C", "GRCh38")
    context = _context(variant, "+", "A" * 30, "G" * 30)
    monkeypatch.setattr("isovar.protein_sequence_creator.reference_contexts_for_variant", lambda *a, **k: [context])
    reads = [_read(prefix, "C", "G" * 30, f"{prefix[0]}{i}", "+")
             for prefix in ("A" * 30, "T" * 30) for i in range(2)]

    class CountedCreator(ProteinSequenceCreator):
        calls = 0

        def assess_variant_reads(self, *args, **kwargs):
            self.calls += 1
            return super().assess_variant_reads(*args, **kwargs)

    creator = CountedCreator(max_transcript_mismatches=0)
    evidence = ReadEvidence.from_variant_and_allele_reads(variant, reads)
    assessment, proteins = _assess_and_select(creator, variant, evidence, None)
    assert creator.calls == 1 and proteins
    exported = assessment.to_dict(EvidenceSets(["sample", "library"]))
    assert exported["frame_state"] == "rna_frame_candidates"
    assert {row["frame_state"] for row in exported["rna_candidates"]} == {
        "rna_frame_candidates", "no_retained_rna_frame_match"}
    # This supplied ReferenceContext has a frame but no named transcript path.
    assert all(row["frame_transcript_ids"] == [] for row in exported["rna_candidates"])
    # Selection never consumes or erases the assessment.
    assert creator.protein_sequences_from_assessment(assessment) == proteins


def test_legacy_creator_hooks_remain_authoritative(monkeypatch):
    variant = no_cds(monkeypatch)
    evidence = ReadEvidence.from_variant_and_allele_reads(variant, [])
    sentinel = object()

    class CustomCreator(ProteinSequenceCreator):
        def sorted_protein_sequences_for_variant(self, **kwargs):
            return [sentinel]

    assessment, proteins = _assess_and_select(CustomCreator(), variant, evidence, None)
    assert assessment is None and proteins == [sentinel]


def test_result_clone_pickle_and_export_keep_unframed_rna(monkeypatch):
    variant = no_cds(monkeypatch)
    reads = [_read("ACG" * 10, "C", "GTA" * 10, f"read{i}", "+") for i in range(2)]
    reads = [AlleleRead(read.prefix, read.allele, read.suffix, name=read.name,
                         source_alignments=((("rg", read.name, 0), (0, 69, "61M", False)),)) for read in reads]
    evidence = ReadEvidence.from_variant_and_allele_reads(variant, reads)
    assessment = ProteinSequenceCreator().assess_variant_reads(variant, reads)
    result = IsovarResult(variant, evidence, None, rna_assessment=assessment)
    cloned = result.clone_with_updates(protein_sequence_settings=ProteinSequenceCreator().settings())
    assert cloned.rna_assessment is assessment
    restored = pickle.loads(pickle.dumps(cloned))
    export = export_protein_hypotheses([restored], sample_id="sample", source="library")
    row, = export["events"]
    assert not row["protein_hypotheses"]
    assert row["rna_assessment"]["frame_state"] == "no_coding_reference_context"
    assert row["rna_assessment"]["rna_candidates"]
    evidence_ids = set(export["evidence_sets"])
    assert all(r["rna_support"]["evidence_set_id"] in evidence_ids for r in row["rna_assessment"]["rna_candidates"])
    del restored.rna_assessment  # pre-assessment pickle semantics
    row, = export_protein_hypotheses([restored], sample_id="sample", source="library")["events"]
    assert row["rna_assessment"] is None
