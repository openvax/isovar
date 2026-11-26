"""Focused unit and regression tests for Varcode 13 adoption (#741)."""

import json
from types import SimpleNamespace

from pyensembl import cached_release
import varcode
from varcode import (
    EffectCandidate,
    GermlineContext,
    MutantTranscript,
    StructuralVariant,
    Variant,
    parse_symbolic_alt,
)
from varcode.effect_hypotheses import (
    ClassifiedOutcome,
    EffectHypothesis,
    RealizedEffectCandidate,
    RealizedTranscriptProduct,
)
from varcode.effects import MutationEffect

from isovar.effect_prediction import effect_protein_comparisons
from isovar.sv_rna_reconciliation import (
    hypotheses_from_varcode,
    reconcile_sv_rna_hypotheses,
)


def _synthetic_effect(model, *, protein="MK", transcript_id="ENST_SYN", desc="synthetic"):
    return SimpleNamespace(
        mutant_transcript=model,
        mutant_protein_sequence=protein,
        transcript=SimpleNamespace(id=transcript_id),
        short_description=desc,
    )


def test_hypotheses_from_varcode_preserves_candidate_and_mutant_transcript_metadata():
    model = MutantTranscript(
        reference_segments=(),
        cdna_sequence="ATGAAATAA",
        mutant_protein_sequence="MK",
        annotator_name="dna",
        evidence={
            "novel_junctions": (("7", 116771655, 116774880, "+"),),
            "back_junctions": (("7", 55173922, 55201354, "+"),),
            "wildtype_junctions": (("7", 116771655, 116771848, "+"),),
            "mate_search_region": ("6", 44218908, 44243907, "+"),
            "exit_frame": 2,
        },
    )
    inner = _synthetic_effect(model)
    germline_like = SimpleNamespace(__str__=lambda self: "chr7:140753335C>T")
    cis_cand = EffectCandidate(
        effect=inner,
        source="varcode.dna",
        evidence={"n_merged_outcomes": 1},
        mechanisms=("exon_skip",),
        phase=((germline_like, "cis"),),
        uncertainty=("phasing", "splicing"),
    )
    trans_cand = EffectCandidate(
        effect=inner,
        source="varcode.dna",
        evidence={"n_merged_outcomes": 1},
        mechanisms=("exon_skip",),
        phase=((germline_like, "trans"),),
        uncertainty=("phasing", "splicing"),
    )

    rows = hypotheses_from_varcode([SimpleNamespace(candidates=(cis_cand, trans_cand))])
    assert len(rows) == 2
    assert len({row["hypothesis_id"] for row in rows}) == 2
    assert {row["phase"] for row in rows} == {
        ((str(germline_like), "cis"),),
        ((str(germline_like), "trans"),),
    }
    for row in rows:
        assert row["sequence"] == "ATGAAATAA"
        assert row["mechanisms"] == ("exon_skip",)
        assert row["uncertainty"] == ("phasing", "splicing")
        assert row["is_complete"] is True
        assert row["novel_junctions"] == (("7", 116771655, 116774880, "+"),)
        assert row["back_junctions"] == (("7", 55173922, 55201354, "+"),)
        assert row["wildtype_junctions"] == (("7", 116771655, 116771848, "+"),)
        assert row["mate_region"] == ("6", 44218908, 44243907, "+")
        assert row["exit_frame"] == 2
        assert row["source"] == {
            "producer": "varcode.dna",
            "varcode_version": varcode.__version__,
        }
        # Strict JSON serialization used by sv_rna_audit must succeed.
        json.dumps(row, sort_keys=True, allow_nan=False)


def test_hypotheses_from_varcode_open_template_materialization_unit():
    open_model = MutantTranscript(
        reference_segments=(),
        cdna_sequence=None,
        mutant_protein_sequence=None,
        annotator_name="dna",
        evidence={
            "sequence_status": "open_search_template",
            "wildtype_junctions": (("7", 101, 200, "+"),),
            "mate_search_region": ("6", 500, 505, "-"),
            "exit_frame": 1,
            "inserted_sequence": "GGGG",
            "template_segments": (
                ("five_prime_host", "7", 95, 100, "+"),
                ("inserted_sequence", "7", 100, 100, "+"),
                ("mate_search_window", "6", 500, 505, "-"),
            ),
        },
    )
    incomplete_without_segments = MutantTranscript(
        reference_segments=(),
        cdna_sequence=None,
        mutant_protein_sequence=None,
        annotator_name="dna",
        evidence={"sequence_status": "unresolved_splice"},
    )
    complete_model = MutantTranscript(
        reference_segments=(),
        cdna_sequence="ATGCCC",
        mutant_protein_sequence="MP",
        annotator_name="dna",
        evidence={
            "template_segments": (("five_prime_host", "7", 95, 100, "+"),),
        },
    )
    effect = SimpleNamespace(
        candidates=(
            EffectCandidate(
                effect=_synthetic_effect(open_model, protein=None, desc="open"),
                source="varcode.dna",
                uncertainty=("structural",),
            ),
            EffectCandidate(
                effect=_synthetic_effect(incomplete_without_segments, protein=None, desc="nosig"),
                source="varcode.dna",
                mechanisms=("cryptic",),
                uncertainty=("splicing",),
            ),
            EffectCandidate(
                effect=_synthetic_effect(complete_model, protein="MP", desc="complete"),
                source="varcode.dna",
                mechanisms=("normal",),
            ),
        )
    )

    without_provider = {r["description"]: r for r in hypotheses_from_varcode([effect])}
    assert without_provider["open"]["sequence"] is None
    assert without_provider["open"]["is_complete"] is False
    assert without_provider["open"]["mate_region"] == ("6", 500, 505, "-")
    assert without_provider["open"]["wildtype_junctions"] == (("7", 101, 200, "+"),)
    assert without_provider["open"]["exit_frame"] == 1
    assert without_provider["nosig"]["sequence"] is None
    assert without_provider["nosig"]["is_complete"] is False
    assert without_provider["complete"]["sequence"] == "ATGCCC"
    assert without_provider["complete"]["is_complete"] is True

    calls = []

    def provider(contig, start, end):
        calls.append((contig, start, end))
        if contig == "7":
            return "acgtac"
        if contig == "6":
            return "tttaaa"
        raise AssertionError(contig)

    with_provider = {
        r["description"]: r
        for r in hypotheses_from_varcode([effect], sequence_provider=provider)
    }
    # Complete model keeps its spliced cDNA; open model materializes 5'->3'
    # pre-mRNA template ("ACGTAC" + "GGGG" + RC("TTTAAA") == "ACGTACGGGGTTTAAA").
    assert calls == [("7", 95, 100), ("6", 500, 505)]
    assert with_provider["open"]["sequence"] == "ACGTACGGGGTTTAAA"
    assert with_provider["open"]["is_complete"] is False
    assert with_provider["nosig"]["sequence"] is None
    assert with_provider["complete"]["sequence"] == "ATGCCC"
    assert with_provider["complete"]["is_complete"] is True

    # Materialized open template is directly reconcilable with observed RNA windows.
    window = "ACGTACGGGGTTTAAA"
    assessment = dict(
        event_id="open-event",
        sample_id="sample",
        source="library",
        reference_name="test",
        reference_models=[],
        limitations=[],
        parameters=dict(junction_flank=4),
        paths=[
            dict(
                path_id="observed",
                sequence=window,
                junctions=[
                    dict(
                        window_sequence=window,
                        relation="breakpoint_junction",
                        exact_sequence_support=dict(complete_fragments=[["rg", "f1"]]),
                    )
                ],
            )
        ],
    )
    reconciled = reconcile_sv_rna_hypotheses(assessment, [with_provider["open"]])
    open_row = next(
        row for row in reconciled["hypotheses"]
        if row["hypothesis"]["hypothesis_id"] == with_provider["open"]["hypothesis_id"]
    )
    assert open_row["sequence_state"] == "materialized"
    assert open_row["status"] == "indistinguishable_in_observed_RNA"


def test_real_varcode13_splice_kdd_phase_and_open_templates():
    genome = cached_release(87)

    # 1. + strand exon-skipping splice junction (MET exon 14 donor +1)
    met = genome.transcript_by_id("ENST00000397752")
    ex13, ex14, ex15 = met.exons[12], met.exons[13], met.exons[14]
    met_eff = Variant("7", ex14.end + 1, "G", "A", genome=genome).effect_on_transcript(met)
    met_rows = hypotheses_from_varcode([met_eff])
    skip_met = next(r for r in met_rows if "exon_skip" in r["mechanisms"])
    assert skip_met["uncertainty"] == ("splicing",)
    assert skip_met["is_complete"] is True
    assert skip_met["novel_junctions"] == (("7", ex13.end + 1, ex15.start - 1, "+"),)
    assert skip_met["novel_junctions"][0][1] <= skip_met["novel_junctions"][0][2]

    # 2. - strand exon-skipping splice junction (TP53 exon 5 donor +1): ascending coords
    tp53 = genome.transcript_by_id("ENST00000269305")
    ex4, ex5, ex6 = tp53.exons[3], tp53.exons[4], tp53.exons[5]
    tp53_eff = Variant("17", ex5.start - 1, "C", "T", genome=genome).effect_on_transcript(tp53)
    tp53_rows = hypotheses_from_varcode([tp53_eff])
    skip_tp53 = next(r for r in tp53_rows if "exon_skip" in r["mechanisms"])
    assert skip_tp53["novel_junctions"] == (("17", ex6.end + 1, ex4.start - 1, "-"),)
    assert skip_tp53["novel_junctions"][0][1] <= skip_tp53["novel_junctions"][0][2]

    # 3. EGFR-KDD tandem duplication back-junction (exon 25 -> exon 18)
    egfr = genome.transcript_by_id("ENST00000275493")
    ex18, ex25 = egfr.exons[17], egfr.exons[24]
    kdd_eff = StructuralVariant(
        "7", ex18.start - 500, "DUP", end=ex25.end + 150, genome=genome
    ).effect_on_transcript(egfr)
    (kdd_row,) = hypotheses_from_varcode([kdd_eff])
    assert kdd_row["is_complete"] is True
    assert kdd_row["novel_junctions"] == ()
    assert kdd_row["back_junctions"] == (("7", ex18.start + 1, ex25.end - 1, "+"),)
    assert kdd_row["wildtype_junctions"]
    assert kdd_row["exit_frame"] == 0

    # 4. Unphased germline phase preservation (BRAF V600E + adjacent germline SNV)
    braf = genome.transcript_by_id("ENST00000288602")
    v_braf = Variant("7", 140753336, "A", "T", genome=genome)
    g_braf = Variant("7", 140753335, "C", "T", genome=genome)
    braf_eff = v_braf.effect_on_transcript(braf, germline=GermlineContext.from_variants([g_braf]))
    braf_rows = hypotheses_from_varcode([braf_eff])
    assert len(braf_rows) == 2
    assert all(r["uncertainty"] == ("phasing",) for r in braf_rows)
    assert {state for r in braf_rows for _, state in r["phase"]} == {"cis", "trans"}
    for row in braf_rows:
        json.dumps(row, sort_keys=True, allow_nan=False)

    # 5. Open rearrangement template: large SV insertion & gene-to-intergenic translocation
    big_insert = "GCTAGC" * 80
    sv_ins = StructuralVariant(
        "7", ex14.end + 400, "INS", info={"SVINSSEQ": big_insert}, genome=genome
    )
    ins_eff = sv_ins.effect_on_transcript(met)
    (ins_unresolved,) = hypotheses_from_varcode([ins_eff])
    assert ins_unresolved["is_complete"] is False
    assert ins_unresolved["sequence"] is None
    assert ins_unresolved["exit_frame"] in (0, 1, 2)
    assert ("7", ex14.end + 1, ex15.start - 1, "+") in ins_unresolved["wildtype_junctions"]

    def stub_provider(contig, start, end):
        assert start <= end
        return "A" * (end - start + 1)

    (ins_materialized,) = hypotheses_from_varcode([ins_eff], sequence_provider=stub_provider)
    assert ins_materialized["is_complete"] is False
    assert ins_materialized["sequence"] is not None
    assert big_insert in ins_materialized["sequence"]

    sv_bnd = parse_symbolic_alt("17", ex5.start - 40, "N", "]10:40000000]N", genome=genome)
    bnd_eff = sv_bnd.effect_on_transcript(tp53)
    (bnd_unresolved,) = hypotheses_from_varcode([bnd_eff])
    assert bnd_unresolved["is_complete"] is False
    assert bnd_unresolved["sequence"] is None
    assert bnd_unresolved["mate_region"] == ("10", 40_000_000 - 25_000 + 1, 40_000_000, "-")
    assert bnd_unresolved["mate_region"][1] <= bnd_unresolved["mate_region"][2]
    assert bnd_unresolved["exit_frame"] in (0, 1, 2)

    def strand_provider(contig, start, end):
        return ("A" * (end - start)) + "C"

    (bnd_materialized,) = hypotheses_from_varcode([bnd_eff], sequence_provider=strand_provider)
    assert bnd_materialized["is_complete"] is False
    assert bnd_materialized["sequence"].startswith("G")


def test_effect_protein_comparisons_parity_between_dna_and_transcript_model():
    # Unit check: EffectCandidate with _paired_outcomes preserves distinct patient baselines.
    eff = MutationEffect(variant=None)
    eff._original_protein_sequence = "MAYDEFG"
    eff.mutant_protein_sequence = "MACDEFG"
    cis_outcome = ClassifiedOutcome(
        EffectHypothesis(evidence={"phase": "cis"}),
        eff,
        RealizedTranscriptProduct(None, None, "MAWDEFG"),
        RealizedTranscriptProduct(None, None, "MACDEFG"),
    )
    trans_outcome = ClassifiedOutcome(
        EffectHypothesis(evidence={"phase": "trans"}),
        eff,
        RealizedTranscriptProduct(None, None, "MAYDEFG"),
        RealizedTranscriptProduct(None, None, "MACDEFG"),
    )
    dna_cand = EffectCandidate(effect=eff, source="varcode.dna", uncertainty=("phasing",))
    object.__setattr__(dna_cand, "_paired_outcomes", (cis_outcome, trans_outcome))
    eff.candidates = (dna_cand,)
    assert list(effect_protein_comparisons(eff)) == [
        (eff, "MAWDEFG", "MACDEFG"),
        (eff, "MAYDEFG", "MACDEFG"),
    ]

    # Real comparison parity between default "dna" and explicit "transcript_model"
    # on a phased BRAF V600E + missense germline SNV in the same codon (V600M),
    # producing distinct patient baselines (V600 vs M600) and mutants (V600E vs M600K).
    genome = cached_release(87)
    braf = genome.transcript_by_id("ENST00000288602")
    v_braf = Variant("7", 140753336, "A", "T", genome=genome)
    g_braf = Variant("7", 140753337, "C", "T", genome=genome)
    ctx = GermlineContext.from_variants([g_braf])

    dna_eff = v_braf.effect_on_transcript(braf, germline=ctx)
    tm_eff = v_braf.effect_on_transcript(braf, germline=ctx, annotator="transcript_model")
    assert all(isinstance(c, EffectCandidate) for c in dna_eff.candidates)
    assert all(isinstance(c, RealizedEffectCandidate) for c in tm_eff.candidates)

    dna_rows = [(e.short_description, b, m) for e, b, m in effect_protein_comparisons(dna_eff)]
    tm_rows = [(e.short_description, b, m) for e, b, m in effect_protein_comparisons(tm_eff)]
    assert dna_rows == tm_rows
    assert len({b for _, b, _ in dna_rows}) == 2
    assert len({m for _, _, m in dna_rows}) == 2
