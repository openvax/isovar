"""Changed residues require a verified comparison, not a marked codon interval."""

import pytest

from isovar.protein_window_comparison import compare_protein_windows, reference_window_stop


@pytest.mark.parametrize("reference,observed,interval,delta,count,field,positions", [
    ("MAKG", "MAKG", (2, 3), 0, 0, "substituted_positions", []),
    ("MAKG", "MAQG", (2, 3), 0, 1, "substituted_positions", [2]),
    ("MAKKG", "MAKG", (2, 2), 3, 1, "deleted_reference_positions", [2]),
    ("MAKG", "MAKKG", (2, 3), -3, 1, "inserted_positions", [2]),
])
def test_identity_substitution_and_focal_gaps(reference, observed, interval, delta, count, field, positions):
    result = compare_protein_windows(reference, observed, interval, allele_length_delta=delta, reference_stop=False)
    assert result["status"] == ("amino_acid_divergent" if count else "reference_identical")
    assert result["changed_residue_count"] == count
    assert result[field] == positions


def test_a_partial_reference_does_not_invent_an_insertion():
    result = compare_protein_windows("MA", "MAKG", (2, 3), reference_stop=False)
    assert result["status"] == "unresolved" and result["changed_residue_count"] is None
    assert result["inserted_positions"] == []
    assert result["unresolved_alignment_columns"] == [2, 3]


def test_frameshift_tail_counts_actual_column_differences():
    result = compare_protein_windows("MAKG", "MAQG", (2, 4), allele_length_delta=1,
                                     frameshift=True, reference_stop=False)
    assert result["changed_residue_count"] == 1
    assert result["substituted_positions"] == [2]


def test_stop_change_is_separate_from_a_residue_count():
    result = compare_protein_windows("MAKG", "MAKG", (4, 4), stop=True, reference_stop=False)
    assert result["status"] == "amino_acid_divergent"
    assert result["changed_residue_count"] == 0 and result["stop_changed"]


@pytest.mark.parametrize("reference,observed", [(None, "MAKG"), ("MAXG", "MAXG"), ("MAKG", "MAXG")])
def test_missing_or_ambiguous_residues_do_not_establish_identity(reference, observed):
    result = compare_protein_windows(reference, observed, (2, 3), reference_stop=False)
    assert result["status"] == "unresolved" and result["changed_residue_count"] is None


def test_known_divergence_does_not_hide_unknown_residue_count():
    result = compare_protein_windows("MAK", "MQKG", (1, 2), reference_stop=False)
    assert result["status"] == "amino_acid_divergent" and result["changed_residue_count"] is None
    assert result["substituted_positions"] == [1]


def test_shifted_stop_is_not_reference_identity():
    result = compare_protein_windows("MA", "MAKG", (2, 4), frameshift=True, stop=True, reference_stop=True)
    assert result["status"] == "amino_acid_divergent"
    assert result["stop_changed"] and result["changed_residue_count"] is None


def test_missing_next_reference_codon_does_not_establish_a_stop_gain():
    assert reference_window_stop("MAKG", False, 4, True) is None
    assert reference_window_stop("MAKGQ", False, 4, True) is False
    assert reference_window_stop("MAKG", True, 4, True) is True
    result = compare_protein_windows("MAKG", "MAKG", (4, 4), stop=True, reference_stop=None)
    assert result["status"] == "unresolved" and result["stop_changed"] is None


@pytest.mark.parametrize("background,protein,annotation_status,count", [
    ("GCT", "MAKG", "reference_identical", 0),
    ("GAT", "MDKG", "amino_acid_divergent", 1),
])
def test_synonymous_focal_variant_with_an_independent_rna_background(background, protein, annotation_status, count):
    from varcode import Variant
    from isovar import ProteinSequence, export_protein_hypotheses
    from isovar.allele_read import AlleleRead
    from isovar.isovar_result import IsovarResult
    from isovar.read_evidence import ReadEvidence
    from isovar.reference_context import ReferenceContext
    from isovar.translation import Translation
    from isovar.variant_orf import VariantORF
    from isovar.variant_sequence import VariantSequence
    from isovar.protein_reporting import protein_report_record
    from .test_protein_hypotheses import TRANSCRIPT

    variant = Variant("1", 110, "A", "G", normalize_contig_names=False)
    prefix = "CCATG" + background + "AA"
    suffix = "GGCTAA"
    read = AlleleRead(prefix=prefix, allele="G", suffix=suffix, name="original",
                      source_alignments=((("rg", "original", 0), (0, 100, "17M", False)),))
    sequence = VariantSequence(prefix, "G", suffix, [read])
    context = ReferenceContext(strand="+", sequence_before_variant_locus="CCATGGCTAA",
        sequence_at_variant_locus="A", sequence_after_variant_locus=suffix, offset_to_first_complete_codon=2,
        contains_start_codon=True, overlaps_start_codon=True, contains_five_prime_utr=True,
        amino_acids_before_variant="MA", variant=variant, transcripts=[TRANSCRIPT])
    orf = VariantORF(cdna_sequence=sequence.sequence, offset_to_first_complete_codon=2,
        variant_cdna_interval_start=10, variant_cdna_interval_end=11,
        reference_cdna_sequence_before_variant="CCATGGCTAA", reference_cdna_sequence_after_variant=suffix,
        num_mismatches_before_variant=int(background != "GCT"), num_mismatches_after_variant=0)
    translation = Translation(amino_acids=protein, contains_mutation=True, mutation_start_idx=2, mutation_end_idx=3,
        ends_with_stop_codon=True, frameshift=False, untrimmed_variant_sequence=sequence,
        reference_context=context, variant_orf=orf)
    evidence = ReadEvidence(trimmed_base1_start=110, trimmed_ref="A", trimmed_alt="G", ref_reads=[],
                            alt_reads=[read], other_reads=[])
    result = IsovarResult(variant=variant, read_evidence=evidence, predicted_effect=None,
                          sorted_protein_sequences=[ProteinSequence.from_translations([translation])])
    exported = export_protein_hypotheses([result], sample_id="s", source="rna")
    product, = exported["events"][0]["protein_hypotheses"]
    summary = protein_report_record(product, "synonymous")
    assert summary["comparison_status"] == annotation_status
    assert summary["changed_residue_count"] == count
    assert summary["focal_comparison_status"] == "reference_identical"
    assert summary["window_fragments"] == summary["single_read_fragments"] == 1
    comparison = product["translations"][0]["reference_comparison"]
    assert comparison["reference_allele_counterfactual_protein"] == protein
    assert comparison["amino_acid_origins"].get("1") == ("unexplained" if background == "GAT" else None)


def test_conflicting_reference_models_remain_unresolved():
    from isovar.protein_reporting import protein_comparison_status
    translations = [dict(reference_comparison=dict(annotation_comparison=dict(status=state),
                        focal_allele_comparison=dict(status=state)))
                    for state in ("reference_identical", "amino_acid_divergent")]
    protein = dict(translations=translations)
    assert protein_comparison_status(protein) == protein_comparison_status(protein, focal=True) == "unresolved"
