"""Vaccine context is useful only when it overlaps the edit and retains support."""

from itertools import permutations
from types import SimpleNamespace

import pytest
from varcode import Variant

from isovar.cli.protein_sequence_args import (
    make_protein_sequences_arg_parser, protein_sequence_creator_from_args,
)
from isovar.protein_sequence import ProteinSequence
from isovar.protein_sequence_creator import ProteinSequenceCreator
from isovar.protein_sequence_helpers import mutant_peptide_window_count, sort_protein_sequences
from isovar.read_evidence import ReadEvidence
from isovar.variant_sequence import VariantSequence
from isovar.variant_sequence_creator import VariantSequenceCreator

from .test_protein_sequences import _reads, _translation_with_reads
from .test_translation_regressions import _context, _genomic_allele, _read


def protein(length, support, start=None, end=None, contains=True, letter="A", source_read_count=1):
    start = length // 2 if start is None else start
    end = start + 1 if end is None else end
    translation = _translation_with_reads(
        letter * length, _reads(letter, support, source_read_count), "AAACTTT")
    translation.mutation_start_idx = start
    translation.mutation_end_idx = end
    translation.contains_mutation = contains
    return ProteinSequence.from_translations([translation])


def test_window_count_against_exhaustive_independent_enumeration():
    # Every valid half-open edit interval, including deletion junctions and
    # terminal edits, for all lengths through the default target. Enumerate
    # starts independently instead of repeating the production formula.
    for length in range(51):
        candidate = SimpleNamespace(amino_acids="A" * length, contains_mutation=True)
        for peptide_length in (1, 2, 3, 8, 9, 15, 25, 30, 49, 51):
            windows = range(max(0, length - peptide_length + 1))
            for start in range(length + 1):
                for end in range(start, length + 1):
                    candidate.mutation_start_idx, candidate.mutation_end_idx = start, end
                    expected = sum(
                        (i < start < i + peptide_length) if start == end else
                        any(start <= residue < end for residue in range(i, i + peptide_length))
                        for i in windows)
                    assert mutant_peptide_window_count(candidate, peptide_length) == expected
        candidate.contains_mutation = False
        assert mutant_peptide_window_count(candidate, 25) == 0


@pytest.mark.parametrize("length,start,end,peptide,expected", [
    (49, 24, 25, 25, 25), (49, 0, 1, 25, 1), (49, 48, 49, 25, 1),
    (25, 12, 13, 25, 1), (24, 12, 13, 25, 0),
    (48, 24, 24, 25, 24), (49, 24, 24, 1, 0),
    (51, 24, 27, 25, 27), (49, 24, 49, 25, 25),
])
def test_window_examples(length, start, end, peptide, expected):
    candidate = protein(length, 2, start, end)
    assert mutant_peptide_window_count(candidate, peptide) == expected


@pytest.mark.parametrize("preference,fraction,expected_length", [
    ("balanced", .9, 49), ("balanced", .95, 9),
    ("balanced", 1, 9), ("balanced", 0, 49),
    ("support", .9, 9), ("context", 1, 49),
])
def test_support_context_tradeoff(preference, fraction, expected_length):
    candidates = [protein(9, 121), protein(49, 111)]
    for order in permutations(candidates):
        ranked = sort_protein_sequences(order, preference, 25, fraction)
        assert len(ranked) == 2
        assert len(ranked[0]) == expected_length


def test_budget_boundary_and_raw_alignment_counts_do_not_inflate_support():
    short, boundary, weak = protein(9, 100), protein(25, 85), protein(49, 84, source_read_count=10)
    assert weak.num_supporting_fragments == 84
    ranked = sort_protein_sequences([weak, short, boundary], "balanced")
    assert ranked == [boundary, short, weak]
    assert sort_protein_sequences([weak, short, boundary], "context")[0] is weak


@pytest.mark.parametrize("fraction,best,retained", [(0.14, 50, 7), (0.28, 25, 7), (0.56, 25, 14)])
def test_fraction_boundary_is_not_rejected_by_float_multiplication(fraction, best, retained):
    short, longer = protein(9, best), protein(25, retained)
    assert sort_protein_sequences([short, longer], "balanced", min_support_fraction=fraction)[0] is longer


def test_nonmutant_candidates_neither_win_nor_set_support_budget():
    mutant = protein(25, 2)
    nonmutant = protein(49, 1000, contains=False)
    assert sort_protein_sequences([nonmutant, mutant], "balanced")[0] is mutant
    assert sort_protein_sequences([nonmutant, mutant], "context")[0] is mutant
    assert sort_protein_sequences([nonmutant], "balanced") == [nonmutant]
    assert sort_protein_sequences([], "balanced") == []


def test_partial_context_is_preferred_only_within_budget():
    short, partial, weak = protein(9, 100), protein(24, 85), protein(25, 84)
    assert sort_protein_sequences([weak, short, partial], "balanced") == [partial, short, weak]


def test_length_without_extra_mutant_windows_does_not_beat_support():
    long_offcenter = protein(49, 90, start=0)
    complete = protein(25, 100)
    assert mutant_peptide_window_count(long_offcenter, 25) == 1
    assert sort_protein_sequences([long_offcenter, complete], "balanced")[0] is complete


@pytest.mark.parametrize("preference", ["balanced", "context"])
def test_content_ties_are_independent_of_input_order(preference):
    candidates = [protein(49, 3, letter=letter) for letter in "ACG"]
    for order in permutations(candidates):
        assert [p.amino_acids for p in sort_protein_sequences(order, preference)] == [
            "G" * 49, "C" * 49, "A" * 49]


@pytest.mark.parametrize("value", [-1, 1.01, float("nan"), float("inf"), -float("inf"), True, "0.9", None, 10 ** 1000])
def test_invalid_support_fractions(value):
    with pytest.raises(ValueError, match="fraction"):
        ProteinSequenceCreator(min_protein_sequence_support_fraction=value)
    with pytest.raises(ValueError, match="fraction"):
        sort_protein_sequences([], min_support_fraction=value)


@pytest.mark.parametrize("value", [0, -1, 2.5, True, "25", None])
def test_invalid_peptide_lengths(value):
    with pytest.raises(ValueError, match="positive integer"):
        ProteinSequenceCreator(protein_context_peptide_length=value)
    with pytest.raises(ValueError, match="positive integer"):
        mutant_peptide_window_count(protein(25, 2), value)


@pytest.mark.parametrize("value", [-1, 3.5, True, "49"])
def test_invalid_context_lengths(value):
    with pytest.raises(ValueError, match="non-negative integer"):
        ProteinSequenceCreator(protein_sequence_length=value)


@pytest.mark.parametrize("floor", [-1, 2.5, 2.0, True, "5", None, float("nan"), float("inf")])
@pytest.mark.parametrize("creator_class", [ProteinSequenceCreator, VariantSequenceCreator])
def test_invalid_absolute_coverage_floors_are_rejected_early(floor, creator_class):
    with pytest.raises(ValueError, match="min_variant_sequence_coverage must be a non-negative integer"):
        creator_class(min_variant_sequence_coverage=floor)


def test_explicit_zero_absolute_floor_remains_supported():
    assert ProteinSequenceCreator(min_variant_sequence_coverage=0).min_variant_sequence_coverage == 0
    assert VariantSequenceCreator(min_variant_sequence_coverage=0).min_variant_sequence_coverage == 0


@pytest.mark.parametrize("value", ["longest", None, "BALANCED", 1])
def test_invalid_preference(value):
    with pytest.raises(ValueError, match="preference"):
        ProteinSequenceCreator(protein_sequence_preference=value)


@pytest.mark.parametrize("peptide", [1, 2, 9, 15, 25, 30, 100])
def test_automatic_context_target(peptide):
    creator = ProteinSequenceCreator(protein_context_peptide_length=peptide)
    assert creator.protein_sequence_length == 2 * peptide - 1
    assert creator.protein_sequence_preference == "balanced"
    assert creator.min_protein_sequence_support_fraction == .85
    assert ProteinSequenceCreator(protein_sequence_length=35,
                                  protein_context_peptide_length=peptide).protein_sequence_length == 35


def test_context_ladder_is_bounded_and_includes_requested_target():
    assert ProteinSequenceCreator().candidate_context_lengths() == [20, 25, 29, 33, 37, 41, 45, 49]
    for target in list(range(101)) + [1000, 1000000]:
        for peptide in (1, 9, 25, 100):
            creator = ProteinSequenceCreator(protein_sequence_length=target,
                                             protein_context_peptide_length=peptide)
            lengths = creator.candidate_context_lengths()
            assert lengths == sorted(set(lengths))
            assert len(lengths) <= 8
            assert max(lengths) == target
            if target > peptide:
                assert peptide in lengths
                assert min(20, target) in lengths
            else:
                assert lengths == [target]
            legacy = ProteinSequenceCreator(protein_sequence_length=target, protein_sequence_preference="support")
            assert legacy.candidate_context_lengths() == [target]
            assert legacy._cdna_sequence_length == 3 * target + 3


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("phase", [0, 1, 2])
@pytest.mark.parametrize("assembly", [False, True])
@pytest.mark.parametrize("peptide_length", [9, 15, 25, 30, 40])
def test_default_centers_single_residue_in_every_codon_phase(strand, phase, assembly, peptide_length, monkeypatch):
    variant = Variant("1", 100, _genomic_allele("G", strand), _genomic_allele("C", strand), "GRCh38")
    prefix, suffix = "ACG" * 40 + "A" * phase, "G" * 120
    reads = [_read(prefix, "C", suffix, str(i), strand) for i in range(3)]
    context = _context(variant, strand, prefix, suffix)
    monkeypatch.setattr("isovar.protein_sequence_creator.reference_contexts_for_variant", lambda *a, **k: [context])
    creator = ProteinSequenceCreator(variant_sequence_assembly=assembly, protein_context_peptide_length=peptide_length)
    evidence = ReadEvidence(ref_reads=[], alt_reads=reads, other_reads=[],
                            trimmed_base1_start=100, trimmed_ref=variant.ref, trimmed_alt=variant.alt)
    top, = creator.sorted_protein_sequences_for_variant(variant, evidence)
    assert len(top) == 2 * peptide_length - 1
    assert (top.mutation_start_idx, top.mutation_end_idx) == (peptide_length - 1, peptide_length)
    assert mutant_peptide_window_count(top, peptide_length) == peptide_length
    assert top.num_supporting_fragments == 3
    assert not top.frameshift
    assert not top.ends_with_stop_codon
    assert top.amino_acids[peptide_length - 1] == ("R", "T", "N")[phase]  # CGG / ACG / AAC
    # The same read objects support multiple contexts without multiplication.
    assert top.supporting_reads == set(reads)


def test_context_ladder_keeps_short_orf_match_without_relaxing_mismatch_filter(monkeypatch):
    variant = Variant("1", 100, "G", "C", "GRCh38")
    # Far-upstream mismatches reject the long context but not the central one.
    reads = [_read("T" * 44 + "A" * 31, "C", "G" * 75, str(i), "+") for i in range(3)]
    context = _context(variant, "+", "A" * 75, "G" * 75)
    monkeypatch.setattr("isovar.protein_sequence_creator.reference_contexts_for_variant", lambda *a, **k: [context])
    creator = ProteinSequenceCreator(max_transcript_mismatches=0)
    assert creator.translate_variant_reads(variant, iter(reads))
    assert creator.translate_variant_reads(variant, iter(())) == []
    legacy = ProteinSequenceCreator(max_transcript_mismatches=0, protein_sequence_preference="support")
    assert legacy.translate_variant_reads(variant, reads) == []
    for translation in creator.translate_variant_reads(variant, reads):
        assert translation.num_mismatches_before_variant == 0
        assert len(translation.reference_cdna_sequence_before_variant) >= 10


def test_exact_rna_deduplication_retains_original_reads_and_conflicting_haplotypes(monkeypatch):
    reads = [_read("AAA", "C", "GGG", str(i), "+") for i in range(2)]
    calls = []

    def candidates(self, variant, original_reads):
        assert original_reads == reads
        calls.append(self.preferred_sequence_length)
        return [VariantSequence("AAA", "C", "GGG", [reads[len(calls) % 2]]),
                VariantSequence("TAA", "C", "GGG", [reads[0]])]

    monkeypatch.setattr("isovar.variant_sequence_creator.VariantSequenceCreator.reads_to_variant_sequences", candidates)
    sequences = ProteinSequenceCreator().variant_sequences_from_reads(None, iter(reads))
    assert len(calls) == 8
    assert len(sequences) == 2
    assert sequences[0].reads == frozenset(reads)
    assert sequences[1].reads == frozenset([reads[0]])
    assert ProteinSequenceCreator().variant_sequences_from_reads(None, []) == []


def test_support_preference_preserves_original_single_scale_candidate_order(monkeypatch):
    reads = [_read("AAA", "C", "GGG", str(i), "+") for i in range(2)]
    original = [VariantSequence("TAA", "C", "GGG", reads), VariantSequence("AAA", "C", "GGG", reads)]
    monkeypatch.setattr("isovar.variant_sequence_creator.VariantSequenceCreator.reads_to_variant_sequences",
                        lambda *a: original)
    creator = ProteinSequenceCreator(protein_sequence_preference="support", protein_sequence_length=20)
    assert creator.variant_sequences_from_reads(None, reads) is original


@pytest.mark.parametrize("strand", ["+", "-"])
def test_default_context_respects_real_stop_and_absolute_coverage(strand, monkeypatch):
    variant = Variant("1", 100, _genomic_allele("G", strand), _genomic_allele("C", strand), "GRCh38")
    prefix, suffix = "ACG" * 40 + "AA", "GGG" * 5 + "TAA" + "GGG" * 30
    reads = [_read(prefix, "C", suffix, str(i), strand) for i in range(2)]
    context = _context(variant, strand, prefix, suffix)
    monkeypatch.setattr("isovar.protein_sequence_creator.reference_contexts_for_variant", lambda *a, **k: [context])
    creator = ProteinSequenceCreator()
    evidence = ReadEvidence(ref_reads=[], alt_reads=reads, other_reads=[],
                            trimmed_base1_start=100, trimmed_ref=variant.ref, trimmed_alt=variant.alt)
    top, = creator.sorted_protein_sequences_for_variant(variant, evidence)
    assert top.amino_acids == "T" * 24 + "N" + "G" * 5
    assert top.ends_with_stop_codon
    assert mutant_peptide_window_count(top, 25) == 6
    assert creator.translate_variant_reads(variant, reads[:1]) == []


def test_output_cap_applies_after_new_ranking(monkeypatch):
    short, long = protein(9, 100), protein(49, 90)
    monkeypatch.setattr(ProteinSequenceCreator, "translate_variant_reads", lambda *a, **k: short.translations + long.translations)
    evidence = ReadEvidence(ref_reads=[], alt_reads=[], other_reads=[],
                            trimmed_base1_start=100, trimmed_ref="G", trimmed_alt="C")
    assert ProteinSequenceCreator().sorted_protein_sequences_for_variant(None, evidence)[0].amino_acids == long.amino_acids
    assert len(ProteinSequenceCreator(max_protein_sequences_per_variant=2).sorted_protein_sequences_for_variant(
        None, evidence)) == 2


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("peptide_length", [15, 25, 30])
@pytest.mark.parametrize("floor", [2, 3, 5, 101])
def test_adaptive_context_obeys_absolute_floor_despite_high_compatible_support(
        strand, peptide_length, floor, monkeypatch):
    """100 compatible names must not make a twice-covered tail look deep."""
    variant = Variant("1", 100, _genomic_allele("G", strand), _genomic_allele("C", strand), "GRCh38")
    prefix, suffix = "ACG" * 40 + "AA", "G" * 120
    reads = [_read(prefix, "C", suffix, str(i), strand) for i in range(2)]
    reads += [_read(prefix[-38:], "C", suffix[:38], str(i), strand) for i in range(2, 100)]
    context = _context(variant, strand, prefix, suffix)
    monkeypatch.setattr("isovar.protein_sequence_creator.reference_contexts_for_variant", lambda *a, **k: [context])
    creator = ProteinSequenceCreator(protein_context_peptide_length=peptide_length,
                                     min_variant_sequence_coverage=floor)
    evidence = ReadEvidence(ref_reads=[], alt_reads=reads, other_reads=[],
                            trimmed_base1_start=100, trimmed_ref=variant.ref, trimmed_alt=variant.alt)
    proteins = creator.sorted_protein_sequences_for_variant(variant, evidence)
    if floor > 100:
        assert proteins == []
        return
    top, = proteins
    assert top.num_supporting_fragments == 100
    assert len(top) == (2 * peptide_length - 1 if floor == 2 else 25)
    assert all(t.untrimmed_variant_sequence.min_coverage() >= floor for t in top.translations)
    assert len(evidence.alt_reads) == len(evidence.alt_read_names) == 100
    # K=30/floor=5 genuinely has no full 30mer; never pad or relax the floor.
    assert mutant_peptide_window_count(top, peptide_length) == max(0, len(top) - peptide_length + 1)


def parse_args(*options):
    return make_protein_sequences_arg_parser().parse_args(["--vcf", "dummy.vcf", "--bam", "dummy.bam", *options])


@pytest.mark.parametrize("peptide_length", [15, 25, 30])
@pytest.mark.parametrize("floor", [2, 5, 10])
def test_cli_adaptive_context_and_both_support_controls(peptide_length, floor):
    creator = protein_sequence_creator_from_args(parse_args(
        "--protein-context-peptide-length", str(peptide_length),
        "--min-variant-sequence-coverage", str(floor),
        "--min-protein-sequence-support-fraction", ".85"))
    assert creator.protein_sequence_length == 2 * peptide_length - 1
    assert creator.min_variant_sequence_coverage == floor
    assert creator.min_protein_sequence_support_fraction == .85


def test_cli_defaults_and_explicit_preferences():
    creator = protein_sequence_creator_from_args(parse_args())
    assert creator.protein_sequence_length == 49
    assert creator.protein_context_peptide_length == 25
    assert creator.min_protein_sequence_support_fraction == .85
    assert creator.min_variant_sequence_coverage == 2
    creator = protein_sequence_creator_from_args(parse_args(
        "--protein-context-peptide-length", "30", "--protein-sequence-preference", "context",
        "--min-protein-sequence-support-fraction", ".8", "--count-mismatches-after-variant"))
    assert creator.protein_sequence_length == 59
    assert creator.protein_sequence_preference == "context"
    assert creator.min_protein_sequence_support_fraction == .8
    assert creator.count_mismatches_after_variant
    legacy = protein_sequence_creator_from_args(parse_args(
        "--protein-sequence-length", "20", "--protein-sequence-preference", "support"))
    assert legacy._cdna_sequence_length == 63
    assert legacy.candidate_context_lengths() == [20]


def test_vaxrank_peptide_size_and_explicit_context_override():
    args = parse_args()
    args.vaccine_peptide_length = None
    assert protein_sequence_creator_from_args(args).protein_sequence_length == 49
    args.vaccine_peptide_length = 27
    creator = protein_sequence_creator_from_args(args)
    assert creator.protein_context_peptide_length == 27
    assert creator.protein_sequence_length == 53
    args.protein_sequence_length = 35
    assert protein_sequence_creator_from_args(args).protein_sequence_length == 35
    args.protein_context_peptide_length = 25
    assert protein_sequence_creator_from_args(args).protein_context_peptide_length == 25
    del args.protein_context_peptide_length
    del args.protein_sequence_preference
    del args.min_protein_sequence_support_fraction
    del args.count_mismatches_after_variant
    assert protein_sequence_creator_from_args(args).protein_context_peptide_length == 27
    assert not protein_sequence_creator_from_args(args).count_mismatches_after_variant


def test_translation_cli_uses_same_policy_and_mismatch_options(monkeypatch):
    from isovar.cli import isovar_translations

    args = parse_args("--protein-context-peptide-length", "30", "--count-mismatches-after-variant")
    seen = []
    monkeypatch.setattr(isovar_translations, "read_evidence_generator_from_args", lambda args: [])

    def translate(self, evidence):
        seen.append(self)
        return iter(())

    monkeypatch.setattr(ProteinSequenceCreator, "translate_variants", translate)
    assert list(isovar_translations.translations_generator_from_args(args)) == []
    assert seen[0].protein_sequence_length == 59
    assert seen[0].protein_sequence_preference == "balanced"
    assert seen[0].count_mismatches_after_variant
