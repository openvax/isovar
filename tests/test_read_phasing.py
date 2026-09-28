# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import doctest

import pytest
from types import SimpleNamespace

from varcode import MolecularPhaseResolver, Variant

import isovar.read_phasing
from isovar import IsovarReadPhasing, run_isovar
from isovar.allele_read import AlleleRead
from isovar.isovar_result import IsovarResult
from isovar.phasing import annotate_phased_variants
from isovar.read_collector import ReadCollector
from isovar.read_evidence import ReadEvidence

from .common import eq_
from .genomes_for_testing import grcm38
from .mock_objects import MockAlignmentFile, make_pysam_read
from .testing_helpers import data_path


def _make_result(variant, alt_read_names, phased_partners=(), ref_read_names=()):
    read_evidence = ReadEvidence(
        trimmed_base1_start=variant.start,
        trimmed_ref=variant.ref,
        trimmed_alt=variant.alt,
        ref_reads=[
            AlleleRead(prefix="A", allele=variant.ref, suffix="T", name=name)
            for name in ref_read_names
        ],
        alt_reads=[
            AlleleRead(prefix="A", allele=variant.alt, suffix="T", name=name)
            for name in alt_read_names
        ],
        other_reads=[],
    )
    return IsovarResult(
        variant=variant,
        read_evidence=read_evidence,
        predicted_effect=None,
        phased_variants_in_supporting_reads=set(phased_partners),
    )


def test_has_evidence_true_for_variant_with_alt_reads():
    v = Variant("1", 10, "A", "C", normalize_contig_names=False)
    phasing = IsovarReadPhasing([_make_result(v, alt_read_names={"r1", "r2"})])
    assert phasing.has_evidence(v)


def test_has_evidence_false_when_no_alt_reads():
    v = Variant("1", 10, "A", "C", normalize_contig_names=False)
    phasing = IsovarReadPhasing([_make_result(v, alt_read_names=set())])
    assert not phasing.has_evidence(v)


def test_has_evidence_false_for_unknown_variant():
    v = Variant("1", 10, "A", "C", normalize_contig_names=False)
    unknown = Variant("1", 99, "G", "T", normalize_contig_names=False)
    phasing = IsovarReadPhasing([_make_result(v, alt_read_names={"r1"})])
    assert not phasing.has_evidence(unknown)


def test_partners_in_cis_returns_phased_variants():
    v1 = Variant("1", 10, "A", "C", normalize_contig_names=False)
    v2 = Variant("1", 11, "G", "T", normalize_contig_names=False)
    v3 = Variant("1", 12, "T", "G", normalize_contig_names=False)
    phasing = IsovarReadPhasing([
        _make_result(v1, alt_read_names={"r1"}, phased_partners={v2, v3}),
        _make_result(v2, alt_read_names={"r1"}, phased_partners={v1}),
        _make_result(v3, alt_read_names={"r1"}, phased_partners={v1}),
    ])
    eq_(phasing.partners_in_cis(v1), (v2, v3))
    eq_(phasing.partners_in_cis(v2), (v1,))


def test_partners_in_cis_empty_when_no_partners():
    v = Variant("1", 10, "A", "C", normalize_contig_names=False)
    phasing = IsovarReadPhasing([_make_result(v, alt_read_names={"r1"})])
    eq_(phasing.partners_in_cis(v), ())


def test_partners_in_cis_empty_for_unknown_variant():
    v = Variant("1", 10, "A", "C", normalize_contig_names=False)
    unknown = Variant("2", 50, "A", "G", normalize_contig_names=False)
    phasing = IsovarReadPhasing([_make_result(v, alt_read_names={"r1"})])
    eq_(phasing.partners_in_cis(unknown), ())


def test_repr_includes_variant_count():
    v1 = Variant("1", 10, "A", "C", normalize_contig_names=False)
    v2 = Variant("1", 11, "G", "T", normalize_contig_names=False)
    phasing = IsovarReadPhasing([
        _make_result(v1, alt_read_names={"r1"}),
        _make_result(v2, alt_read_names={"r2"}),
    ])
    assert repr(phasing) == "IsovarReadPhasing(2 variants)"


def test_docstring_examples_run():
    results = doctest.testmod(isovar.read_phasing, verbose=False)
    assert results.failed == 0, (
        "doctest failures in isovar.read_phasing: %d/%d" % (
            results.failed, results.attempted))
    assert results.attempted > 0, "no doctest examples were collected"


def test_partners_in_cis_is_deterministically_ordered():
    anchor = Variant("1", 10, "A", "C", normalize_contig_names=False)
    later = Variant("1", 20, "A", "C", normalize_contig_names=False)
    earlier = Variant("1", 5, "A", "C", normalize_contig_names=False)
    other_chrom = Variant("2", 1, "A", "C", normalize_contig_names=False)
    phasing = IsovarReadPhasing([
        _make_result(
            anchor,
            alt_read_names={"r1"},
            phased_partners={later, earlier, other_chrom}),
    ])
    eq_(phasing.partners_in_cis(anchor), (earlier, later, other_chrom))


# ---- pipeline smoke tests against the b16 RNA fixture --------------------
#
# Runs the full Isovar pipeline on the same b16 BAM/VCF used by
# tests/test_mutant_transcript_source.py and feeds the results through
# IsovarReadPhasing. The b16 fixture has four variants on four different
# chromosomes, so no two variants are ever co-observed on the same read
# -- every result's partner set is empty. These tests therefore verify
# wiring (the adapter consumes real `run_isovar` output without errors
# and routes the fields it claims to route), not non-empty-partner
# behavior. The non-empty-partner path is exercised by the hand-rolled
# tests above and the annotate_phased_variants integration test below.


@pytest.fixture(scope="module")
def b16_results():
    return run_isovar(
        variants=data_path("data/b16.f10/b16.vcf"),
        alignment_file=data_path("data/b16.f10/b16.combined.sorted.bam"))


def test_b16_has_evidence_matches_alt_fragment_count(b16_results):
    phasing = IsovarReadPhasing(b16_results)
    # `has_evidence` is defined as `num_alt_fragments > 0`; verify it
    # against the real IsovarResult.num_alt_fragments rather than the
    # adapter's internal state.
    for result in b16_results:
        eq_(
            phasing.has_evidence(result.variant),
            result.num_alt_fragments > 0)


def test_b16_partners_are_subset_of_input_variants(b16_results):
    phasing = IsovarReadPhasing(b16_results)
    input_variants = {result.variant for result in b16_results}
    for result in b16_results:
        for partner in phasing.partners_in_cis(result.variant):
            assert partner in input_variants, (
                "partner %s not present in the input result set" % (partner,))


def test_b16_partners_match_phased_variants_field(b16_results):
    phasing = IsovarReadPhasing(b16_results)
    for result in b16_results:
        eq_(
            set(phasing.partners_in_cis(result.variant)),
            set(result.phased_variants_in_supporting_reads))


def test_b16_phasing_is_symmetric(b16_results):
    # Symmetry is the contract we promise in the module docstring; this
    # is the integration-level check that it holds end-to-end against
    # real data.
    phasing = IsovarReadPhasing(b16_results)
    for result in b16_results:
        for partner in phasing.partners_in_cis(result.variant):
            assert result.variant in phasing.partners_in_cis(partner), (
                "asymmetric phasing: %s -> %s but not back" % (
                    result.variant, partner))


# ---- integration against the real upstream phasing logic -----------------
#
# `annotate_phased_variants` is the function `run_isovar` uses to
# populate `phased_variants_in_supporting_reads`. The b16 fixture
# doesn't produce co-phased variants, so to verify
# IsovarReadPhasing against the real upstream logic we hand-roll
# IsovarResult objects with overlapping read names, run them through
# `annotate_phased_variants`, and then exercise the adapter on the
# annotated output. This catches drift if either the upstream phasing
# rules or the adapter's interpretation of the field changes.


def test_adapter_consistent_with_annotate_phased_variants():
    v1 = Variant("1", 10, "A", "C", normalize_contig_names=False)
    v2 = Variant("1", 20, "G", "T", normalize_contig_names=False)
    v3 = Variant("1", 30, "T", "A", normalize_contig_names=False)
    # v1 & v2 share reads r12a, r12b; v2 & v3 share read r23.
    # No reads bridge v1 to v3 directly.
    annotated = annotate_phased_variants(
        [
            _make_result(v1, alt_read_names={"r12a", "r12b"}),
            _make_result(v2, alt_read_names={"r12a", "r12b", "r23"}),
            _make_result(v3, alt_read_names={"r23"}),
        ],
        min_shared_fragments_for_phasing=1,
    )
    phasing = IsovarReadPhasing(annotated)

    # Adapter must expose every partner that the upstream pipeline
    # populated, with no additions or omissions.
    for result in annotated:
        eq_(
            set(phasing.partners_in_cis(result.variant)),
            set(result.phased_variants_in_supporting_reads))

    # Spot-check the expected topology: v2 is partnered with both v1
    # and v3; v1 and v3 partner only with v2 (not with each other,
    # since no read spans v1 -> v3).
    eq_(set(phasing.partners_in_cis(v1)), {v2})
    eq_(set(phasing.partners_in_cis(v2)), {v1, v3})
    eq_(set(phasing.partners_in_cis(v3)), {v2})

    # All three have evidence.
    for v in (v1, v2, v3):
        assert phasing.has_evidence(v)


def test_threshold_filter_propagates_through_adapter():
    # The adapter is passive: it surfaces whatever
    # annotate_phased_variants populated. This test asserts that when
    # the upstream threshold filters out a sub-threshold pair, the
    # adapter reports them as un-partnered (i.e. the filter
    # propagates rather than being silently undone).
    v1 = Variant("1", 10, "A", "C", normalize_contig_names=False)
    v2 = Variant("1", 20, "G", "T", normalize_contig_names=False)
    # Only one read bridges v1 and v2 -- below the threshold of 2.
    annotated = annotate_phased_variants(
        [
            _make_result(v1, alt_read_names={"r12", "r1only"}),
            _make_result(v2, alt_read_names={"r12", "r2only"}),
        ],
        min_shared_fragments_for_phasing=2,
    )
    phasing = IsovarReadPhasing(annotated)
    eq_(phasing.partners_in_cis(v1), ())
    eq_(phasing.partners_in_cis(v2), ())


# ---- synthetic somatic + germline phasing scenario -----------------------
#
# Builds a synthetic mini-BAM in-memory and walks the full
# ReadCollector -> annotate_phased_variants -> IsovarReadPhasing
# stack. The scenario mirrors the canonical motivating use case for
# openvax/isovar#182 / openvax/varcode#268: a somatic SNV whose codon
# overlaps a germline SNP, where reading direction alone can't tell
# whether the two variants are in cis or trans -- only RNA-read
# co-occurrence can.
#
# Imagined locus (mouse, b16-adjacent): chr1:1-10 reference is
# `ACCTTGATCG`. A somatic SNV sits at chr1:4 (T>G), and an imagined
# germline SNP sits at chr1:8 (T>A). The reads below are the kind of
# RNA fragments Isovar's read collector would harvest from a tumor
# BAM at this locus.
#
#   ref       :  A C C T T G A T C G
#   pos (1-based) 1 2 3 4 5 6 7 8 9 10
#   somatic@4 :        ^                      T>G
#   germline@8:                    ^          T>A
#
# Three reads carry both alts (cis evidence), two carry only the
# somatic alt, two carry only the germline alt. After phasing with
# the default min_shared_fragments threshold of 2, the somatic and
# germline calls should be reported as partners.


def test_two_input_variants_with_competing_alleles_remain_unphased():
    chromosome = "1"
    somatic = Variant(
        chromosome, 4, "T", "G", grcm38, normalize_contig_names=False)
    germline = Variant(
        chromosome, 8, "T", "A", grcm38, normalize_contig_names=False)

    # Reads cover positions 1-10 (reference_start=0, 10M cigar). The
    # MD tag encodes mismatches relative to the assumed reference
    # `ACCTTGATCG`:
    #   * both alts:    ACCGTGAACG -> mismatches at 0-based 3 and 7 -> MD=3T3T2
    #   * somatic only: ACCGTGATCG -> mismatch at 0-based 3        -> MD=3T6
    #   * germline only: ACCTTGAACG -> mismatch at 0-based 7        -> MD=7T2
    def both_alts(name):
        return make_pysam_read(
            seq="ACCGTGAACG", cigar="10M", mdtag="3T3T2",
            name=name, reference_start=0)

    def somatic_only(name):
        return make_pysam_read(
            seq="ACCGTGATCG", cigar="10M", mdtag="3T6",
            name=name, reference_start=0)

    def germline_only(name):
        return make_pysam_read(
            seq="ACCTTGAACG", cigar="10M", mdtag="7T2",
            name=name, reference_start=0)

    reads = [
        both_alts("frag-both-1"),
        both_alts("frag-both-2"),
        both_alts("frag-both-3"),
        somatic_only("frag-somatic-1"),
        somatic_only("frag-somatic-2"),
        germline_only("frag-germline-1"),
        germline_only("frag-germline-2"),
    ]
    samfile = MockAlignmentFile(references=(chromosome,), reads=reads)
    collector = ReadCollector()

    somatic_evidence = collector.read_evidence_for_variant(
        variant=somatic, alignment_file=samfile)
    germline_evidence = collector.read_evidence_for_variant(
        variant=germline, alignment_file=samfile)

    # Sanity-check the synthetic-data construction itself before
    # exercising the adapter -- if the read collector doesn't pick
    # the alts up as expected, downstream assertions would be
    # misleading.
    eq_(
        somatic_evidence.alt_read_names,
        {"frag-both-1", "frag-both-2", "frag-both-3",
         "frag-somatic-1", "frag-somatic-2"})
    eq_(
        germline_evidence.alt_read_names,
        {"frag-both-1", "frag-both-2", "frag-both-3",
         "frag-germline-1", "frag-germline-2"})

    annotated = annotate_phased_variants([
        IsovarResult(
            variant=somatic,
            read_evidence=somatic_evidence,
            predicted_effect=None),
        IsovarResult(
            variant=germline,
            read_evidence=germline_evidence,
            predicted_effect=None),
    ])
    phasing = IsovarReadPhasing(annotated)

    assert phasing.has_evidence(somatic)
    assert phasing.has_evidence(germline)
    # Both variants are input results, so use the four-combination rule.
    # (3 both, 2 first alone, 2 second alone) is ambiguous despite sharing
    # three alternate fragments. This is not the matched-germline path.
    assert phasing.in_cis(somatic, germline) is None
    eq_(phasing.partners_in_cis(somatic), ())
    eq_(phasing.partners_in_cis(germline), ())
    assert all(result.phase_group_from_supporting_reads is None for result in annotated)


def test_b16_in_cis_is_unknown_across_chromosomes(b16_results):
    """No fragment covers two loci, so no pair may be reported as trans."""
    phasing = IsovarReadPhasing(b16_results)
    variants = [result.variant for result in b16_results]
    for v1 in variants:
        for v2 in variants:
            if v1 != v2:
                assert phasing.in_cis(v1, v2) is None


# ---- in_cis: cis and trans from fragments covering both loci --------------

V1 = Variant("1", 10, "A", "C", normalize_contig_names=False)
V2 = Variant("1", 20, "G", "T", normalize_contig_names=False)


def test_in_cis_true_when_fragments_carry_both_alts():
    phasing = IsovarReadPhasing([
        _make_result(V1, alt_read_names={"f1", "f2"}),
        _make_result(V2, alt_read_names={"f1", "f2"}),
    ])
    assert phasing.in_cis(V1, V2) is True
    assert phasing.in_cis(V2, V1) is True


def test_in_cis_false_when_each_alt_comes_with_the_other_reference():
    phasing = IsovarReadPhasing([
        _make_result(V1, alt_read_names={"f1", "f2"}, ref_read_names={"f3", "f4"}),
        _make_result(V2, alt_read_names={"f3", "f4"}, ref_read_names={"f1", "f2"}),
    ])
    assert phasing.in_cis(V1, V2) is False
    assert phasing.in_cis(V2, V1) is False


def test_a_subclone_nested_in_the_others_clone_is_cis():
    # #393: V1 arose on V2's copy in a subclone of V2's cells. Every V1
    # molecule carries V2; V2's other cells give (V1 ref, V2 alt) fragments,
    # which outnumber the both-alt ones but aren't trans evidence.
    nested = ["n%d" % i for i in range(5)]
    subclone = IsovarReadPhasing([
        _make_result(V1, alt_read_names={"b1", "b2"}, ref_read_names=set(nested)),
        _make_result(V2, alt_read_names={"b1", "b2", *nested}),
    ])
    assert subclone.in_cis(V1, V2) is True
    assert subclone.in_cis(V2, V1) is True
    # Nested the other way round.
    ancestor = IsovarReadPhasing([
        _make_result(V1, alt_read_names={"b1", "b2", *nested}),
        _make_result(V2, alt_read_names={"b1", "b2"}, ref_read_names=set(nested)),
    ])
    assert ancestor.in_cis(V1, V2) is True


def test_in_cis_none_when_only_one_alt_is_seen_without_the_other():
    # V2 in trans, or on V1's copy in a subclone these fragments missed.
    phasing = IsovarReadPhasing([
        _make_result(V1, alt_read_names={"f1", "f2"}),
        _make_result(V2, alt_read_names={"f3"}, ref_read_names={"f1", "f2"}),
    ])
    assert phasing.in_cis(V1, V2) is None


def test_in_cis_none_when_the_smaller_direction_is_below_the_minimum():
    # Ten fragments with V1's alt allele alone, one with V2's: the old summed
    # rule said trans, but one fragment doesn't separate the pair.
    alone = ["a%d" % i for i in range(10)]
    phasing = IsovarReadPhasing([
        _make_result(V1, alt_read_names=set(alone), ref_read_names={"b1"}),
        _make_result(V2, alt_read_names={"b1"}, ref_read_names=set(alone)),
    ])
    assert phasing.in_cis(V1, V2) is None


def _phasing_of_table(both, first, second, neither, **kwargs):
    """Phasing of V1 and V2 over fragments with (alt, alt), (alt, ref), (ref, alt) and (ref, ref)."""
    def names(prefix, count):
        return {"%s%d" % (prefix, i) for i in range(count)}
    together, alone1, alone2, reference = (
        names("c", both), names("a", first), names("b", second), names("z", neither))
    return IsovarReadPhasing([
        _make_result(V1, alt_read_names=together | alone1, ref_read_names=alone2 | reference),
        _make_result(V2, alt_read_names=together | alone2, ref_read_names=alone1 | reference),
    ], **kwargs)


def test_errors_on_the_others_alt_fragments_are_not_cis():
    # #410: V2 unexpressed, or a subclone on the other copy; 1% errors put its
    # alt allele on two of V1's 200 alt fragments.
    assert _phasing_of_table(2, 198, 0, 20).in_cis(V1, V2) is None
    assert _phasing_of_table(2, 200, 1, 100).in_cis(V1, V2) is None
    # Two in 200 is more than errors at a rate of 0.1%.
    assert _phasing_of_table(2, 198, 0, 20, phasing_error_rate=0.001).in_cis(V1, V2) is True


def test_trans_despite_errors_on_the_dominant_copy():
    # V1's copy amplified: its errors give 10 both-alt fragments in 1000.
    assert _phasing_of_table(10, 990, 6, 94).in_cis(V1, V2) is False


def test_a_small_subclone_on_the_other_copy_is_trans():
    assert _phasing_of_table(0, 50, 3, 40).in_cis(V1, V2) is False


def test_all_four_combinations_have_no_single_lineage():
    assert _phasing_of_table(3, 10, 2, 5).in_cis(V1, V2) is None


def test_nested_cis_with_few_reference_fragments():
    # A subclone in a quarter of V1's alt fragments, few (ref, ref) fragments.
    assert _phasing_of_table(3, 10, 0, 5).in_cis(V1, V2) is True
    assert _phasing_of_table(2, 5, 0, 50).in_cis(V1, V2) is True


def test_nested_cis_either_way_round():
    assert _phasing_of_table(3, 0, 10, 8).in_cis(V1, V2) is True
    assert _phasing_of_table(3, 10, 0, 8).in_cis(V1, V2) is True


def test_errors_leak_into_a_combination_from_both_neighbours():
    # Clonal pairs at depth with 1% errors: 100 fragments of each alt allele
    # alone are errors from the both-alt and neither fragments together.
    assert _phasing_of_table(3500, 100, 100, 6300).in_cis(V1, V2) is True
    assert _phasing_of_table(70, 3500, 3500, 3000).in_cis(V1, V2) is False


def test_one_stray_fragment_does_not_veto():
    assert _phasing_of_table(3, 1, 1, 0).in_cis(V1, V2) is True


def test_error_prone_reads_need_their_error_rate():
    # Trans with B's copy barely expressed, from reads with 5% errors: at the
    # default 1% those errors look like cis.
    assert _phasing_of_table(3, 47, 2, 48).in_cis(V1, V2) is True
    assert _phasing_of_table(3, 47, 2, 48, phasing_error_rate=0.05).in_cis(V1, V2) is None


def test_binomial_tail():
    from isovar.phasing import _binomial_tail
    assert _binomial_tail(0, 5, 0.01) == 1.0
    assert _binomial_tail(6, 5, 0.01) == 0.0
    assert abs(_binomial_tail(1, 1, 0.01) - 0.01) < 1e-15
    # 1 - P(0) - P(1) for Binomial(200, 0.01)
    assert abs(_binomial_tail(2, 200, 0.01) - (1 - 0.99 ** 200 - 200 * 0.01 * 0.99 ** 199)) < 1e-12
    assert _binomial_tail(3000, 200000, 0.01) < 1e-90


def test_phasing_rates_are_probabilities():
    for name in ("phasing_error_rate", "max_p_value_for_phasing"):
        for rate in (0, 1, -0.1):
            with pytest.raises(ValueError, match=name):
                IsovarReadPhasing([], **{name: rate})


@pytest.mark.parametrize("rate, p_value, expected", [
    (0.01, 0.05, None),
    (0.001, 0.05, True),
    (0.001, 0.01, None),
])
def test_run_isovar_uses_phasing_thresholds_on_collected_reads(monkeypatch, rate, p_value, expected):
    # Hold translation out of this miniature reference: exercise the actual
    # collector -> result -> phase annotation path on the #410 fragment table.
    monkeypatch.setattr("isovar.main.top_varcode_effect", lambda **kwargs: None)
    creator = SimpleNamespace(
        sorted_protein_sequences_for_variant=lambda **kwargs: [],
        settings=lambda: {})
    v1 = Variant("1", 4, "T", "G", grcm38, normalize_contig_names=False)
    v2 = Variant("1", 8, "T", "A", grcm38, normalize_contig_names=False)
    reads = []
    for label, count, sequence, md in (
            ("both", 2, "ACCGTGAACG", "3T3T2"),
            ("first", 198, "ACCGTGATCG", "3T6"),
            ("neither", 20, "ACCTTGATCG", "10")):
        reads.extend(make_pysam_read(
            seq=sequence, cigar="10M", mdtag=md, name="%s%d" % (label, i), reference_start=0)
            for i in range(count))
    thresholds = dict(phasing_error_rate=rate, max_p_value_for_phasing=p_value)
    results = run_isovar(
        [v1, v2], MockAlignmentFile(references=("1",), reads=reads),
        protein_sequence_creator=creator, **thresholds)
    assert [r.num_alt_fragments for r in results] == [200, 2]
    phasing = IsovarReadPhasing(results, **thresholds)
    assert phasing.in_cis(v1, v2) is expected
    assert phasing.partners_in_cis(v1) == ((v2,) if expected else ())
    assert (results[0].phase_group_from_supporting_reads is not None) == (expected is True)


def test_in_cis_none_without_fragments_covering_both_loci():
    phasing = IsovarReadPhasing([
        _make_result(V1, alt_read_names={"f1", "f2"}),
        _make_result(V2, alt_read_names={"f3", "f4"}),
    ])
    assert phasing.partners_in_cis(V1) == ()
    assert phasing.in_cis(V1, V2) is None


def test_in_cis_needs_min_shared_fragments():
    results = [
        _make_result(V1, alt_read_names={"f1", "f2", "f3"}),
        _make_result(V2, alt_read_names={"f1"}, ref_read_names={"f2"}),
    ]
    assert IsovarReadPhasing(results).in_cis(V1, V2) is None
    single = IsovarReadPhasing(
        results[:1] + [_make_result(V2, alt_read_names={"f1"})],
        min_shared_fragments_for_phasing=1)
    assert single.in_cis(V1, V2) is True


def test_matched_germline_edit_in_assembly_is_cis():
    germline = Variant("1", 12, "T", "G", normalize_contig_names=False)
    unseen = Variant("1", 14, "T", "G", normalize_contig_names=False)
    somatic = SimpleNamespace(
        variant=V1, num_alt_fragments=2,
        phased_variants_in_supporting_reads=set(),
        germline_variants_in_top_protein_sequence={germline})
    phasing = IsovarReadPhasing([somatic])
    assert phasing.in_cis(V1, germline) is True
    assert phasing.in_cis(germline, V1) is True
    assert phasing.in_cis(V1, unseen) is None
    assert phasing.in_cis(germline, unseen) is None


def test_varcode_resolver_uses_in_cis_for_trans():
    phasing = IsovarReadPhasing([
        _make_result(V1, alt_read_names={"f1", "f2"}, ref_read_names={"f3", "f4"}),
        _make_result(V2, alt_read_names={"f3", "f4"}, ref_read_names={"f1", "f2"}),
    ])
    assert MolecularPhaseResolver(phasing).in_cis(V1, V2) is False
