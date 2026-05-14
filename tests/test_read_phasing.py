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
from varcode import Variant

import isovar.read_phasing
from isovar import IsovarReadPhasing, run_isovar
from isovar.allele_read import AlleleRead
from isovar.isovar_result import IsovarResult
from isovar.read_evidence import ReadEvidence

from .common import eq_
from .testing_helpers import data_path


def _make_result(variant, alt_read_names, phased_partners=()):
    read_evidence = ReadEvidence(
        trimmed_base1_start=variant.start,
        trimmed_ref=variant.ref,
        trimmed_alt=variant.alt,
        ref_reads=[],
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


# ---- end-to-end test against the b16 RNA fixture -------------------------
#
# Runs the full Isovar pipeline on the same b16 BAM/VCF used by
# tests/test_varcode_adapter.py. The b16 fixture is small (a handful of
# variants), so partner sets may be empty -- the structural assertions
# below pass either way and would still catch a regression where the
# adapter drops, fabricates, or misroutes phasing data.


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
