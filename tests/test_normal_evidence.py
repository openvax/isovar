"""Reference confidence is explicit; neither missing calls nor caller copies provide it."""

import pysam
import pytest
from varcode import Variant

from isovar.normal_evidence import (
    NormalEvidencePolicy, assess_normal_record, normal_source_state, reconcile_normal_states,
)


def record(*, genotype=(0, 0), alleles=("A", "<NON_REF>"), start=99, stop=110,
           quality=40, depth=12, minimum_depth=10, passed=True):
    header = pysam.VariantHeader()
    header.contigs.add("chr1", length=1000)
    header.add_sample("normal.COLUMN")
    for name, number, kind in (("GT", 1, "String"), ("GQ", 1, "Integer"),
                               ("DP", 1, "Integer"), ("MIN_DP", 1, "Integer")):
        header.formats.add(name, number, kind, name)
    header.info.add("END", 1, "Integer", "End of reference block")
    row = header.new_record(contig="chr1", start=start, stop=stop, alleles=alleles)
    row.samples["normal.COLUMN"]["GT"] = genotype
    if quality is not None:
        row.samples["normal.COLUMN"]["GQ"] = quality
    if depth is not None:
        row.samples["normal.COLUMN"]["DP"] = depth
    if minimum_depth is not None:
        row.samples["normal.COLUMN"]["MIN_DP"] = minimum_depth
    if passed:
        row.filter.add("PASS")
    return row


TARGET = Variant("1", 105, "A", "G")
POLICY = NormalEvidencePolicy()


def assess(row, target=TARGET):
    return assess_normal_record(row, "normal.COLUMN", target, POLICY)


def state(rows, *, confidence=True):
    return normal_source_state([assess(row) for row in rows], (104, 105), reference_confidence=confidence)


def test_reference_block_spans_internal_locus():
    assert state([record()]) == "confidently_reference"


@pytest.mark.parametrize("kwargs", [dict(genotype=(None, None)), dict(genotype=(0, None)),
    dict(quality=None), dict(quality=19), dict(minimum_depth=None), dict(minimum_depth=0), dict(passed=False)])
def test_no_call_or_missing_confidence_never_becomes_reference(kwargs):
    evidence = assess(record(**kwargs))
    assert evidence["reference_interval"] is None and evidence["rejection_reasons"]
    assert state([record(**kwargs)]) == "uncalled"


def test_absence_states_remain_distinct():
    assert state([], confidence=True) == "uncalled"
    assert state([], confidence=False) == "variant_only_absence"
    assert reconcile_normal_states([]) == "unavailable"


def test_multiallelic_genotype_does_not_call_every_alt():
    row = record(genotype=(0, 2), alleles=("A", "G", "T"), start=104, stop=105)
    evidence = assess(row)
    assert evidence["called_alternate_indexes"] == [2]
    assert not evidence["target_allele_called"]
    assert state([row]) == "other_allele_observed"
    assert assess(row, Variant("1", 105, "A", "T"))["target_allele_called"]


def test_reference_alt_conflicts_are_preserved():
    alt = record(genotype=(0, 1), alleles=("A", "G"), start=104, stop=105)
    assert state([alt]) == "alternate_observed"
    assert state([record(), alt]) == "conflicting_calls"
    assert reconcile_normal_states(["confidently_reference", "alternate_observed"]) == "conflicting_calls"
    assert reconcile_normal_states(["confidently_reference", "other_allele_observed"]) == "conflicting_calls"


def test_duplicate_processing_products_do_not_outvote_conflicting_evidence():
    assert reconcile_normal_states(["confidently_reference"] * 9 + ["alternate_observed"]) == "conflicting_calls"


def test_structural_variant_reference_genotype_does_not_call_every_base_reference():
    assert assess(record(alleles=("A", "<DEL>")))["reference_interval"] is None


def test_whole_query_interval_requires_contiguous_reference_blocks():
    rows = [assess(record(start=99, stop=104)), assess(record(start=105, stop=110))]
    assert normal_source_state(rows, (103, 107), reference_confidence=True) == "uncalled"
    rows.append(assess(record(start=104, stop=105)))
    assert normal_source_state(rows, (103, 107), reference_confidence=True) == "confidently_reference"


def test_incompatible_explicit_reference_bases_do_not_establish_confidence():
    evidence = assess(record(start=104, stop=105, alleles=("C", "<NON_REF>")))
    assert evidence["rejection_reasons"] == ["reference_allele_mismatch"]
    assert evidence["reference_interval"] is None


@pytest.mark.parametrize("sample_filter,accepted", [("PASS", True), (("PASS",), True),
                                                  (("LowGQX", "HighDP"), False)])
def test_sample_filters_support_vcf_scalar_and_vector_fields(sample_filter, accepted):
    row = record()
    row.header.formats.add("FT", 1 if isinstance(sample_filter, str) else ".", "String", "Genotype filter")
    row.samples["normal.COLUMN"]["FT"] = sample_filter
    assert assess(row)["accepted"] is accepted


def test_strelka_indel_depth_is_declared_and_does_not_replace_block_minimum():
    row = record(genotype=(0, 1), start=104, stop=106, alleles=("AG", "A"), depth=None)
    row.header.formats.add("DPI", 1, "Integer", "Depth at the position preceding the indel")
    row.samples["normal.COLUMN"]["DPI"] = 12
    target = Variant("1", 105, "AG", "A")
    assert not assess_normal_record(row, "normal.COLUMN", target, POLICY)["accepted"]
    policy = NormalEvidencePolicy(indel_depth_field="DPI")
    result = assess_normal_record(row, "normal.COLUMN", target, policy)
    assert result["accepted"] and result["target_allele_called"]
    assert (result["depth_field"], result["depth"]) == ("DPI", 12)
    block = record(minimum_depth=None)
    block.header.formats.add("DPI", 1, "Integer", "Indel depth")
    block.samples["normal.COLUMN"]["DPI"] = 100
    assert not assess_normal_record(block, "normal.COLUMN", TARGET, policy)["accepted"]


@pytest.mark.parametrize("value", [-1, True, 1.5, float("nan"), "20"])
def test_thresholds_require_nonnegative_integers(value):
    with pytest.raises(ValueError, match="nonnegative integers"):
        NormalEvidencePolicy(minimum_depth=value)
