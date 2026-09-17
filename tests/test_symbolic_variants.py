"""Symbolic SVs must never enter literal-allele coordinate/translation paths."""
from types import SimpleNamespace

import pytest
from varcode import StructuralVariant, Variant

from isovar import ProteinSequenceCreator, ReadCollector, run_isovar
from isovar.reference_context_helpers import reference_contexts_for_variant
from isovar.variant_helpers import (
    base0_interval_for_variant, base0_interval_for_variant_fields,
    interbase_range_affected_by_variant_on_transcript, trim_variant, trim_variant_fields,
)


@pytest.mark.parametrize("kind", ["DEL", "DUP", "INV", "INS", "BND"])
@pytest.mark.parametrize("legacy_placeholder", [False, True])
def test_structural_variant_is_rejected_before_interpreting_alleles(kind, legacy_placeholder):
    variant = StructuralVariant("1", 100, kind, end=200, ref="C", genome="GRCh38")
    if legacy_placeholder:
        # Older Varcode exposed a literal-looking ALT for structural events.
        # Preserve that guard without requiring current Varcode to do so.
        variant = SimpleNamespace(is_structural=True, sv_type=kind,
                                  start=100, ref="C", alt="A")
    with pytest.raises(ValueError, match="structural variants"):
        trim_variant(variant)
    with pytest.raises(ValueError, match="structural variants"):
        base0_interval_for_variant(variant)
    with pytest.raises(ValueError, match="structural variants"):
        interbase_range_affected_by_variant_on_transcript(variant, None)
    with pytest.raises(ValueError, match="structural variants"):
        ReadCollector().read_evidence_for_variant(variant, None)
    with pytest.raises(ValueError, match="structural variants"):
        ProteinSequenceCreator().translate_variant_reads(variant, [])
    with pytest.raises(ValueError, match="structural variants"):
        ProteinSequenceCreator().variant_sequences_from_reads(variant, [])
    with pytest.raises(ValueError, match="structural variants"):
        reference_contexts_for_variant(variant)
    # Fail before opening this nonexistent alignment, also for one-shot inputs.
    with pytest.raises(ValueError, match="structural variants"):
        run_isovar(iter([variant]), "/nonexistent/unsupported-sv.bam")


@pytest.mark.parametrize("alt", ["<DEL>", "<DUP:TANDEM>", "C[2:200[", "]2:200]C", "C.", ".C", "*", "."])
def test_symbolic_fields_are_rejected(alt):
    variant = SimpleNamespace(start=100, ref="C", alt=alt)
    for function in (trim_variant, base0_interval_for_variant):
        with pytest.raises(ValueError, match="literal nucleotide"):
            function(variant)
    for function in (trim_variant_fields, base0_interval_for_variant_fields):
        with pytest.raises(ValueError, match="literal nucleotide"):
            function(100, "C", alt)


@pytest.mark.parametrize("ref,alt,interval", [
    ("C", "A", (99, 100)),
    ("C" + "G" * 120, "C", (100, 220)),
    ("C", "C" + "G" * 120, (100, 100)),
])
def test_literal_alleles_including_long_indels_still_work(ref, alt, interval):
    assert base0_interval_for_variant(Variant("1", 100, ref, alt, genome="GRCh38")) == interval
