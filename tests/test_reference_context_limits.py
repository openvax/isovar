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

import pytest
from varcode import Variant

from isovar.protein_sequence_creator import ProteinSequenceCreator
from isovar.variant_orf_helpers import match_variant_sequence_to_reference_context
from isovar.variant_sequence import VariantSequence
from isovar.variant_sequence_creator import VariantSequenceCreator

from .test_translation_regressions import _context, _genomic_allele, _read


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("assembly", [False, True])
@pytest.mark.parametrize("replacement", [False, True])
@pytest.mark.parametrize("alt_length", [1, 3, 42, 43, 44, 45, 60, 61, 62, 63, 64, 65, 69, 90, 120])
def test_long_alternates_retain_required_context(strand, assembly, replacement, alt_length, monkeypatch):
    """#212: output-length preferences must not erase available ORF evidence."""
    cdna_alt = "A" * alt_length
    cdna_ref = "T" * alt_length if replacement else ""
    variant = Variant("1", 100, _genomic_allele(cdna_ref, strand),
                      _genomic_allele(cdna_alt, strand), "GRCh38")
    prefix, suffix = "ACG" * 4, "G" * 30
    reads = [_read(prefix, cdna_alt, suffix, str(i), strand) for i in range(3)]
    context = _context(variant, strand, prefix, suffix)
    monkeypatch.setattr(
        "isovar.protein_sequence_creator.reference_contexts_for_variant",
        lambda *args, **kwargs: [context])
    creator = ProteinSequenceCreator(variant_sequence_assembly=assembly)

    translation, = creator.translate_variant_reads(variant, reads)

    assert translation.contains_mutation
    assert 0 < len(translation.amino_acids) <= creator.protein_sequence_length
    assert translation.reads == frozenset(reads)
    assert translation.num_mismatches_before_variant == 0
    assert len(translation.reference_cdna_sequence_before_variant) >= 10
    orf = translation.variant_orf
    assert orf.variant_cdna_interval_end - orf.variant_cdna_interval_start == alt_length
    assert orf.cdna_sequence[orf.variant_cdna_interval_start:orf.variant_cdna_interval_end] == cdna_alt
    if alt_length == 45:
        # Prefix trimming leaves three complete ACG codons, followed by
        # fifteen inserted/replaced AAA codons and the retained GGG suffix.
        assert translation.amino_acids == "TTT" + "K" * 15 + "GG"
        assert (translation.mutation_start_idx, translation.mutation_end_idx) == (3, 18)
        assert not translation.frameshift


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("minimum", [3, 10, 25])
@pytest.mark.parametrize("missing_context", [None, "rna", "reference"])
def test_long_insertion_respects_custom_context_minimum(strand, minimum, missing_context, monkeypatch):
    variant = Variant("1", 100, "", _genomic_allele("A" * 69, strand), "GRCh38")
    rna_length = minimum - 1 if missing_context == "rna" else minimum + 2
    reference_length = minimum - 1 if missing_context == "reference" else minimum + 2
    reads = [_read("A" * rna_length, "A" * 69, "G" * 30, str(i), strand) for i in range(3)]
    context = _context(variant, strand, "A" * reference_length, "G" * 30)
    monkeypatch.setattr(
        "isovar.protein_sequence_creator.reference_contexts_for_variant",
        lambda *args, **kwargs: [context])
    creator = ProteinSequenceCreator(min_transcript_prefix_length=minimum)
    translations = creator.translate_variant_reads(variant, reads)

    if missing_context:
        assert translations == []
    else:
        translation, = translations
        assert translation.contains_mutation
        assert len(translation.reference_cdna_sequence_before_variant) >= minimum


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("rna_length", [0, 3, 9, 10, 11, 12])
@pytest.mark.parametrize("reference_length", [0, 3, 9, 10, 11, 12])
@pytest.mark.parametrize("minimum", [0, 3, 10])
def test_match_minimum_uses_actual_shared_prefix(strand, rna_length, reference_length, minimum):
    """#213: enforce the minimum on both sides of prefix normalization."""
    variant = Variant("1", 100, _genomic_allele("G", strand),
                      _genomic_allele("C", strand), "GRCh38")
    read = _read("A" * rna_length, "C", "G" * 20, "read", strand)
    sequence = VariantSequence(read.prefix, read.allele, read.suffix, [read])
    context = _context(variant, strand, "A" * reference_length, "G" * 20)
    orf = match_variant_sequence_to_reference_context(
        sequence, context, min_transcript_prefix_length=minimum, max_transcript_mismatches=0)

    shared_length = min(rna_length, reference_length)
    if shared_length < minimum:
        assert orf is None
    else:
        assert orf is not None
        assert len(orf.reference_cdna_sequence_before_variant) == shared_length
        assert orf.num_mismatches_before_variant == 0


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("core_length", [9, 10, 11])
def test_match_minimum_is_rechecked_after_coverage_trimming(strand, core_length):
    variant = Variant("1", 100, _genomic_allele("G", strand),
                      _genomic_allele("C", strand), "GRCh38")
    long_read = _read("C" + "A" * core_length, "C", "G" * 20, "long", strand)
    short_read = _read("A" * core_length, "C", "G" * 20, "short", strand)
    sequence = VariantSequence(long_read.prefix, long_read.allele, long_read.suffix, [long_read, short_read])
    context = _context(variant, strand, "A" * 12, "G" * 20)

    orf = match_variant_sequence_to_reference_context(
        sequence, context, min_transcript_prefix_length=10,
        max_transcript_mismatches=0, max_trimming_attempts=1)

    if core_length < 10:
        assert orf is None
    else:
        assert orf is not None
        assert len(orf.reference_cdna_sequence_before_variant) == core_length


def test_reference_lookup_budget_is_not_smaller_than_matching_minimum(monkeypatch):
    seen = []

    def lookup(variant, context_size, transcript_id_whitelist):
        seen.append(context_size)
        return []

    monkeypatch.setattr("isovar.protein_sequence_creator.reference_contexts_for_variant", lookup)
    variant = Variant("1", 100, "G", "C", "GRCh38")
    read = _read("A" * 80, "C", "G" * 80, "read", "+")
    ProteinSequenceCreator(min_transcript_prefix_length=70).translate_variant_reads(variant, [read])
    assert len(seen) == 1
    assert seen[0] >= 70


@pytest.mark.parametrize("preferred_length", [0, 1, 3, 62, 63, 64])
@pytest.mark.parametrize("minimum_flank", [0, 3, 10])
def test_sequence_budget_never_creates_negative_flank_limits(preferred_length, minimum_flank):
    creator = VariantSequenceCreator(
        preferred_sequence_length=preferred_length,
        min_flanking_sequence_length=minimum_flank,
        variant_sequence_assembly=True)
    reads = [_read("A" * 20, "C" * 63, "G" * 20, str(i), "+") for i in range(2)]
    variant = Variant("1", 100, "", "C" * 63, "GRCh38")
    sequence, = creator.reads_to_variant_sequences(variant, reads)
    expected_prefix_length = 1 if preferred_length == 64 and minimum_flank == 0 else minimum_flank
    assert sequence.prefix == "A" * expected_prefix_length
    assert sequence.suffix == "G" * minimum_flank
    assert sequence.alt == "C" * 63


@pytest.mark.parametrize("kwargs", [
    {"preferred_sequence_length": -1},
    {"min_flanking_sequence_length": -1},
])
def test_negative_sequence_budgets_are_rejected(kwargs):
    with pytest.raises(ValueError, match="non-negative"):
        VariantSequenceCreator(**kwargs)
