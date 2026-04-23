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

from dataclasses import dataclass
from types import SimpleNamespace

from varcode import Variant

from isovar import ProteinSequence, VarcodeAdapter
from isovar.allele_read import AlleleRead
from isovar.isovar_result import IsovarResult
from isovar.read_evidence import ReadEvidence
from isovar.translation import Translation
from isovar.transcript_assembly_edit import TranscriptAssemblyEdit
from isovar.transcript_edit_helpers import (
    categorize_transcript_assembly_edits_from_translation,
)
from isovar.variant_orf import VariantORF

from .common import eq_


@dataclass(frozen=True, order=True)
class DummyTranscript(object):
    id: str
    name: str
    strand: str = "+"
    genomic_offset: int = 100

    def spliced_offset(self, dna_pos):
        return dna_pos - self.genomic_offset


def make_translation(
        variant,
        transcript,
        reference_prefix="TTAA",
        observed_prefix="TGAA",
        observed_alt="G",
        reference_suffix="CCGG",
        observed_suffix="CCAG"):
    cdna_sequence = observed_prefix + observed_alt + observed_suffix
    variant_orf = VariantORF(
        cdna_sequence=cdna_sequence,
        offset_to_first_complete_codon=0,
        variant_cdna_interval_start=len(observed_prefix),
        variant_cdna_interval_end=len(observed_prefix) + len(observed_alt),
        reference_cdna_sequence_before_variant=reference_prefix,
        reference_cdna_sequence_after_variant=reference_suffix,
        num_mismatches_before_variant=1,
        num_mismatches_after_variant=1,
    )
    reference_context = SimpleNamespace(
        variant=variant,
        transcripts=(transcript,),
    )
    return Translation(
        amino_acids="M",
        contains_mutation=True,
        mutation_start_idx=0,
        mutation_end_idx=1,
        ends_with_stop_codon=False,
        frameshift=False,
        untrimmed_variant_sequence=SimpleNamespace(reads=[]),
        reference_context=reference_context,
        variant_orf=variant_orf,
    )


def make_expected_edit(
        transcript,
        cdna_start,
        cdna_end,
        alt_bases,
        source_variant=None):
    from varcode.mutant_transcript import TranscriptEdit

    return TranscriptAssemblyEdit(
        transcript_id=transcript.id,
        transcript_name=transcript.name,
        edit=TranscriptEdit(
            cdna_start=cdna_start,
            cdna_end=cdna_end,
            alt_bases=alt_bases,
            source_variant=source_variant,
        ),
    )


def make_isovar_result(variant, protein_sequence):
    read_evidence = ReadEvidence(
        trimmed_base1_start=variant.start,
        trimmed_ref=variant.ref,
        trimmed_alt=variant.alt,
        ref_reads=[],
        alt_reads=[
            AlleleRead(prefix="A", allele=variant.alt, suffix="T", name="read-1")
        ],
        other_reads=[],
    )
    return IsovarResult(
        variant=variant,
        read_evidence=read_evidence,
        predicted_effect=None,
        sorted_protein_sequences=[protein_sequence],
    )


def test_categorize_transcript_assembly_edits_from_translation():
    variant = Variant("1", 105, "A", "G", normalize_contig_names=False)
    transcript = DummyTranscript(id="tx-1", name="TX1")
    translation = make_translation(variant, transcript)

    categorized_edits = categorize_transcript_assembly_edits_from_translation(
        translation,
        transcript,
    )

    eq_(
        categorized_edits["known_somatic"],
        (make_expected_edit(transcript, 5, 6, "G", variant),),
    )
    eq_(categorized_edits["known_germline"], ())
    eq_(
        categorized_edits["unexplained"],
        (
            make_expected_edit(transcript, 2, 3, "G"),
            make_expected_edit(transcript, 8, 9, "A"),
        ),
    )


def test_categorize_transcript_assembly_edits_ignores_terminal_suffix_gap():
    variant = Variant("1", 105, "A", "G", normalize_contig_names=False)
    transcript = DummyTranscript(id="tx-1", name="TX1")
    translation = make_translation(
        variant,
        transcript,
        reference_suffix="CCGG",
        observed_suffix="CC",
    )

    categorized_edits = categorize_transcript_assembly_edits_from_translation(
        translation,
        transcript,
    )

    eq_(
        categorized_edits["known_somatic"],
        (make_expected_edit(transcript, 5, 6, "G", variant),),
    )
    eq_(categorized_edits["known_germline"], ())
    eq_(categorized_edits["unexplained"], (make_expected_edit(transcript, 2, 3, "G"),))


def test_protein_sequence_exposes_deduplicated_transcript_assembly_edits():
    variant = Variant("1", 105, "A", "G", normalize_contig_names=False)
    transcript = DummyTranscript(id="tx-1", name="TX1")
    translation = make_translation(variant, transcript)
    protein_sequence = ProteinSequence.from_translations([translation, translation])

    eq_(
        protein_sequence.known_somatic_transcript_edits,
        (make_expected_edit(transcript, 5, 6, "G", variant),),
    )
    eq_(protein_sequence.known_germline_transcript_edits, ())
    eq_(
        protein_sequence.unexplained_transcript_edits,
        (
            make_expected_edit(transcript, 2, 3, "G"),
            make_expected_edit(transcript, 8, 9, "A"),
        ),
    )


def test_varcode_adapter_mutant_transcript_includes_unexplained_edits():
    variant = Variant("1", 105, "A", "G", normalize_contig_names=False)
    transcript = DummyTranscript(id="tx-1", name="TX1")
    translation = make_translation(variant, transcript)
    protein_sequence = ProteinSequence.from_translations([translation])
    provider = VarcodeAdapter([make_isovar_result(variant, protein_sequence)])

    mutant_transcript = provider.mutant_transcript(variant, transcript)

    assert mutant_transcript is not None
    eq_(
        tuple(
            (edit.cdna_start, edit.cdna_end, edit.alt_bases, edit.source_variant)
            for edit in mutant_transcript.edits
        ),
        (
            (2, 3, "G", None),
            (5, 6, "G", variant),
            (8, 9, "A", None),
        ),
    )
