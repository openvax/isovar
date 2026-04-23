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

"""
Adapter from IsovarResult collections to varcode's IsovarAssemblyProvider.

This is intentionally conservative: it exposes transcript-keyed assembled
contigs from existing Translation objects and surfaces already-known phased RNA
variants on those contigs. It does not yet discover germline SNPs or unexplained
contig edits directly from the assembled cDNA sequence.
"""

from collections import namedtuple

from .dna import reverse_complement_dna
from .variant_helpers import interbase_range_affected_by_variant_on_transcript


TrimmedVariant = namedtuple(
    "TrimmedVariant",
    ["start", "ref", "alt", "is_insertion"])


class VarcodeAdapter(object):
    """
    Present IsovarResult objects through varcode's IsovarAssemblyProvider
    protocol without changing run_isovar's existing return type.
    """

    def __init__(self, isovar_results):
        self._entries_by_variant_and_transcript = {}
        for result in isovar_results:
            for protein_sequence in result.sorted_protein_sequences:
                for translation in protein_sequence.translations:
                    for transcript in translation.reference_context.transcripts:
                        key = (result.variant, self._transcript_key(transcript))
                        if key not in self._entries_by_variant_and_transcript:
                            self._entries_by_variant_and_transcript[key] = (
                                result,
                                translation,
                                transcript,
                            )

    @staticmethod
    def _transcript_key(transcript):
        return getattr(transcript, "id", transcript)

    @staticmethod
    def _variant_sort_key(variant):
        return (
            variant.contig,
            variant.start,
            variant.ref,
            variant.alt,
        )

    def _entry(self, variant, transcript):
        if transcript is None:
            return None
        return self._entries_by_variant_and_transcript.get((
            variant,
            self._transcript_key(transcript),
        ))

    @staticmethod
    def _trimmed_variant_from_result(result):
        ref = result.read_evidence.trimmed_ref
        alt = result.read_evidence.trimmed_alt
        return TrimmedVariant(
            start=result.read_evidence.trimmed_base1_start,
            ref=ref,
            alt=alt,
            is_insertion=len(ref) == 0 and len(alt) > 0)

    def _transcript_edit(self, result, transcript):
        from varcode.mutant_transcript import TranscriptEdit

        trimmed_variant = self._trimmed_variant_from_result(result)
        try:
            cdna_start, cdna_end = interbase_range_affected_by_variant_on_transcript(
                trimmed_variant,
                transcript)
        except ValueError:
            return None

        alt_bases = trimmed_variant.alt
        if getattr(transcript, "strand", None) == "-":
            alt_bases = reverse_complement_dna(alt_bases)

        return TranscriptEdit(
            cdna_start=cdna_start,
            cdna_end=cdna_end,
            alt_bases=alt_bases,
            source_variant=result.variant)

    def has_contig(self, variant, transcript):
        return self._entry(variant, transcript) is not None

    def variants_in_contig(self, variant, transcript):
        entry = self._entry(variant, transcript)
        if entry is None:
            return ()

        result, _, _ = entry
        variants = [
            phased_variant
            for phased_variant in result.phased_variants_in_protein_sequence
            if self.has_contig(phased_variant, transcript)
        ]
        return tuple(sorted(variants, key=self._variant_sort_key))

    def mutant_transcript(self, variant, transcript):
        entry = self._entry(variant, transcript)
        if entry is None:
            return None

        from varcode.mutant_transcript import MutantTranscript

        result, translation, matched_transcript = entry
        transcript_edit = self._transcript_edit(result, matched_transcript)
        evidence = {
            "num_supporting_reads": sum(
                getattr(read, "source_read_count", 1)
                for read in translation.reads
            ),
            "num_supporting_fragments": len({read.name for read in translation.reads}),
            "num_cdna_mismatches_before_variant": (
                translation.num_mismatches_before_variant),
            "num_cdna_mismatches_after_variant": (
                translation.num_mismatches_after_variant),
        }
        edits = (transcript_edit,) if transcript_edit is not None else ()
        return MutantTranscript(
            reference_transcript=matched_transcript,
            edits=edits,
            cdna_sequence=translation.cdna_sequence,
            mutant_protein_sequence=translation.amino_acids,
            annotator_name="isovar",
            evidence=evidence)
