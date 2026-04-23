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
variants on those contigs. It surfaces transcript-relative edits observed in
the assembled cDNA, but it does not yet discover germline SNPs as separate
known variants.
"""

from .transcript_edit_helpers import (
    categorize_transcript_assembly_edits_from_translation,
    transcript_assembly_edit_sort_key,
)


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

        _, translation, matched_transcript = entry
        categorized_edits = categorize_transcript_assembly_edits_from_translation(
            translation,
            matched_transcript,
        )
        transcript_assembly_edits = []
        for category in ("known_somatic", "known_germline", "unexplained"):
            transcript_assembly_edits.extend(categorized_edits[category])
        transcript_assembly_edits.sort(key=transcript_assembly_edit_sort_key)
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
        edits = tuple(
            transcript_assembly_edit.edit
            for transcript_assembly_edit in transcript_assembly_edits
        )
        return MutantTranscript(
            reference_transcript=matched_transcript,
            edits=edits,
            cdna_sequence=translation.cdna_sequence,
            mutant_protein_sequence=translation.amino_acids,
            annotator_name="isovar",
            evidence=evidence)
