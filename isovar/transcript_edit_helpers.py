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
Helpers for extracting transcript-relative edits from local Isovar assemblies.
"""

from difflib import SequenceMatcher

from varcode.mutant_transcript import TranscriptEdit

from .dna import reverse_complement_dna
from .transcript_assembly_edit import TranscriptAssemblyEdit
from .variant_helpers import (
    interbase_range_affected_by_variant_on_transcript,
    trim_variant_fields,
)


class TrimmedVariant(object):
    def __init__(self, start, ref, alt):
        self.start = start
        self.ref = ref
        self.alt = alt

    @property
    def is_insertion(self):
        return len(self.ref) == 0 and len(self.alt) > 0


def trimmed_variant_from_variant(variant):
    start, ref, alt = trim_variant_fields(variant.start, variant.ref, variant.alt)
    return TrimmedVariant(start=start, ref=ref, alt=alt)


def transcript_edit_from_variant(variant, transcript):
    trimmed_variant = trimmed_variant_from_variant(variant)
    cdna_start, cdna_end = interbase_range_affected_by_variant_on_transcript(
        trimmed_variant,
        transcript,
    )

    alt_bases = trimmed_variant.alt
    if getattr(transcript, "strand", None) == "-":
        alt_bases = reverse_complement_dna(alt_bases)

    return TranscriptEdit(
        cdna_start=cdna_start,
        cdna_end=cdna_end,
        alt_bases=alt_bases,
        source_variant=variant,
    )


def transcript_edits_from_sequence_diff(
        reference_sequence,
        observed_sequence,
        start_offset,
        ignore_terminal_deletions=False,
        ignore_terminal_insertions=False):
    """
    Convert a local sequence difference into TranscriptEdit objects.

    This is intentionally conservative about trailing differences since Isovar
    assemblies may stop before the reference context does.
    """
    edits = []

    if len(reference_sequence) == len(observed_sequence):
        mismatch_start = None
        for idx, (ref_base, observed_base) in enumerate(zip(
                reference_sequence, observed_sequence)):
            if ref_base != observed_base:
                if mismatch_start is None:
                    mismatch_start = idx
            elif mismatch_start is not None:
                edits.append(
                    TranscriptEdit(
                        cdna_start=start_offset + mismatch_start,
                        cdna_end=start_offset + idx,
                        alt_bases=observed_sequence[mismatch_start:idx],
                        source_variant=None,
                    )
                )
                mismatch_start = None

        if mismatch_start is not None:
            edits.append(
                TranscriptEdit(
                    cdna_start=start_offset + mismatch_start,
                    cdna_end=start_offset + len(reference_sequence),
                    alt_bases=observed_sequence[mismatch_start:],
                    source_variant=None,
                )
            )
        return tuple(edits)

    matcher = SequenceMatcher(
        a=reference_sequence,
        b=observed_sequence,
        autojunk=False,
    )
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == "equal":
            continue

        at_reference_end = i2 == len(reference_sequence)
        at_observed_end = j2 == len(observed_sequence)

        if (
                tag == "delete" and
                ignore_terminal_deletions and
                at_reference_end and
                at_observed_end):
            continue

        if (
                tag == "insert" and
                ignore_terminal_insertions and
                at_reference_end and
                at_observed_end):
            continue

        edits.append(
            TranscriptEdit(
                cdna_start=start_offset + i1,
                cdna_end=start_offset + i2,
                alt_bases=observed_sequence[j1:j2],
                source_variant=None,
            )
        )
    return tuple(edits)


def unexplained_transcript_edits_from_translation(translation, transcript):
    """
    Diff the matched local cDNA against the reference transcript window.
    """
    focal_edit = transcript_edit_from_variant(
        translation.reference_context.variant,
        transcript,
    )

    reference_prefix = translation.reference_cdna_sequence_before_variant
    observed_prefix = translation.cdna_sequence[:translation.variant_cdna_interval_start]
    prefix_start = focal_edit.cdna_start - len(reference_prefix)
    prefix_edits = transcript_edits_from_sequence_diff(
        reference_sequence=reference_prefix,
        observed_sequence=observed_prefix,
        start_offset=prefix_start,
    )

    reference_suffix = translation.variant_orf.reference_cdna_sequence_after_variant
    observed_suffix = translation.cdna_sequence[translation.variant_cdna_interval_end:]
    suffix_edits = transcript_edits_from_sequence_diff(
        reference_sequence=reference_suffix,
        observed_sequence=observed_suffix,
        start_offset=focal_edit.cdna_end,
        ignore_terminal_deletions=True,
        ignore_terminal_insertions=True,
    )
    return prefix_edits + suffix_edits


def _transcript_assembly_edit(transcript, edit):
    return TranscriptAssemblyEdit(
        transcript_id=getattr(transcript, "id", ""),
        transcript_name=getattr(transcript, "name", ""),
        edit=edit,
    )


def _source_variant_sort_key(variant):
    if variant is None:
        return ("", -1, "", "")
    return (
        getattr(variant, "contig", ""),
        getattr(variant, "start", -1),
        getattr(variant, "ref", ""),
        getattr(variant, "alt", ""),
    )


def transcript_assembly_edit_sort_key(transcript_assembly_edit):
    return (
        transcript_assembly_edit.transcript_id,
        transcript_assembly_edit.transcript_name,
        transcript_assembly_edit.cdna_start,
        transcript_assembly_edit.cdna_end,
        transcript_assembly_edit.alt_bases,
        _source_variant_sort_key(transcript_assembly_edit.source_variant),
    )


def categorize_transcript_assembly_edits_from_translation(translation, transcript):
    known_somatic_edit = _transcript_assembly_edit(
        transcript,
        transcript_edit_from_variant(translation.reference_context.variant, transcript),
    )
    unexplained_edits = tuple(
        _transcript_assembly_edit(transcript, edit)
        for edit in unexplained_transcript_edits_from_translation(
            translation,
            transcript,
        )
    )
    return {
        "known_somatic": (known_somatic_edit,),
        "known_germline": (),
        "unexplained": unexplained_edits,
    }
