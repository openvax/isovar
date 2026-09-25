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


def _flanks(translation, transcript):
    """The focal edit and the reference flanks (sequence, transcript offset) around it."""
    focal_edit = transcript_edit_from_variant(
        translation.reference_context.variant,
        transcript,
    )
    reference_prefix = translation.reference_cdna_sequence_before_variant
    reference_suffix = translation.variant_orf.reference_cdna_sequence_after_variant
    return focal_edit, (
        (reference_prefix, focal_edit.cdna_start - len(reference_prefix)),
        (reference_suffix, focal_edit.cdna_end),
    )


def unexplained_transcript_edits_from_translation(translation, transcript):
    """
    Diff the matched local cDNA against the reference transcript window.
    """
    focal_edit, ((reference_prefix, prefix_start), (reference_suffix, suffix_start)) = \
        _flanks(translation, transcript)
    observed_prefix = translation.cdna_sequence[:translation.variant_cdna_interval_start]
    prefix_edits = transcript_edits_from_sequence_diff(
        reference_sequence=reference_prefix,
        observed_sequence=observed_prefix,
        start_offset=prefix_start,
    )

    observed_suffix = translation.cdna_sequence[translation.variant_cdna_interval_end:]
    suffix_edits = transcript_edits_from_sequence_diff(
        reference_sequence=reference_suffix,
        observed_sequence=observed_suffix,
        start_offset=suffix_start,
        ignore_terminal_deletions=True,
        ignore_terminal_insertions=True,
    )
    return prefix_edits + suffix_edits


def _apply_edit(sequence, offset, edit):
    """``sequence`` (starting at transcript ``offset``) with ``edit``, or None outside it."""
    start, end = edit.cdna_start - offset, edit.cdna_end - offset
    if not 0 <= start <= end <= len(sequence):
        return None
    return sequence[:start] + edit.alt_bases + sequence[end:]


def _supplied_edits(transcript, known_variants, focal_variant):
    """Transcript edits of the supplied variants on this transcript's exons."""
    contig = getattr(transcript, "contig", focal_variant.contig)
    start, end = getattr(transcript, "start", 0), getattr(transcript, "end", float("inf"))
    edits = []
    for variant in known_variants.overlapping(contig, start, end):
        if variant == focal_variant:
            continue
        try:
            edit = transcript_edit_from_variant(variant, transcript)
        except ValueError:
            continue
        edits.append((edit, variant))
    return edits


def _is_substitution(edit):
    return 0 < len(edit.alt_bases) == edit.cdna_end - edit.cdna_start


def _attribute_edit(edit, flank, offset, supplied, known_variants):
    """
    Split one observed edit into (edit, supplied variant or None) pieces.

    An edit is a supplied variant when it gives the same flank sequence, which
    also matches an indel shifted within a repeat; the supplied representation
    is reported. A run of substitutions is split into the supplied
    substitutions it contains, and the rest stays unexplained. Candidates of
    conflicting origin leave the edit unexplained.
    """
    observed = _apply_edit(flank, offset, edit)
    matches = sorted(
        ((candidate, variant) for candidate, variant in supplied
         if _apply_edit(flank, offset, candidate) == observed),
        key=lambda pair: _source_variant_sort_key(pair[1]))
    if matches:
        if len({known_variants.origin(variant) for _, variant in matches}) > 1:
            return [(edit, None)]
        candidate, variant = matches[0]
        return [(TranscriptEdit(cdna_start=candidate.cdna_start, cdna_end=candidate.cdna_end,
                                alt_bases=candidate.alt_bases, source_variant=variant), variant)]
    if not _is_substitution(edit):
        return [(edit, None)]
    inside = [(candidate, variant) for candidate, variant in supplied
              if _is_substitution(candidate) and edit.cdna_start <= candidate.cdna_start
              and candidate.cdna_end <= edit.cdna_end
              and edit.alt_bases[candidate.cdna_start - edit.cdna_start:
                                 candidate.cdna_end - edit.cdna_start] == candidate.alt_bases]
    # Overlapping candidates would each explain the same bases; use neither.
    unique = [(candidate, variant) for candidate, variant in inside if not any(
        other is not candidate and other.cdna_start < candidate.cdna_end
        and candidate.cdna_start < other.cdna_end for other, _ in inside)]
    if not unique:
        return [(edit, None)]
    pieces, position = [], edit.cdna_start
    for candidate, variant in sorted(unique, key=lambda pair: pair[0].cdna_start):
        if position < candidate.cdna_start:
            pieces.append((_sub_edit(edit, position, candidate.cdna_start), None))
        pieces.append((TranscriptEdit(cdna_start=candidate.cdna_start, cdna_end=candidate.cdna_end,
                                      alt_bases=candidate.alt_bases, source_variant=variant), variant))
        position = candidate.cdna_end
    if position < edit.cdna_end:
        pieces.append((_sub_edit(edit, position, edit.cdna_end), None))
    return pieces


def _sub_edit(edit, start, end):
    return TranscriptEdit(
        cdna_start=start, cdna_end=end,
        alt_bases=edit.alt_bases[start - edit.cdna_start:end - edit.cdna_start], source_variant=None)


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


def categorize_transcript_assembly_edits_from_translation(translation, transcript, known_variants=None):
    """
    Group an assembled cDNA's differences from ``transcript`` by origin.

    Parameters
    ----------
    translation : Translation
    transcript : pyensembl.Transcript
    known_variants : KnownVariants or None
        Supplied somatic and germline variants. Without them only the
        translation's own variant is known.

    Returns
    -------
    dict
        ``known_somatic`` (the translation's variant, then co-somatic edits
        that are other supplied somatic variants), ``known_germline`` and
        ``unexplained``, each a tuple of TranscriptAssemblyEdit.
    """
    focal_variant = translation.reference_context.variant
    known_somatic = [transcript_edit_from_variant(focal_variant, transcript)]
    known_germline, unexplained = [], []
    observed = unexplained_transcript_edits_from_translation(translation, transcript)
    if known_variants is None or not observed:
        unexplained = list(observed)
    else:
        focal_edit, (prefix, suffix) = _flanks(translation, transcript)
        supplied = _supplied_edits(transcript, known_variants, focal_variant)
        for edit in observed:
            flank, offset = prefix if edit.cdna_end <= focal_edit.cdna_start else suffix
            for piece, variant in _attribute_edit(edit, flank, offset, supplied, known_variants):
                origin = None if variant is None else known_variants.origin(variant)
                {"somatic": known_somatic, "germline": known_germline, None: unexplained}[origin].append(piece)
    return {
        "known_somatic": tuple(_transcript_assembly_edit(transcript, edit) for edit in known_somatic),
        "known_germline": tuple(_transcript_assembly_edit(transcript, edit) for edit in known_germline),
        "unexplained": tuple(_transcript_assembly_edit(transcript, edit) for edit in unexplained),
    }
