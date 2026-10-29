"""Observed allele-centred RNA, kept independently of a reference coding frame."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterable, Literal

from .allele_read import AlleleRead
from .rna_evidence import EvidenceSets, content_identifier
from .transcript_compatibility import annotate_reads_with_transcript_compatibility
from .variant_helpers import require_literal_variant, trim_variant
from .variant_sequence import VariantSequence

if TYPE_CHECKING:
    from pyensembl import Transcript
    from varcode import Variant
    from .protein_sequence_creator import ProteinSequenceCreator
    from .reference_context import ReferenceContext
    from .translation import Translation


AssemblyMode = Literal["configured", "read_contexts"]
AssemblyState = Literal["assembled", "no_alt_reads", "no_retained_rna"]
FrameState = Literal[
    "rna_frame_candidates", "no_alt_reads", "no_coding_reference_context",
    "no_reference_path_compatible_reads", "no_retained_compatible_rna",
    "no_retained_rna_frame_match",
]


@dataclass(frozen=True)
class VariantRNAAssembly:
    """RNA candidates and classified observations for one supplied model set.

    ``compatible_sequences`` obey at least one supplied transcript path;
    ``unmatched_sequences`` retain RNA that obeys none. Neither collection is a
    translation claim. Their original observations may overlap between paths;
    support counts must be combined by identity, never added across candidates.
    ``mode`` records whether the configured assembler or bounded direct read
    contexts were used. The latter does not extend reads by overlap assembly.
    """

    compatible_sequences: tuple[VariantSequence, ...]
    unmatched_sequences: tuple[VariantSequence, ...]
    compatible_reads: tuple[AlleleRead, ...]
    unmatched_reads: tuple[AlleleRead, ...]
    mode: AssemblyMode

    @property
    def state(self) -> AssemblyState:
        if self.compatible_sequences or self.unmatched_sequences:
            return "assembled"
        return "no_retained_rna" if self.compatible_reads or self.unmatched_reads else "no_alt_reads"

    def candidate_records(self, evidence: EvidenceSets, labels=None) -> list[tuple[VariantSequence, dict]]:
        """Return observed sequences and JSON records with scoped contributors.

        Intervals are zero-based, half-open on the genomic-forward RNA string.
        Contributors do not necessarily span the entire candidate sequence.
        Placement-sensitive IDs preserve alternatives from one original read.
        """
        rows = []
        for state, sequences in (("reference_path_compatible", self.compatible_sequences),
                                 ("no_supplied_reference_path", self.unmatched_sequences)):
            for sequence in sequences:
                support = evidence.support(sequence.reads, labels=labels)
                placements = sorted(content_identifier("placement", [read.source_alignments,
                                                                       read.source_alignment_paths])
                                    for read in sequence.reads)
                transcript_ids = sorted(sequence.compatible_transcript_ids or ())
                ident = content_identifier("rna_candidate", [state, sequence.sequence, len(sequence.prefix),
                                                              transcript_ids, support["evidence_set_id"], placements])
                rows.append((sequence, dict(
                    rna_candidate_id=ident, state=state, sequence=sequence.sequence,
                    variant_interval=[len(sequence.prefix), len(sequence.prefix) + len(sequence.alt)],
                    compatible_transcript_ids=transcript_ids, rna_support=support)))
        return rows


@dataclass(frozen=True)
class VariantReadAssessment:
    """A local RNA assembly and its ordinary-CDS translation hypotheses.

    ``translations`` contains all retained frame matches before protein ranking
    or the protein-count cap. A frame match does not establish translation in
    the sample. Missing frames never erase assembled RNA. Objects referenced by
    this frozen container retain their existing Isovar types and identities.
    """

    assembly: VariantRNAAssembly
    reference_contexts: tuple[ReferenceContext, ...]
    translations: tuple[Translation, ...]

    @property
    def frame_state(self) -> FrameState:
        if self.translations:
            return "rna_frame_candidates"
        if self.assembly.state == "no_alt_reads":
            return "no_alt_reads"
        if not self.reference_contexts:
            return "no_coding_reference_context"
        has_paths = any(context.transcripts for context in self.reference_contexts)
        if has_paths and not self.assembly.compatible_reads:
            return "no_reference_path_compatible_reads"
        if not (self.assembly.compatible_sequences if has_paths else self.assembly.unmatched_sequences):
            return "no_retained_compatible_rna"
        return "no_retained_rna_frame_match"

    def to_dict(self, evidence: EvidenceSets, labels=None) -> dict:
        """Export RNA and per-path frame outcomes into an existing evidence scope.

        Reference transcript IDs identify frame hypotheses, not an observed
        start codon. Protein sequences and their interval support remain in the
        protein-hypotheses export. Unassessed legacy results use ``None`` instead
        of manufacturing an empty assessment.
        """
        framed: dict[int, set[str]] = {}
        for translation in self.translations:
            framed.setdefault(id(translation.untrimmed_variant_sequence), set()).update(
                transcript.id for transcript in translation.reference_context.transcripts)
        rows = []
        for sequence, row in self.assembly.candidate_records(evidence, labels):
            has_frame = id(sequence) in framed
            transcript_ids = sorted(framed.get(id(sequence), ()))
            row.update(
                frame_state=("rna_frame_candidates" if has_frame else
                             "no_coding_reference_context" if not self.reference_contexts else
                             "no_reference_path_compatible_reads" if not sequence.compatible_transcript_ids
                             and any(context.transcripts for context in self.reference_contexts) else
                             "no_retained_rna_frame_match"),
                frame_transcript_ids=transcript_ids)
            rows.append(row)
        return dict(
            schema="isovar.variant_read_assessment.v1",
            assembly_state=self.assembly.state, assembly_mode=self.assembly.mode,
            frame_state=self.frame_state, rna_candidates=rows,
            reference_context_count=len(self.reference_contexts),
            translation_count_before_protein_selection=len(self.translations),
            reference_incompatible_support=evidence.support(self.assembly.unmatched_reads, labels=labels),
            limitations=["Local allele-centred RNA only; no whole-transcript discovery claim.",
                         "RNA contributors need not span the entire candidate sequence.",
                         "Reference frame assignment is not observed initiation or translation."])


def assemble_variant_reads(
        variant: Variant, reads: Iterable[AlleleRead], transcripts: Iterable[Transcript], *,
        creator: ProteinSequenceCreator, mode: AssemblyMode = "configured") -> VariantRNAAssembly:
    """Assemble coding-compatible and unmatched ALT RNA using one policy.

    No supplied transcript is required. Each group uses the same assembler,
    coverage floor, context budget and original-placement compatibility rules.
    Reads must carry the nominated trimmed ALT allele. Reference-path labels
    are recomputed for this model set; unmatched reads lose those labels only
    for assembly and are never passed to a supplied coding frame.
    """
    require_literal_variant(variant)
    if mode not in ("configured", "read_contexts"):
        raise ValueError("assembly mode must be 'configured' or 'read_contexts'")
    reads = list(reads)
    alt = trim_variant(variant)[2]
    if any(read.allele != alt for read in reads):
        raise ValueError("Variant assembly requires reads carrying the nominated ALT allele")
    if reads:
        left, right = creator._variant_sequence_creator.flanking_sequence_lengths(len(alt))
        reads = annotate_reads_with_transcript_compatibility(reads, transcripts, left, right)
    compatible = tuple(read for read in reads if read.compatible_transcript_ids)
    unmatched = tuple(read for read in reads if not read.compatible_transcript_ids)
    return VariantRNAAssembly(
        compatible_sequences=tuple(creator.variant_sequences_from_reads(variant, compatible, mode=mode)),
        unmatched_sequences=tuple(creator.variant_sequences_from_reads(
            variant, (read.with_compatible_transcript_ids(None) for read in unmatched), mode=mode)),
        compatible_reads=compatible, unmatched_reads=unmatched, mode=mode)
