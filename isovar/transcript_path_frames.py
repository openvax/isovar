"""Conditional reading frames on observed RNA paths, shared across event types.

Placements and RNA sequence are observations supplied by the caller. Reference
transcripts supply coordinate-anchored frame hypotheses, never missing bases.
This module does not acquire reads, establish molecular linkage or measure
translation. All query and reference intervals are zero-based and half-open.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from functools import cached_property
from types import MappingProxyType
from typing import Literal, NamedTuple

from .fusion import FusionReference
from .genetic_code import standard_genetic_code


Position = tuple[str, int, str]
PathPositions = Sequence[Position | None]
FocalCoverage = Literal[
    "not_assessed", "covered", "partially_covered", "upstream_of_frame_anchor",
    "after_stop", "beyond_observed_translation",
]


def next_position(position: Position) -> Position:
    """Return the next genomic base in this oriented path."""
    contig, base, strand = position
    return contig, base + (1 if strand == "+" else -1), strand


def reverse_position(position: Position) -> Position:
    """Reverse the strand of one placed base without moving it."""
    contig, base, strand = position
    return contig, base, "-" if strand == "+" else "+"


class ExactAnchor(NamedTuple):
    """A contiguous exact RNA/transcript match and its transcript offset."""

    query_start: int
    query_end: int
    transcript_start: int


@dataclass(frozen=True)
class TranscriptPathModel:
    """A versioned transcript with reusable immutable coordinate indexes.

    ``reference`` uses the existing FusionReference model shared with supplied
    ORFs. It currently validates standard genetic code 1 only. Positions follow
    transcript orientation; offsets map those positions back to cDNA indices.
    """

    reference: FusionReference

    def __post_init__(self) -> None:
        if not isinstance(self.reference, FusionReference):
            raise TypeError("A transcript path model requires a FusionReference")

    @cached_property
    def positions(self) -> tuple[Position, ...]:
        """Placed transcript bases, cached in transcript orientation."""
        reference = self.reference
        positions = [(reference.contig, base, reference.strand)
                     for start, end in reference.exons for base in range(start, end)]
        if reference.strand == "-":
            positions.reverse()
        return tuple(positions)

    @cached_property
    def offsets(self) -> Mapping[Position, int]:
        """Read-only genomic placement to transcript-offset lookup."""
        return MappingProxyType({position: index for index, position in enumerate(self.positions)})

    @cached_property
    def protein(self) -> str | None:
        """Supplied reference CDS translation, or None for an absent CDS."""
        reference = self.reference
        return None if reference.cds_start is None else standard_genetic_code.translate(
            reference.sequence[reference.cds_start:reference.cds_end], first_codon_is_start=True)[0]

    def __reduce__(self):
        # Rebuild derived indexes; MappingProxyType is intentionally read-only
        # and is not itself pickleable.
        return type(self), (self.reference,)

    def junctions(self) -> set[tuple[Position, Position]]:
        """Return annotated exon joins in both path orientations."""
        pairs = set()
        for left, right in zip(self.positions, self.positions[1:]):
            if right != next_position(left):
                pairs.update({(left, right), (reverse_position(right), reverse_position(left))})
        return pairs


def exact_anchor_runs(sequence: str, positions: PathPositions, model: TranscriptPathModel) -> list[ExactAnchor]:
    """Find maximal exact, collinear sequence matches to one transcript."""
    if len(sequence) != len(positions):
        raise ValueError("Observed RNA needs one placement per base")
    runs: list[ExactAnchor] = []
    for query, (base, position) in enumerate(zip(sequence, positions)):
        offset = model.offsets.get(position) if position is not None else None
        if offset is None or model.reference.sequence[offset] != base:
            continue
        if runs and runs[-1].query_end == query and runs[-1].transcript_start + query - runs[-1].query_start == offset:
            runs[-1] = runs[-1]._replace(query_end=query + 1)
        else:
            runs.append(ExactAnchor(query, query + 1, offset))
    return runs


def retained_intron_departures(positions: PathPositions, model: TranscriptPathModel) -> dict[int, str]:
    """Find observed genomic continuations across this model's exon ends.

    The model would splice to its next exon, but the observed RNA instead
    continues into the intervening intron on the same strand. This is a path
    interpretation relative to one model, not a mature-isoform/causation claim.
    Unplaced sequence, terminal extensions and genomic jumps do not qualify.
    """
    departures = {}
    for query in range(1, len(positions)):
        previous, current = positions[query - 1:query + 1]
        offset = model.offsets.get(previous) if previous is not None else None
        if previous is None or offset is None or offset + 1 >= len(model.positions) or current != next_position(previous):
            continue
        expected = model.positions[offset + 1]
        if expected != current:
            departures[query] = "retained_intron_entry"
    return departures


@dataclass(frozen=True)
class PathFrameReading:
    """One conditional coding frame over supplied, observed RNA sequence.

    ``translation_end`` includes an observed stop codon when present; amino
    acids exclude it. ``complete_5prime`` means the observed anchor includes
    the supplied CDS start, not that initiation was measured. ``departure``
    is the first positional departure from that model. Substitutions keep the
    frame. ``focal_coverage`` concerns only a caller-supplied allele interval.
    """

    translation_start: int
    translation_end: int
    amino_acids: str
    ends_with_stop_codon: bool
    complete_5prime: bool
    reference_start: int
    departure: int | None
    departure_relation: str | None
    substitutions: tuple[int, ...]
    anchors: tuple[ExactAnchor, ...]
    focal_coverage: FocalCoverage = "not_assessed"


@dataclass(frozen=True)
class TranscriptFrameAssessment:
    """Frame outcomes for one model; ordinary and departing paths stay distinct."""

    model: TranscriptPathModel
    anchors: tuple[ExactAnchor, ...]
    readings: tuple[PathFrameReading, ...] = ()
    collinear_readings: tuple[PathFrameReading, ...] = ()
    other_departures: int = 0
    reason: Literal["no_sequence_anchor", "CDS_unavailable", "no_coding_anchor"] | None = None


def _focal_coverage(reading: PathFrameReading, interval: tuple[int, int] | None) -> FocalCoverage:
    if interval is None:
        return "not_assessed"
    start, end = interval
    if start == end:
        # A deletion boundary requires translated RNA on both sides.
        covered = reading.translation_start < start < reading.translation_end
    else:
        covered = reading.translation_start <= start < end <= reading.translation_end
    if covered:
        return "covered"
    if end <= reading.translation_start:
        return "upstream_of_frame_anchor"
    if start >= reading.translation_end:
        return "after_stop" if reading.ends_with_stop_codon else "beyond_observed_translation"
    return "partially_covered"


def _reading_from_anchor(sequence: str, positions: PathPositions, model: TranscriptPathModel,
                         anchor: ExactAnchor, minimum: int) -> PathFrameReading | None:
    reference = model.reference
    assert reference.cds_start is not None and reference.cds_end is not None
    query_start, query_end, transcript_start = anchor
    low = max(query_start, query_start + reference.cds_start - transcript_start)
    high = min(query_end, query_start + reference.cds_end - 3 - transcript_start)
    if high - low < minimum:
        return None
    start = low + (reference.cds_start - (transcript_start + low - query_start)) % 3
    offset = transcript_start + start - query_start
    protein, stop = standard_genetic_code.translate(sequence[start:], first_codon_is_start=offset == reference.cds_start)
    end = start + 3 * (len(protein) + stop)
    departure, substitutions = None, []
    for query in range(start, end):
        expected = offset + query - start
        if expected >= len(model.positions) or positions[query] not in (None, model.positions[expected]):
            departure = query
        elif sequence[query] != reference.sequence[expected]:
            if positions[query] is None:
                departure = query
            else:
                substitutions.append(query)
        if departure is not None:
            break
    return PathFrameReading(start, end, protein, bool(stop), offset == reference.cds_start,
                            offset, departure, None, tuple(substitutions), (anchor,))


def _retain_longest(readings: dict[int, PathFrameReading], reading: PathFrameReading) -> None:
    """Reads in one frame ending at the same codon share their longest context."""
    previous = readings.get(reading.translation_end)
    if previous is None:
        readings[reading.translation_end] = reading
        return
    longest = reading if reading.translation_start < previous.translation_start else previous
    readings[reading.translation_end] = replace(longest, anchors=previous.anchors + reading.anchors)


def assess_transcript_path(
        sequence: str, positions: PathPositions, model: TranscriptPathModel, *,
        min_anchor_bases: int = 12, permitted_departures: Mapping[int, str] | None = None,
        allow_retained_introns: bool = False,
        focal_interval: tuple[int, int] | None = None) -> TranscriptFrameAssessment:
    """Transfer a supplied CDS frame downstream over an observed RNA path.

    Parameters
    ----------
    sequence : str
        Uppercase observed RNA represented with DNA bases. No reference fill.
    positions : sequence of (contig, position, strand) or None
        One placement per nucleotide, in the path's 5-prime to 3-prime order.
        None means unplaced, not proof of a genomic or transcript match.
    model : TranscriptPathModel
        Explicit versioned transcript and reusable coordinate index.
    min_anchor_bases : int
        Minimum exact collinear coding match, at least one codon.
    permitted_departures : mapping of int to str, optional
        Query offsets whose observed departure has been assessed by the caller,
        such as an SV junction or a literal indel. Other departures are counted
        separately and cannot produce a supported reading.
    allow_retained_introns : bool
        Also permit observed exon-to-intron continuations relative to this model.
    focal_interval : pair of int, optional
        An already established allele's query interval. A deletion is an empty
        boundary. Each reading reports whether the allele lies in translation,
        after a stop, or outside the observed frame. This does not establish
        allele causation, phase, or a changed protein by itself.
    """
    if len(sequence) != len(positions) or not sequence or set(sequence) - set("ACGT"):
        raise ValueError("Observed RNA needs explicit uppercase A/C/G/T bases and one placement per base")
    if any(position is not None and (not isinstance(position, tuple) or len(position) != 3 or
            not isinstance(position[0], str) or not position[0] or type(position[1]) is not int or
            position[1] < 0 or position[2] not in ("+", "-")) for position in positions):
        raise ValueError("Each placement needs a contig, nonnegative position and +/- strand")
    if type(min_anchor_bases) is not int or min_anchor_bases < 3:
        raise ValueError("A frame anchor needs at least three bases")
    if type(allow_retained_introns) is not bool:
        raise ValueError("allow_retained_introns must be a boolean")
    if focal_interval is not None and (len(focal_interval) != 2 or
            any(type(value) is not int for value in focal_interval) or
            not 0 <= focal_interval[0] <= focal_interval[1] <= len(sequence)):
        raise ValueError("The focal interval must lie within the observed RNA")
    departures = dict(permitted_departures or {})
    if any(type(offset) is not int or not 0 <= offset < len(sequence) or not isinstance(relation, str) or not relation
           for offset, relation in departures.items()):
        raise ValueError("Permitted departures need valid query offsets and nonempty relations")
    if allow_retained_introns:
        for offset, relation in retained_intron_departures(positions, model).items():
            departures.setdefault(offset, relation)
    anchors = tuple(anchor for anchor in exact_anchor_runs(sequence, positions, model)
                    if anchor.query_end - anchor.query_start >= min_anchor_bases)
    if not anchors:
        return TranscriptFrameAssessment(model, anchors, reason="no_sequence_anchor")
    if model.reference.cds_start is None:
        return TranscriptFrameAssessment(model, anchors, reason="CDS_unavailable")
    readings: dict[int, PathFrameReading] = {}
    collinear: dict[int, PathFrameReading] = {}
    coding_anchors = other_departures = 0
    for anchor in anchors:
        reading = _reading_from_anchor(sequence, positions, model, anchor, min_anchor_bases)
        if reading is None:
            continue
        coding_anchors += 1
        reading = replace(reading, focal_coverage=_focal_coverage(reading, focal_interval))
        if reading.departure is None:
            _retain_longest(collinear, reading)
        elif reading.departure in departures:
            reading = replace(reading, departure_relation=departures[reading.departure])
            _retain_longest(readings, reading)
        else:
            other_departures += 1
    return TranscriptFrameAssessment(
        model, anchors, tuple(sorted(readings.values(), key=lambda reading: reading.translation_start)),
        tuple(sorted(collinear.values(), key=lambda reading: reading.translation_start)), other_departures,
        None if coding_anchors else "no_coding_anchor")


def _legacy_readings(assessment: TranscriptFrameAssessment) -> dict | None:
    """Adapt the shared typed result to the existing SV grouping contract."""
    if assessment.reason == "no_sequence_anchor":
        return None
    anchors = [list(anchor) for anchor in assessment.anchors]
    if assessment.reason:
        return dict(reason=assessment.reason, anchors=anchors)
    return dict(readings=[dict(
        translation_start=reading.translation_start, translation_end=reading.translation_end,
        amino_acids=reading.amino_acids, ends_with_stop_codon=reading.ends_with_stop_codon,
        complete_5prime=reading.complete_5prime, departure=reading.departure,
        substitutions=list(reading.substitutions), anchors=[list(anchor) for anchor in reading.anchors])
        for reading in assessment.readings], reference_reading=bool(assessment.collinear_readings),
        other_departures=assessment.other_departures, anchors=anchors)


def sv_frame_evidence(sequence: str, positions: PathPositions, departures: Mapping[int, str],
                      models: Sequence[TranscriptPathModel], min_anchor: int,
                      reference_peptides: set[str], lengths: Sequence[int]) -> dict:
    """Group translations across all anchored models; keep every alternative."""
    groups, unresolved, reference_readings, other_departures = {}, [], [], 0
    for model in models:
        result = _legacy_readings(assess_transcript_path(
            sequence, positions, model, min_anchor_bases=min_anchor, permitted_departures=departures))
        if result is None:
            continue
        reference = model.reference
        if "reason" in result:
            unresolved.append(dict(transcript_id=reference.transcript_id, annotation=reference.annotation,
                                   reason=result["reason"], anchors=result["anchors"]))
            continue
        if result["reference_reading"]:
            reference_readings.append(reference.transcript_id)
        other_departures += result["other_departures"]
        for reading in result["readings"]:
            key = (reading["translation_start"], reading["amino_acids"], reading["ends_with_stop_codon"])
            if key not in groups:
                protein = reading["amino_acids"]
                groups[key] = dict(
                    translation_start=reading["translation_start"], translation_end=reading["translation_end"],
                    amino_acids=protein, ends_with_stop_codon=reading["ends_with_stop_codon"],
                    transcript_ids=[], frame_evidence=[], competing_noncoding_models=[],
                    candidate_peptides=[dict(sequence=protein[i:i + k], protein_interval=[i, i + k])
                                        for k in lengths for i in range(len(protein) - k + 1)
                                        if protein[i:i + k] not in reference_peptides],
                    peptide_novelty="absent_from_supplied_reference_proteins", translation_observed=False)
            groups[key]["transcript_ids"].append(reference.transcript_id)
            groups[key]["frame_evidence"].append(dict(
                transcript_id=reference.transcript_id, annotation=reference.annotation,
                basis="exact_collinear_CDS_anchor", anchors=reading["anchors"],
                departure=reading["departure"], departure_relation=departures[reading["departure"]],
                frame_preserving_substitutions=reading["substitutions"],
                complete_5prime=reading["complete_5prime"],
                upstream_frame_assumed=not reading["complete_5prime"]))
    translations = sorted(groups.values(), key=lambda g: (g["translation_start"], g["amino_acids"]))
    for translation in translations:
        translation["departure_relations"] = sorted({e["departure_relation"] for e in translation["frame_evidence"]})
        translation["complete_5prime"] = any(e["complete_5prime"] for e in translation["frame_evidence"])
        # A noncoding model matching the same frame-anchor sequence is an
        # alternative explanation, not permission to drop the protein.
        spans = [a for e in translation["frame_evidence"] for a in e["anchors"]]
        translation["competing_noncoding_models"] = sorted(
            u["transcript_id"] for u in unresolved
            if any(a[0] < s[1] and s[0] < a[1] for a in u["anchors"] for s in spans))
    if not translations:
        # An exact coding anchor whose every reading leaves the model at an
        # indel or error is still a gene anchor, not the absence of one.
        status = ("reference_protein_only" if reference_readings else "unresolved" if unresolved
                  else "departs_elsewhere_only" if other_departures else "no_gene_anchor")
    elif len(translations) == 1 and not translations[0]["competing_noncoding_models"]:
        status = "translated"
    else:
        status = "ambiguous"
    return dict(frame_status=status, translations=translations, unresolved_models=unresolved,
                reference_readings=sorted(reference_readings), readings_departing_elsewhere=other_departures)

