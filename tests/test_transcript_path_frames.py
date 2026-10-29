"""Observed intronic paths must retain frame/allele uncertainty on both strands."""

from dataclasses import FrozenInstanceError, asdict, replace
import pickle

import pytest

from isovar.fusion import FusionReference
from isovar.transcript_path_frames import (
    TranscriptPathModel, assess_transcript_path, exact_anchor_runs,
    next_position, retained_intron_departures,
)


def model(strand="+"):
    return TranscriptPathModel(FusionReference(
        "T", "fixture", "models-v1", "chr1", strand, ((100, 112), (200, 218)),
        "ATG" + "GCC" * 8 + "TAA", 0, 30))


def retention(strand="+"):
    reference = model(strand)
    boundary = next(i for i in range(1, len(reference.positions))
                    if reference.positions[i] != next_position(reference.positions[i - 1]))
    sequence = reference.reference.sequence[:boundary] + "GAA" * 3 + "TAA" + "ACC" * 4
    start = reference.positions[0][1]
    step = 1 if strand == "+" else -1
    positions = tuple(("chr1", start + step * i, strand) for i in range(len(sequence)))
    return reference, sequence, positions, boundary


@pytest.mark.parametrize("strand", ["+", "-"])
def test_retention_frame_and_stop_before_focal_allele(strand):
    reference, sequence, positions, boundary = retention(strand)
    assert retained_intron_departures(positions, reference) == {boundary: "retained_intron_entry"}
    disallowed = assess_transcript_path(sequence, positions, reference, min_anchor_bases=6)
    assert not disallowed.readings and disallowed.other_departures == 1
    assessed = assess_transcript_path(sequence, positions, reference, min_anchor_bases=6,
                                     allow_retained_introns=True, focal_interval=(boundary + 15, boundary + 16))
    reading, = assessed.readings
    assert reading.departure == boundary and reading.departure_relation == "retained_intron_entry"
    assert reading.amino_acids == "M" + "A" * (boundary // 3 - 1) + "EEE"
    assert reading.translation_end == boundary + 12 and reading.ends_with_stop_codon
    assert reading.complete_5prime and reading.focal_coverage == "after_stop"
    # Same path, but a genuinely translated focal interval.
    covered = assess_transcript_path(sequence, positions, reference, min_anchor_bases=6,
                                    allow_retained_introns=True, focal_interval=(boundary + 1, boundary + 2))
    assert covered.readings[0].focal_coverage == "covered"


@pytest.mark.parametrize("strand", ["+", "-"])
def test_unplaced_and_terminal_extensions_are_not_retention(strand):
    reference, sequence, positions, boundary = retention(strand)
    unplaced = positions[:boundary] + (None,) + positions[boundary + 1:]
    assert retained_intron_departures(unplaced, reference) == {}
    terminal = reference.positions + (next_position(reference.positions[-1]),)
    assert retained_intron_departures(terminal, reference) == {}


@pytest.mark.parametrize("strand", ["+", "-"])
def test_cryptic_jump_requires_explicit_departure_and_preserves_partial_frame(strand):
    reference, sequence, positions, boundary = retention(strand)
    # A positional jump is not evidence of contiguous intron retention.
    shifted = positions[:boundary] + tuple((c, p + 30, s) for c, p, s in positions[boundary:])
    result = assess_transcript_path(sequence[3:], shifted[3:], reference, min_anchor_bases=6,
                                   permitted_departures={boundary - 3: "observed_cryptic_junction"},
                                   focal_interval=(boundary, boundary + 1))
    reading, = result.readings
    assert not reading.complete_5prime and reading.reference_start == 3
    assert reading.departure_relation == "observed_cryptic_junction" and reading.focal_coverage == "covered"
    blocked = assess_transcript_path(sequence, shifted, reference, min_anchor_bases=6, allow_retained_introns=True)
    assert not blocked.readings and blocked.other_departures == 1


def test_reference_collinear_substitution_remains_a_conditional_frame():
    reference = model()
    sequence = reference.reference.sequence[:15] + "T" + reference.reference.sequence[16:]
    result = assess_transcript_path(sequence, reference.positions, reference,
                                   min_anchor_bases=6, focal_interval=(15, 16))
    assert not result.readings
    reading, = result.collinear_readings
    assert reading.focal_coverage == "covered" and 15 in reading.substitutions
    assert reading.amino_acids != reference.protein
    assert reading.complete_5prime and reading.ends_with_stop_codon


@pytest.mark.parametrize("interval,expected", [
    ((10, 10), "covered"), ((0, 0), "upstream_of_frame_anchor"),
    ((24, 24), "after_stop"), ((23, 25), "partially_covered"),
])
def test_focal_interval_distinguishes_deletion_boundaries_and_partial_coverage(interval, expected):
    reference, sequence, positions, _ = retention()
    result = assess_transcript_path(sequence, positions, reference, min_anchor_bases=6,
                                   allow_retained_introns=True, focal_interval=interval)
    assert result.readings[0].focal_coverage == expected


def test_unfinished_codon_does_not_establish_allele_translation():
    reference, sequence, positions, boundary = retention()
    end = boundary + 7
    result = assess_transcript_path(sequence[:end], positions[:end], reference, min_anchor_bases=6,
                                   allow_retained_introns=True, focal_interval=(end - 1, end))
    reading, = result.readings
    assert not reading.ends_with_stop_codon and reading.focal_coverage == "beyond_observed_translation"


def test_noncoding_and_missing_anchor_outcomes_and_model_immutability():
    reference = model()
    noncoding = TranscriptPathModel(replace(reference.reference, cds_start=None, cds_end=None))
    result = assess_transcript_path(reference.reference.sequence, reference.positions, noncoding)
    assert result.reason == "CDS_unavailable" and not result.readings
    assert assess_transcript_path("T" * 30, reference.positions, reference).reason == "no_sequence_anchor"
    with pytest.raises(TypeError):
        reference.offsets[reference.positions[0]] = 99
    with pytest.raises(FrozenInstanceError):
        reference.protein = "fake"
    with pytest.raises(ValueError, match="placement"):
        exact_anchor_runs("AAA", (), reference)
    copied = pickle.loads(pickle.dumps(reference))
    assert copied.positions == reference.positions and copied.offsets == reference.offsets
    assert asdict(result)["model"] == {"reference": asdict(noncoding.reference)}


def test_ambiguous_bases_and_invalid_placements_fail_before_translation():
    reference = model()
    sequence = reference.reference.sequence
    with pytest.raises(ValueError, match="explicit"):
        assess_transcript_path(sequence[:15] + "N" + sequence[16:], reference.positions, reference)
    for position in [("chr1", True, "+"), ("chr1", -1, "+"), ("chr1", 1, "?"), ("", 1, "+")]:
        with pytest.raises(ValueError, match="placement"):
            assess_transcript_path(sequence, (position,) + reference.positions[1:], reference)


def test_utr_anchor_does_not_manufacture_a_coding_frame():
    reference = TranscriptPathModel(FusionReference(
        "UTR", "fixture", "models-v1", "chr1", "+", ((100, 139),),
        "CCC" * 3 + "ATG" + "GCC" * 8 + "TAA", 9, 39))
    result = assess_transcript_path(reference.reference.sequence[:9], reference.positions[:9], reference,
                                   min_anchor_bases=6)
    assert result.reason == "no_coding_anchor" and not result.readings and not result.collinear_readings


@pytest.mark.parametrize("changes", [
    {"min_anchor_bases": True}, {"min_anchor_bases": 2}, {"allow_retained_introns": 1},
    {"focal_interval": (True, 4)}, {"focal_interval": (4, 2)}, {"focal_interval": (0, 99)},
    {"permitted_departures": {True: "join"}}, {"permitted_departures": {5: ""}},
])
def test_invalid_frame_policy_is_explicit(changes):
    reference = model()
    with pytest.raises(ValueError):
        assess_transcript_path(reference.reference.sequence, reference.positions, reference, **changes)
