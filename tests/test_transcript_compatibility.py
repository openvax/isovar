from types import SimpleNamespace

from isovar.allele_read import AlleleRead
from isovar.transcript_compatibility import (
    annotate_reads_with_transcript_compatibility,
)

from .common import eq_


def transcript(transcript_id, exon_intervals):
    return SimpleNamespace(id=transcript_id, exon_intervals=exon_intervals)


SHORT = transcript("short", [(101, 120), (201, 220)])
LONG = transcript("long", [(101, 150), (201, 220)])


def classified_ids(reference_blocks, splice_junctions=()):
    read = AlleleRead(
        prefix="AAAA",
        allele="C",
        suffix="TTTT",
        name="read",
        reference_blocks=reference_blocks,
        splice_junctions=splice_junctions)
    classified = annotate_reads_with_transcript_compatibility(
        [read], [SHORT, LONG])
    return classified[0].compatible_transcript_ids


def test_read_before_discriminating_junction_remains_ambiguous():
    # Base-0 half-open [110, 119) lies within both alternative exons.
    eq_(classified_ids(((0, 8, 110, 118),)), frozenset({"short", "long"}))


def test_junction_read_selects_compatible_splice_path():
    eq_(
        classified_ids(((0, 4, 116, 120), (4, 8, 200, 204)), ((120, 200),)),
        frozenset({"short"}))
    eq_(
        classified_ids(((0, 4, 146, 150), (4, 8, 200, 204)), ((150, 200),)),
        frozenset({"long"}))


def test_contiguous_read_crossing_alternative_donor_selects_long_exon():
    # The read passes the short transcript's donor without a reference skip.
    eq_(classified_ids(((0, 8, 116, 124),)), frozenset({"long"}))


def test_unaligned_read_metadata_does_not_choose_a_transcript():
    eq_(classified_ids(()), frozenset({"short", "long"}))


def test_alignment_beyond_terminal_transcript_boundary_is_compatible():
    terminal = transcript("terminal", [(101, 200)])
    read = AlleleRead(
        prefix="AAAA",
        allele="C",
        suffix="TTTT",
        name="read",
        reference_blocks=((0, 50, 75, 125),))

    annotated, = annotate_reads_with_transcript_compatibility(
        [read], [terminal])

    eq_(annotated.compatible_transcript_ids, frozenset({"terminal"}))


def test_deletion_gap_within_one_exon_is_not_a_splice_conflict():
    single_exon = transcript("deletion", [(101, 200)])
    read = AlleleRead(
        prefix="A" * 10,
        allele="",
        suffix="T" * 10,
        name="read",
        reference_blocks=((0, 10, 120, 130), (10, 20, 145, 155)))

    annotated, = annotate_reads_with_transcript_compatibility(
        [read], [single_exon])

    eq_(annotated.compatible_transcript_ids, frozenset({"deletion"}))
