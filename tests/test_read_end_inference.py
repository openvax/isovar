"""End inference is additive; derived trimming must preserve alignment semantics."""

from dataclasses import asdict
from itertools import product
import json
from pathlib import Path

import pysam
import pytest

from isovar import (
    Adapter, AlleleRead, ReadCollector, ReadEndProfile, infer_read_ends,
    read_sequence_view_from_alignment,
)
from isovar.dna import reverse_complement_dna as reverse_complement
from isovar.read_end_inference import adapter_matches, read_end_profiles_from_json
from isovar.cli.rna_args import read_collector_from_args
from isovar.cli.main_args import make_isovar_arg_parser
from .mock_objects import MockAlignmentFile, make_pysam_read


MOTIF = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
BODY = "GCTACGTCGCTG"
PROFILE = ReadEndProfile("test-library", adapters=(Adapter("adapter", MOTIF),))


@pytest.mark.parametrize("observed", [
    MOTIF, MOTIF[:16], MOTIF[:13],
    MOTIF[:8] + "T" + MOTIF[9:],  # substitution
    MOTIF[:8] + "C" + MOTIF[8:],  # insertion
    MOTIF[:8] + MOTIF[9:],        # deletion
    MOTIF[:5] + "T" + MOTIF[6:12] + MOTIF[13:22] + "AC" + MOTIF[22:],  # mixed errors
])
@pytest.mark.parametrize("reverse", [False, True])
def test_noisy_and_partial_adapter_coordinate_mapping(observed, reverse):
    sequence = BODY + observed
    if reverse:
        sequence = reverse_complement(sequence)
    view = infer_read_ends(sequence, profile=PROFILE, reverse=reverse, trim_adapters=True)
    assert view.sequence == (reverse_complement(BODY) if reverse else BODY)
    assert view.original_sequence == sequence
    assert view.sequenced_interval(0, len(BODY)) == (0, len(BODY))
    assert view.annotations[0].kind == "adapter"
    assert view.original_qualities is None


def test_left_adapter_partial_suffix_and_custom_orientation():
    profile = ReadEndProfile("custom", adapters=(Adapter("left", MOTIF, end="left"),))
    assert infer_read_ends(MOTIF[-16:] + BODY, profile=profile, trim_adapters=True).sequence == BODY
    bidirectional = ReadEndProfile("cdna", adapters=PROFILE.adapters, both_orientations=True)
    assert infer_read_ends(reverse_complement(BODY + MOTIF), profile=bidirectional,
                           trim_adapters=True).sequence == reverse_complement(BODY)


def test_short_chance_matches_and_unconfigured_adapter_are_not_trimmed():
    for sequence, profile in ((BODY + MOTIF[:8], PROFILE), (BODY + MOTIF, None)):
        view = infer_read_ends(sequence, profile=profile, trim_adapters=True)
        assert view.sequence == sequence
        assert not view.annotations


@pytest.mark.parametrize("tail", ["A" * 26, "A" * 14 + "C" + "A" * 14])
def test_poly_a_is_heuristic_and_opt_in(tail):
    sequence = BODY + tail
    annotated = infer_read_ends(sequence)
    assert annotated.sequence == sequence
    assert annotated.annotations[0].kind == "poly_a"
    trimmed = infer_read_ends(sequence, trim_poly_a=True)
    assert trimmed.sequence == BODY
    reverse = infer_read_ends(reverse_complement(sequence), trim_poly_a=True)
    assert reverse.sequence == reverse_complement(BODY)


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("leading_tail", ["", "A" * 20])
def test_long_tail_cannot_rescue_a_scan_through_the_read_body(reverse, leading_tail):
    sequence = leading_tail + BODY + "A" * 100
    expected = BODY
    if reverse:
        sequence, expected = reverse_complement(sequence), reverse_complement(expected)
    view = infer_read_ends(sequence, trim_poly_a=True)
    assert view.sequence == expected


@pytest.mark.parametrize("reverse", [False, True])
def test_long_tail_preserves_genuine_opposite_soft_clip(reverse):
    sequence = BODY + "A" * 100
    cigar = "6S6M100S"
    if reverse:
        sequence, cigar = reverse_complement(sequence), "100S6M6S"
    read = make_pysam_read(sequence, cigar)
    read.flag = 16 if reverse else 0
    view = read_sequence_view_from_alignment(read, trim_poly_a=True)
    assert view.sequence == (reverse_complement(BODY) if reverse else BODY)


@pytest.mark.parametrize("tail", ["A" * 11 + "C", "A" * 5 + "C" + "A" * 18])
def test_local_tail_support_tolerates_terminal_and_internal_errors(tail):
    assert infer_read_ends(BODY + tail, trim_poly_a=True).sequence == BODY


def test_combined_adapter_and_tail_never_stitches_across_internal_sequence():
    sequence = BODY + "A" * 25 + MOTIF
    assert infer_read_ends(sequence, profile=PROFILE, trim_poly_a=True).sequence == sequence
    assert infer_read_ends(sequence, profile=PROFILE, trim_adapters=True).sequence == BODY + "A" * 25
    view = infer_read_ends(sequence, profile=PROFILE, trim_adapters=True, trim_poly_a=True)
    assert view.sequence == BODY
    assert {a.kind for a in view.annotations} == {"adapter", "poly_a"}


def test_window_saturated_tail_is_annotated_not_partially_trimmed():
    profile = ReadEndProfile(end_window=20)
    view = infer_read_ends(BODY + "A" * 40, profile=profile, trim_poly_a=True)
    assert view.sequence == BODY + "A" * 40
    assert view.annotations[-1].ambiguous


def test_equally_supported_adapter_boundaries_are_not_forced():
    motif = "ACGCTGCATCGT"
    sequence = BODY + motif + "GCGC" + motif
    profile = ReadEndProfile("ambiguous", adapters=(Adapter("repeat", motif),))
    view = infer_read_ends(sequence, profile=profile, trim_adapters=True)
    assert view.sequence == sequence
    assert len(view.annotations) == 2
    assert all(a.ambiguous for a in view.annotations)


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("aligned", [False, True])
def test_tied_adapter_starts_at_one_end_preserve_flanking_bases(reverse, aligned):
    sequence = BODY + MOTIF[1:]
    if reverse:
        sequence = reverse_complement(sequence)
    if aligned:
        cigar = f"6M{len(sequence) - 6}S" if not reverse else f"{len(sequence) - 6}S6M"
        read = make_pysam_read(sequence, cigar)
        read.flag = 16 if reverse else 0
        view = read_sequence_view_from_alignment(read, PROFILE, trim_adapters=True)
    else:
        view = infer_read_ends(sequence, profile=PROFILE, reverse=reverse, trim_adapters=True)
    assert view.sequence == sequence
    adapters = [a for a in view.annotations if a.kind == "adapter"]
    assert len(adapters) == 2
    assert all(a.ambiguous and a.errors == 1 for a in adapters)
    assert {a.end if reverse else a.start for a in adapters} == (
        {len(sequence) - 11, len(sequence) - 12} if reverse else {11, 12})


@pytest.mark.parametrize("side", ["left", "right"])
def test_adapter_locations_match_exhaustive_edit_distance(side):
    # Independent, deliberately small oracle checks every substring boundary,
    # including tied starts, tied ends and end-window coordinate offsets.
    def distance(first, second):
        row = list(range(len(second) + 1))
        for i, a in enumerate(first, 1):
            next_row = [i]
            for j, b in enumerate(second, 1):
                next_row.append(min(row[j] + 1, next_row[-1] + 1, row[j - 1] + (a != b)))
            row = next_row
        return row[-1]

    adapter = Adapter("short", "ACGT", side, min_overlap=4, max_error_rate=0.3)
    for letters in product("ACGT", repeat=5):
        sequence = "GC" + "".join(letters) + "CT"
        first, last = (0, 6) if side == "left" else (len(sequence) - 6, len(sequence))
        scores = {(start, end): distance(adapter.sequence, sequence[start:end])
                  for start in range(first, last) for end in range(start + 1, last + 1)}
        best = min(scores.values())
        expected = {span for span, score in scores.items() if score == best and score <= 1}
        observed = {(hit.start, hit.end) for hit in adapter_matches(sequence, adapter, window=6)}
        assert observed == expected, (sequence, side)


def test_adapter_conflict_cannot_enable_tail_trimming_across_it():
    sequence = BODY + "A" * 20 + MOTIF
    read = make_pysam_read(sequence, "38M27S")
    view = read_sequence_view_from_alignment(read, PROFILE, trim_adapters=True, trim_poly_a=True)
    assert view.sequence == sequence


def test_empty_sequences_and_quality_validation():
    assert infer_read_ends("", profile=PROFILE, trim_adapters=True, trim_poly_a=True).sequence == ""
    with pytest.raises(ValueError, match="quality"):
        infer_read_ends("ACGT", qualities=[30])


def test_available_qualities_slice_in_lockstep_and_remain_immutable():
    qualities = list(range(len(BODY) + 20))
    view = infer_read_ends(BODY + "A" * 20, qualities, trim_poly_a=True)
    assert view.qualities == tuple(qualities[:len(BODY)])
    qualities[0] = 99
    assert view.original_qualities[0] == 0


@pytest.mark.parametrize("cigar", ["12M26S", "12M26I", "38M", "12M12I14S"])
def test_tail_trimming_preserves_mapped_and_inserted_bases(cigar):
    read = make_pysam_read(BODY + "A" * 26, cigar)
    before = read.to_dict() if read.header is not None else (read.seq, read.cigarstring, read.qual)
    view = read_sequence_view_from_alignment(read, trim_poly_a=True)
    expected_end = max(len(BODY), read.query_alignment_end)
    assert view.end == expected_end
    assert view.sequence == read.query_sequence[:expected_end]
    assert (read.seq, read.cigarstring, read.qual) == before


@pytest.mark.parametrize("clip", [False, True])
@pytest.mark.parametrize("compact", [False, True])
@pytest.mark.parametrize("gap", ["D", "N"])
def test_collector_slices_all_arrays_once_and_keeps_cigar_allele(clip, compact, gap):
    sequence = "T" * 18 + "ACG" + "TT" + "CG" + "GTACG" + "A" * 18
    read = make_pysam_read(sequence, "18S3M2I2M3%s5M18S" % gap)
    read.query_qualities = None
    collector = ReadCollector(use_soft_clipped_bases=clip, trim_poly_a=True)
    locus, = collector.collect_locus_reads_with_optional_compaction(
        MockAlignmentFile(["1"], [read]), "1", 3, 3, compact=compact)
    allele = AlleleRead.from_locus_read(locus)
    assert allele.allele == "TT"
    assert allele.prefix == "ACG" and allele.suffix == "CGGTACG"
    assert len(locus.sequence) == len(locus.quality_scores) == 12
    assert all(q is None for q in locus.quality_scores)
    view = allele.source_read_views[0][1]
    assert view.original_sequence == sequence and view.original_qualities is None
    assert view.original_interval(locus.read_base0_start_inclusive,
                                  locus.read_base0_end_exclusive) == (21, 23)
    assert allele.with_compatible_transcript_ids({"t1"}).source_read_views == allele.source_read_views


def test_terminal_insertion_does_not_become_a_tail_or_empty_allele():
    read = make_pysam_read(BODY + "A" * 26, "12M12I14S")
    collector = ReadCollector(use_soft_clipped_bases=True, trim_poly_a=True)
    locus = collector.locus_read_from_pysam_aligned_segment(read, 12, 12)
    assert AlleleRead.from_locus_read(locus).allele == "A" * 12
    assert locus.sequence == BODY + "A" * 12


def test_soft_clip_exclusion_and_annotation_use_one_retained_coordinate_map():
    read = make_pysam_read("CGCGG" + BODY + "A" * 20, "5S12M20S")
    collector = ReadCollector(infer_read_ends=True)
    locus = collector.locus_read_from_pysam_aligned_segment(read, 2, 3)
    view = locus.source_read_views[0][1]
    assert view.sequence == locus.sequence == BODY
    assert view.original_interval(locus.read_base0_start_inclusive,
                                  locus.read_base0_end_exclusive) == (7, 8)


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("hard_clip", [False, True])
def test_trimmed_supplementary_paths_keep_original_query_intervals(reverse, hard_clip):
    from .test_chimeric_phasing import record, link, results, assert_placement_phasing

    first = record()
    cigar = ("30M30H" if hard_clip else "30M30S") if reverse else ("30H30M" if hard_clip else "30S30M")
    second = record(200, cigar, 2048 | (16 if reverse else 0))
    original = link(first, second)
    before = [r.to_string() for r in original]
    plain = results(original, use_soft_clipped_bases=True)
    trimmed = results(original, use_soft_clipped_bases=True, trim_poly_a=True)
    for a, b in zip(plain, trimmed):
        assert a.alt_reads[0].source_alignment_paths == b.alt_reads[0].source_alignment_paths
        assert a.alt_reads[0].source_alignments == b.alt_reads[0].source_alignments
    assert_placement_phasing(trimmed, 1, {"split"})
    assert [r.to_string() for r in original] == before


def test_adapter_cannot_override_aligned_bases_or_insertions():
    sequence = BODY + MOTIF
    for cigar in (f"{len(sequence)}M", f"12M{len(MOTIF)}I", f"18M{len(MOTIF) - 6}S"):
        read = make_pysam_read(sequence, cigar)
        assert read_sequence_view_from_alignment(read, PROFILE, trim_adapters=True).sequence == sequence


def test_hard_clip_reverse_mapping_and_failed_upstream_estimate():
    read = make_pysam_read("T" * 20 + BODY, "7H20S12M9H")
    read.flag = 16
    read.set_tag("pt", -1)
    view = read_sequence_view_from_alignment(read, trim_poly_a=True)
    assert view.sequence == BODY
    assert view.hard_clips == (7, 9)
    assert view.sequenced_interval(0, 12) == (9, 21)
    assert view.upstream_poly_a_length == -1  # Not a zero-length biological tail.
    with pytest.raises(ValueError):
        view.original_interval(0, 13)


def test_mate_merge_retains_separate_original_views():
    first = make_pysam_read("T" * 18 + "ACCGTG", "18S6M", name="pair")
    second = make_pysam_read("CGTGCC" + "A" * 18, "6M18S", name="pair", reference_start=2)
    first.flag, second.flag = 65, 129
    collector = ReadCollector(use_soft_clipped_bases=True, trim_poly_a=True)
    merged, = collector.get_locus_reads(MockAlignmentFile(["1"], [first, second]), "1", 3, 4)
    assert merged.sequence == "ACCGTGCC"
    assert len(merged.source_read_views) == 2
    assert {v.original_sequence for _, v in merged.source_read_views} == {first.seq, second.seq}
    assert len(AlleleRead.from_locus_read(merged).source_read_views) == 2


def test_profiles_are_read_group_and_mate_specific():
    profile = ReadEndProfile("library", adapters=(Adapter("r2", MOTIF, mate=2),))
    collector = ReadCollector(read_end_profile={"known": profile}, trim_adapters=True)
    read = make_pysam_read(BODY + MOTIF, f"12M{len(MOTIF)}S")
    read.set_tag("RG", "known")
    read.flag = 65
    assert collector.read_sequence_view(read).sequence == read.seq
    read.flag = 129
    assert collector.read_sequence_view(read).sequence == BODY
    read.set_tag("RG", "unknown")
    assert collector.read_sequence_view(read).sequence == read.seq


def test_collection_infers_ends_only_after_locus_and_read_filters(monkeypatch):
    retained = [make_pysam_read(BODY, "12M", name=name, mapq=30, reference_start=100)
                for name in ("reference", "alternate")]
    retained[1].query_sequence = "GC" + "A" + BODY[3:]
    retained[1].query_qualities = [30] * len(BODY)
    duplicate = make_pysam_read(BODY, "12M", name="duplicate", mapq=30, reference_start=100)
    duplicate.flag = 1024
    rejected = [duplicate,
                make_pysam_read(BODY, "12M", name="distant", mapq=30, reference_start=300),
                make_pysam_read(BODY, "12M", name="low-mapq", mapq=0, reference_start=100),
                make_pysam_read(BODY, "2M5N10M", name="spliced-out", mapq=30, reference_start=100),
                make_pysam_read(BODY, "2M10S", name="no-allele", mapq=30, reference_start=100)]
    collector = ReadCollector(infer_read_ends=True, use_duplicate_reads=False, min_mapping_quality=20)
    seen, fetches = [], []
    original = collector.read_sequence_view

    def annotate(read):
        seen.append(read.query_name)
        return original(read)

    def fetch(*args):
        fetches.append(args)
        return retained + rejected

    bam = MockAlignmentFile(["1"], [])
    monkeypatch.setattr(bam, "fetch", fetch)
    monkeypatch.setattr(collector, "read_sequence_view", annotate)
    loci = collector.get_locus_reads(bam, "1", 102, 103)
    assert fetches == [("1", 101, 104)]
    assert seen == ["reference", "alternate"]
    assert {locus.name for locus in loci} == set(seen)


def test_cli_profile_roundtrip_and_shared_defaults(tmp_path):
    path = tmp_path / "profile.json"
    path.write_text(json.dumps({"read_groups": {"known": asdict(PROFILE)}}))
    args = make_isovar_arg_parser().parse_args([
        "--bam", "unused", "--read-end-profile", str(path), "--trim-adapters", "--trim-poly-a"])
    collector = read_collector_from_args(args)
    assert collector.read_end_profile == {"known": PROFILE}
    assert collector.trim_adapters and collector.trim_poly_a
    assert read_end_profiles_from_json(path) == collector.read_end_profile
    assert ReadEndProfile.from_dict(asdict(PROFILE)).fingerprint == PROFILE.fingerprint
    assert ReadEndProfile(PROFILE.name, PROFILE.version).fingerprint != PROFILE.fingerprint
    assert infer_read_ends(BODY, profile=PROFILE).profile[-1] == PROFILE.fingerprint
    with pytest.raises(ValueError, match="explicit"):
        ReadCollector(trim_adapters=True)


@pytest.mark.parametrize("kwargs", [
    {"sequence": "ACGT", "min_overlap": 12}, {"sequence": "ACGTNNNACGTA"},
    {"sequence": MOTIF, "max_error_rate": float("nan")}, {"sequence": MOTIF, "end": "internal"},
])
def test_invalid_adapter_profile_fails_explicitly(kwargs):
    with pytest.raises(ValueError):
        Adapter("invalid", **kwargs)


def test_actual_sid_records_preserved_and_tails_do_not_trim_mapped_a_bases():
    with pysam.AlignmentFile(Path(__file__).parent / "data/read_ends/osteosarc.sam") as bam:
        reads = list(bam)
    profile = ReadEndProfile("Illumina-TruSeq-R2", adapters=(
        Adapter("TruSeq-R2", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", mate=2),))
    collector = ReadCollector(read_end_profile=profile, trim_adapters=True, trim_poly_a=True)
    for read in reads:
        original = read.to_string()
        view = collector.read_sequence_view(read)
        assert view.original_sequence == read.query_sequence
        assert view.start <= read.query_alignment_start and view.end >= read.query_alignment_end
        if read.query_name.startswith("A01231"):
            assert view.end == 114
            assert view.annotations[0].errors == 1
            assert "AGATCGGAAGAGC" not in view.sequence
        elif read.query_name.startswith("231d"):
            assert view.end == 1002  # 5 aligned A bases are NOT removed.
            assert view.annotations[-1].start == 997
        elif read.query_name.startswith("a81a"):
            assert view.end == 2743  # Unexplained CT is NOT silently trimmed.
        else:
            assert view.original_qualities is None
            assert view.sequence == read.query_sequence
        assert read.to_string() == original
