"""Exact collection comparisons and bounded, non-aliasing coordinate sharing."""

from functools import lru_cache
import gc
import json
from pathlib import Path
import random
from types import SimpleNamespace
import weakref

import pysam
import pytest

from isovar.allele_read_helpers import allele_reads_from_locus_reads
from isovar.read_collector import ReadCollector
import isovar.read_collector as collector_module
from tests.mock_objects import MockAlignmentFile, make_pysam_read


class UnsharedCollector(ReadCollector):
    """Frozen copy of the Isovar 1.8.2 collection loop, omitting only logging.

    This is the differential baseline, not a mirror of the current collector.
    Do not edit it to follow later ``ReadCollector.get_locus_reads`` changes:
    a mismatch against it means collection behavior changed relative to 1.8.2,
    and must be reviewed as such rather than hidden by updating this copy.
    """

    def get_locus_reads(self, alignment_file, chromosome, base0_start_inclusive,
                        base0_end_exclusive, trimmed_base1_start=None, trimmed_ref=None, trimmed_alt=None):
        reads = []
        before = max(0, base0_start_inclusive - 1)
        after = base0_end_exclusive + 1
        for segment in alignment_file.fetch(chromosome, before, after):
            if segment.get_overlap(before, after) == 0:
                continue
            read = self.locus_read_from_pysam_aligned_segment(
                segment, base0_start_inclusive=base0_start_inclusive,
                base0_end_exclusive=base0_end_exclusive, trimmed_base1_start=trimmed_base1_start,
                trimmed_ref=trimmed_ref, trimmed_alt=trimmed_alt)
            if read is not None:
                reads.append(read)
        return self._merge_overlapping_locus_reads(reads) if self.merge_overlapping_fragments else reads


def test_frozen_182_baseline_loop_is_unchanged():
    """Editing the baseline must be a deliberate, reviewed change to this pin."""
    import hashlib
    import inspect

    source = inspect.getsource(UnsharedCollector.get_locus_reads)
    assert hashlib.sha256(source.encode()).hexdigest() == "1265dd36850d3bc525328201947ed2013b7df6f83b3b0b7fe58aedaeef0626af"


def assert_same_reads(left, right):
    assert type(left) is type(right) is list
    assert left == right  # ValueObject compares every field and class, in order.
    assert allele_reads_from_locus_reads(left) == allele_reads_from_locus_reads(right)


CORPUS = Path(__file__).parent / "data/osteosarc/expansion/corpus"
CASES = json.loads((CORPUS / "manifest.json").read_text())["cases"]


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["case_id"])
@pytest.mark.parametrize("mode", ["defaults", "primary_only"])
@pytest.mark.parametrize("merge", [False, True])
def test_original_locus_and_allele_reads_match_182(case, mode, merge):
    from varcode import Variant
    from isovar.variant_helpers import base0_interval_for_variant_fields, trim_variant

    record = case["variant"]
    position, ref, alt = trim_variant(Variant(record["chrom"], record["pos"], record["ref"], record["alt"]))
    start, end = base0_interval_for_variant_fields(position, ref, alt)
    results = []
    for cls in (UnsharedCollector, ReadCollector):
        with pysam.AlignmentFile(CORPUS / case["bam" if mode == "defaults" else "primary_bam"]) as bam:
            chromosome = cls._infer_chromosome_name(record["chrom"], bam.references)
            results.append(cls(merge_overlapping_fragments=merge).get_locus_reads(
                bam, chromosome, start, end, position, ref, alt))
    assert_same_reads(*results)


@pytest.mark.parametrize("seed", range(12))
@pytest.mark.parametrize("merge,soft_clips,secondary,duplicates,min_mapq", [
    (False, False, True, False, 1), (True, False, True, False, 1),
    (False, True, False, True, 0), (True, True, False, True, 0),
    (False, True, True, True, 30), (True, True, True, True, 30),
])
def test_randomized_flags_qualities_mates_splices_and_indels(
        seed, merge, soft_clips, secondary, duplicates, min_mapq):
    rng = random.Random(seed)
    reads = []
    cigars = ["60M", "20M2I38M", "20M2D40M", "20M100N40M", "3S54M3S", "2H60M2H", "20=2X38="]
    for i in range(35):
        read = make_pysam_read(
            seq="".join(rng.choices("ACGTN", k=60)), cigar=rng.choice(cigars),
            name=f"fragment-{rng.randrange(12)}", mapq=rng.choice([0, 1, 20, 30, 60, 255]),
            baseq=rng.choice([0, 10, 30, 40]), reference_start=10000 + rng.randrange(10))
        read.flag = rng.choice([0, 1 | 64, 1 | 128, 256, 1024, 2048])
        if i % 13 == 0:
            read.query_qualities = None
        reads.append(read)
    options = dict(merge_overlapping_fragments=merge, use_soft_clipped_bases=soft_clips,
                   use_secondary_alignments=secondary, use_duplicate_reads=duplicates, min_mapping_quality=min_mapq)
    for start, end, position, ref, alt in [(10025, 10026, 10026, "A", "C"),
                                          (10025, 10025, 10025, "", "AC"),
                                          (10025, 10027, 10026, "AA", ""),
                                          (10025, 10027, 10026, "AA", "C")]:
        sam = MockAlignmentFile(["1"], reads)
        assert_same_reads(*(cls(**options).get_locus_reads(sam, "1", start, end, position, ref, alt)
                            for cls in (UnsharedCollector, ReadCollector)))


def repeated_reads(n=3, length=20):
    return MockAlignmentFile(["1"], [make_pysam_read(
        "A" * length, f"{length}M", name=f"read-{i}", reference_start=10000) for i in range(n)])


def test_only_coordinate_integers_are_shared_and_hook_objects_survive():
    observed = []

    class ObservedCollector(ReadCollector):
        def locus_read_from_pysam_aligned_segment(self, *args, **kwargs):
            read = super().locus_read_from_pysam_aligned_segment(*args, **kwargs)
            observed.append((read, read.reference_positions, read.quality_scores))
            return read

    reads = ObservedCollector().get_locus_reads(repeated_reads(), "1", 10005, 10006)
    for read, (original, positions, qualities) in zip(reads, observed):
        assert read is original
        assert read.reference_positions is positions
        assert read.quality_scores is qualities
        assert positions == list(range(10000, 10020))
        assert read.source_read_count == 1
    assert len({id(r) for r in reads}) == 3
    assert len({id(r.reference_positions) for r in reads}) == 3
    for column in zip(*(r.reference_positions for r in reads)):
        assert len({id(p) for p in column}) == 1
    reads[0].reference_positions[0] = 999
    reads[0].quality_scores[0] = 1
    assert reads[1].reference_positions[0] == 10000
    assert reads[1].quality_scores[0] == 30


@pytest.mark.parametrize("positions_type", [tuple, type("CustomList", (list,), {})])
def test_custom_position_sequences_from_hooks_are_not_modified(positions_type):
    originals = []

    class CustomCollector(ReadCollector):
        def locus_read_from_pysam_aligned_segment(self, *args, **kwargs):
            read = super().locus_read_from_pysam_aligned_segment(*args, **kwargs)
            read.reference_positions = positions_type(read.reference_positions)
            originals.append(read.reference_positions)
            return read

    reads = CustomCollector().get_locus_reads(repeated_reads(), "1", 10005, 10006)
    assert all(r.reference_positions is p and type(p) is positions_type for r, p in zip(reads, originals))
    assert reads[0].reference_positions[0] is not reads[1].reference_positions[0]


def test_custom_integer_objects_from_hooks_keep_their_type_and_identity():
    class Coordinate(int):
        pass

    originals = []

    class CustomCollector(ReadCollector):
        def locus_read_from_pysam_aligned_segment(self, *args, **kwargs):
            read = super().locus_read_from_pysam_aligned_segment(*args, **kwargs)
            read.reference_positions[0] = Coordinate(read.reference_positions[0])
            originals.append(read.reference_positions[0])
            return read

    reads = CustomCollector().get_locus_reads(repeated_reads(), "1", 10005, 10006)
    assert all(r.reference_positions[0] is p and type(p) is Coordinate for r, p in zip(reads, originals))


def test_coordinate_cache_is_bounded_and_does_not_survive_collection(monkeypatch):
    caches, sizes = [], []

    def observed_cache(**options):
        def decorate(function):
            wrapped = lru_cache(**options)(function)
            caches.append(weakref.ref(wrapped))
            return wrapped
        return decorate

    class ObservedCollector(ReadCollector):
        @classmethod
        def _merge_overlapping_locus_reads(cls, reads):
            sizes.append(caches[-1]().cache_info())
            return super()._merge_overlapping_locus_reads(reads)

    monkeypatch.setattr(collector_module, "lru_cache", observed_cache)
    collector = ObservedCollector(merge_overlapping_fragments=True)
    length = 65540
    sam = repeated_reads(n=2, length=length)
    result = collector.get_locus_reads(sam, "1", 10005, 10006)
    assert_same_reads(result, UnsharedCollector(merge_overlapping_fragments=True).get_locus_reads(
        sam, "1", 10005, 10006))
    assert sizes[0].currsize == sizes[0].maxsize == 65536
    assert len(caches) == 1 and caches[0]() is None
    assert collector.__dict__ == ReadCollector(merge_overlapping_fragments=True).__dict__


def test_coordinate_cache_is_released_after_a_conversion_error(monkeypatch):
    caches = []

    def observed_cache(**options):
        def decorate(function):
            wrapped = lru_cache(**options)(function)
            caches.append(weakref.ref(wrapped))
            return wrapped
        return decorate

    class FailingCollector(ReadCollector):
        def locus_read_from_pysam_aligned_segment(self, read, *args, **kwargs):
            if read.query_name == "read-1":
                assert caches[-1]().cache_info().currsize == 20
                raise ValueError("conversion failure")
            return super().locus_read_from_pysam_aligned_segment(read, *args, **kwargs)

    monkeypatch.setattr(collector_module, "lru_cache", observed_cache)
    with pytest.raises(ValueError, match="conversion failure"):
        FailingCollector().get_locus_reads(repeated_reads(), "1", 10005, 10006)
    gc.collect()
    assert len(caches) == 1 and caches[0]() is None


def test_public_evidence_path_still_calls_each_existing_subclass_hook():
    calls = []

    class HookedCollector(ReadCollector):
        def locus_reads_overlapping_variant(self, *args, **kwargs):
            calls.append("variant")
            return super().locus_reads_overlapping_variant(*args, **kwargs)

        def get_locus_reads(self, *args, **kwargs):
            calls.append("locus")
            return super().get_locus_reads(*args, **kwargs)

        def locus_read_from_pysam_aligned_segment(self, segment, *args, **kwargs):
            calls.append(segment.query_name)
            return super().locus_read_from_pysam_aligned_segment(segment, *args, **kwargs)

        @classmethod
        def _merge_overlapping_locus_reads(cls, reads):
            assert type(reads) is list
            calls.append("merge")
            return super()._merge_overlapping_locus_reads(reads)

    variant = SimpleNamespace(contig="1", start=10006, ref="A", alt="C", gene_names=[])
    collector = HookedCollector(merge_overlapping_fragments=True)
    result = collector.read_evidence_for_variant(variant, repeated_reads())
    assert calls == ["variant", "locus", "read-0", "read-1", "read-2", "merge"]
    assert len(result.ref_reads) == 3
    assert result.alt_reads == result.other_reads == []


def test_separate_calls_do_not_share_a_persistent_coordinate_cache():
    collector = ReadCollector()
    first = collector.get_locus_reads(repeated_reads(), "1", 10005, 10006)
    second = collector.get_locus_reads(repeated_reads(), "1", 10005, 10006)
    assert_same_reads(first, second)
    assert first[0].reference_positions[0] is not second[0].reference_positions[0]
