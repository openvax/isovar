"""#264: placements of one segment are not extra reads or mates."""

import pytest
from varcode import Variant

from isovar.isovar_result import IsovarResult
from isovar.read_collector import ReadCollector
from isovar.read_evidence import ReadEvidence
from isovar.variant_sequence_creator import VariantSequenceCreator
from tests.mock_objects import MockAlignmentFile, make_pysam_read


def record(name="template", start=100, flag=0, group=None):
    read = make_pysam_read("A" * 60, "60M", name=name, mapq=3, reference_start=start)
    read.flag = flag
    if group is not None:
        read.set_tag("RG", group)
    return read


def evidence(records, **options):
    variant = Variant("1", 135, "T", "A", genome="GRCh38")
    reads = ReadCollector(**options).allele_reads_overlapping_variant(
        variant, MockAlignmentFile(["1"], records))
    grouped = ReadEvidence.from_variant_and_allele_reads(variant, reads)
    return variant, grouped, IsovarResult(variant, grouped, None)


@pytest.mark.parametrize("alternative_flag", [0, 256, 2048, 256 | 2048])
@pytest.mark.parametrize("merge", [False, True])
@pytest.mark.parametrize("assembly", [False, True])
def test_alternative_placements_never_extend_or_inflate_support(alternative_flag, merge, assembly):
    records = [record(name, start, flag)
               for name in ["template1", "template2"]
               for start, flag in [(100, 0), (118, alternative_flag)]]
    variant, grouped, result = evidence(records, merge_overlapping_fragments=merge)
    assert result.num_alt_reads == result.num_total_reads == 2
    assert result.num_alt_fragments == 2
    from isovar.filtering import evaluate_threshold_filters
    assert not evaluate_threshold_filters(result, {"min_num_alt_reads": 3})["min_num_alt_reads"]
    sequences = VariantSequenceCreator(
        preferred_sequence_length=180, variant_sequence_assembly=assembly,
    ).reads_to_variant_sequences(variant, grouped.alt_reads)
    assert sequences
    assert max(map(len, sequences)) == 60
    assert all(max(s.coverage()) == 2 for s in sequences)


def test_complementary_primary_mates_still_merge():
    _, grouped, result = evidence([record(flag=65), record(start=118, flag=129)])
    assert len(grouped.alt_reads) == 1
    assert len(grouped.alt_reads[0]) == 78
    assert result.num_alt_reads == 2
    assert result.num_alt_fragments == 1


@pytest.mark.parametrize("second_flag", [65, 65 | 256, 65 | 2048, 129 | 256, 129 | 2048, 193])
def test_only_complementary_primary_mates_collapse(second_flag):
    _, grouped, _ = evidence([record(flag=65), record(start=118, flag=second_flag)])
    assert len(grouped.alt_reads) == 2
    assert {len(r) for r in grouped.alt_reads} == {60}


def test_read_groups_separate_fragments_without_rewriting_names():
    _, grouped, result = evidence([record(flag=65, group="a"), record(flag=129, group="b")])
    assert len(grouped.alt_reads) == 2
    assert result.num_alt_reads == result.num_alt_fragments == 2
    assert result.alt_read_names == {"template"}


def test_secondary_only_evidence_is_retained():
    _, grouped, result = evidence([record(flag=256)])
    assert len(grouped.alt_reads) == result.num_alt_reads == 1
    assert evidence([record(flag=256)], use_secondary_alignments=False)[2].num_alt_reads == 0


def test_conflicting_alleles_on_one_segment_are_not_definitive_alt_support():
    first, second = record(), record(flag=256)
    second.query_sequence = "A" * 34 + "T" + "A" * 25
    second.query_qualities = first.query_qualities
    _, grouped, result = evidence([first, second])
    assert not grouped.alt_reads and not grouped.ref_reads
    assert len(grouped.other_reads) == 2
    assert result.num_other_reads == result.num_total_reads == 1


@pytest.mark.parametrize("assembly", [False, True])
def test_duplicate_alignment_records_count_once_everywhere(assembly):
    from isovar.dataframe_helpers import allele_counts_dataframe
    from isovar.read_identity import count_reads

    records = [record(flag=flag) for flag in [0, 256, 2048]]
    variant, grouped, result = evidence(records)
    assert result.num_alt_reads == result.num_total_reads == 1
    frame = allele_counts_dataframe([(variant, grouped)])
    assert frame.iloc[0].num_alt_reads == frame.iloc[0].num_alt_fragments == 1
    sequences = VariantSequenceCreator(
        min_variant_sequence_coverage=1, preferred_sequence_length=180,
        variant_sequence_assembly=assembly,
    ).reads_to_variant_sequences(variant, grouped.alt_reads)
    assert sequences and all(count_reads(s.reads) == 1 for s in sequences)
    assert all(max(s.coverage()) == 1 for s in sequences)
    assert not VariantSequenceCreator(min_variant_sequence_coverage=2).reads_to_variant_sequences(
        variant, grouped.alt_reads)


def test_identical_cdna_on_competing_alignment_paths_stays_separate():
    reads = [record(), record(flag=256)]
    reads[1].cigarstring = "40M1D20M"
    variant, grouped, _ = evidence(reads)
    for ordered in [grouped.alt_reads, grouped.alt_reads[::-1]]:
        sequences = VariantSequenceCreator(
            preferred_sequence_length=180, min_variant_sequence_coverage=1,
        ).reads_to_variant_sequences(variant, ordered)
        assert len(sequences) == 2
        assert len({s.sequence for s in sequences}) == 1
        assert all(len(s._source_alignments) == 1 and max(s.coverage()) == 1 for s in sequences)


def test_independent_secondary_evidence_can_extend_another_read():
    variant, grouped, _ = evidence([record(name="first"), record(name="second", start=118, flag=256)])
    sequences = VariantSequenceCreator(
        preferred_sequence_length=180, min_variant_sequence_coverage=1,
    ).reads_to_variant_sequences(variant, grouped.alt_reads)
    assert max(map(len, sequences)) == 78


def test_merged_pair_and_retained_single_view_do_not_double_count_protein_support():
    from isovar.protein_sequence import ProteinSequence
    from isovar.variant_sequence import VariantSequence
    from tests.test_protein_sequences import _translation_with_reads

    _, grouped, result = evidence([record(flag=65), record(start=118, flag=129), record(flag=65 | 256)])
    assert result.num_alt_reads == 2 and result.num_alt_fragments == 1
    merged = max(grouped.alt_reads, key=len)
    sequence = VariantSequence(merged.prefix, merged.allele, merged.suffix, grouped.alt_reads)
    assert set(sequence.coverage()) == {1}
    protein = ProteinSequence.from_translations([_translation_with_reads(
        amino_acids="MKH", reads=grouped.alt_reads, cdna_sequence="ATGAAACAC")])
    assert protein.num_supporting_reads == 2 and protein.num_supporting_fragments == 1


def test_ambiguous_collapsed_pair_does_not_count_its_mate_again_as_definitive():
    first, mate = record(flag=65), record(start=118, flag=129)
    competing, mate_copy = record(flag=65|256), record(start=118, flag=129|256)
    competing.query_sequence = "A"*34 + "T" + "A"*25
    competing.query_qualities = first.query_qualities
    _, grouped, result = evidence([first, mate, competing, mate_copy])
    assert not grouped.alt_reads and not grouped.ref_reads
    assert result.num_total_reads == result.num_other_reads == 2
