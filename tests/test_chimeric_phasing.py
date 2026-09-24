"""#286: supplementary pieces, paired mates, and alternatives are distinct."""

from pathlib import Path

import pysam
import pytest
from varcode import Variant

from isovar.chimeric_alignment import (
    _alignment, _matches_reported, source_alignment_paths_from_pysam,
)
from isovar.read_collector import ReadCollector, _CompactLocusRead
from isovar.read_identity import (
    alignment_constraints, compatible_alignments, count_reads, fragment_ids,
    source_alignments_from_pysam,
)
from isovar.allele_read import AlleleRead
from isovar.phasing import annotate_phased_variants
from tests.mock_objects import MockAlignmentFile
from tests.test_phasing import DummyProteinSequence, assert_placement_phasing
from isovar.isovar_result import IsovarResult


HEADER = pysam.AlignmentHeader.from_references(["1", "2"], [10000, 10000])


def record(start=100, cigar="30M30S", flag=0, name="split", group="a", contig="1"):
    read = pysam.AlignedSegment(HEADER)
    read.query_name = name
    read.reference_id = HEADER.get_tid(contig)
    read.reference_start = start
    read.cigarstring = cigar
    read.flag = flag
    read.query_sequence = "A" * read.infer_query_length()
    read.query_qualities = [30] * len(read.query_sequence)
    read.mapping_quality = 60
    read.set_tag("RG", group)
    return read


def sa_entry(read, cigar=None):
    return "%s,%d,%s,%s,60,0;" % (
        read.reference_name, read.reference_start + 1, "-" if read.is_reverse else "+",
        cigar or read.cigarstring.replace("H", "S"))


def link(*reads):
    for read in reads:
        read.set_tag("SA", "".join(sa_entry(other) for other in reads if other is not read))
    return reads


def results(records, loci=(("1", 115), ("1", 215)), **options):
    collector = ReadCollector(**options)
    values = []
    for contig, position in loci:
        variant = Variant(contig, position, "T", "A", genome="GRCh38")
        # MockAlignmentFile does not filter fetches by contig.
        bam = MockAlignmentFile(HEADER.references, [r for r in records if r.reference_name == contig])
        evidence = collector.read_evidence_for_variant(variant, bam)
        protein = DummyProteinSequence({r.name for r in evidence.alt_reads}, "K", (), (), ())
        protein.supporting_reads = set(evidence.alt_reads)
        values.append(IsovarResult(variant, evidence, None, sorted_protein_sequences=[protein]))
    return values


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("hard_clip", [False, True])
@pytest.mark.parametrize("paired", [False, True])
def test_reciprocal_split_path_phases_once_not_as_mates(reverse, hard_clip, paired):
    first = record(flag=65 if paired else 0)
    cigar = ("30M30H" if hard_clip else "30M30S") if reverse else ("30H30M" if hard_clip else "30S30M")
    second = record(200, cigar, 2048 | (16 if reverse else 0) | (65 if paired else 0))
    inputs = results(link(first, second))
    assert all(r.num_alt_reads == r.num_alt_fragments == 1 for r in inputs)
    reads = inputs[0].alt_reads + inputs[1].alt_reads
    assert count_reads(reads) == len(fragment_ids(reads)) == 1
    assert not compatible_alignments(alignment_constraints(reads[:1]), reads[1:])
    for ordered in (inputs, inputs[::-1]):
        assert_placement_phasing(ordered, 1, {"split"})
        assert_placement_phasing(ordered, 2, set())


def test_reciprocal_cross_contig_path_and_two_fragment_threshold():
    reads = []
    for name in ("split1", "split2"):
        reads.extend(link(record(name=name), record(200, "30H30M", 2048, name, contig="2")))
    inputs = results(reads, (("1", 115), ("2", 215)))
    assert_placement_phasing(inputs, 2, {"split1", "split2"})
    assert_placement_phasing(inputs, 3, set())


@pytest.mark.parametrize("secondary", [False, True])
def test_two_supplementary_pieces_can_share_an_unobserved_representative(secondary):
    flag = 256 if secondary else 0
    first, middle, last = link(record(cigar="20M40S", flag=flag),
                               record(200, "20H20M20H", 2048 | flag),
                               record(300, "40H20M", 2048 | flag))
    # The representative need not overlap either requested variant.
    assert_placement_phasing(results([middle, last], (("1", 210), ("1", 310))), 1, {"split"})


@pytest.mark.parametrize("flags", [(0, 0), (0, 256), (0, 256 | 2048), (256, 2048)])
def test_sa_does_not_turn_alternatives_into_supplementary_pieces(flags):
    inputs = results(link(record(flag=flags[0]), record(200, "30S30M", flags[1])))
    assert_placement_phasing(inputs, 1, set())


def test_valid_secondary_chimeric_hypothesis_remains_available():
    inputs = results(link(record(flag=256), record(200, "30H30M", 256 | 2048)))
    assert_placement_phasing(inputs, 1, {"split"})


@pytest.mark.parametrize("change", ["missing", "one_sided", "wrong_start", "wrong_strand", "other_path"])
def test_missing_or_conflicting_path_links_do_not_phase(change):
    first, second = link(record(), record(200, "30S30M", 2048))
    if change in ("missing", "one_sided"):
        second.set_tag("SA", None)
        if change == "missing":
            first.set_tag("SA", None)
    elif change == "wrong_start":
        first.set_tag("SA", "1,202,+,30S30M,60,0;")
    elif change == "wrong_strand":
        first.set_tag("SA", "1,201,-,30M30S,60,0;")
    else:
        second.set_tag("SA", "1,301,+,30M30S,60,0;")
    assert_placement_phasing(results([first, second]), 1, set())


@pytest.mark.parametrize("tag", ["", "bad", "1,201,+,30S30M,60;", "1,0,+,30S30M,60,0;",
                                "1,201,?,30S30M,60,0;", "1,201,+,30S30M,256,0;",
                                "1,201,+,30S30M,60,-1;", "missing,201,+,30S30M,60,0;",
                                "1,201,+,30S31M,60,0;", "1,201,+,30M30S,60,0;",
                                "1,201,+,0M30S30M,60,0;", "1,201,+,10M10S40M,60,0;"])
def test_malformed_sa_is_uncertain_without_discarding_allele_evidence(tag):
    first, second = link(record(), record(200, "30S30M", 2048))
    first.set_tag("SA", tag)
    inputs = results([first, second])
    assert all(r.num_alt_fragments == 1 for r in inputs)
    assert_placement_phasing(inputs, 1, set())


def test_sa_alone_never_creates_an_observed_partner():
    first, _ = link(record(), record(200, "30H30M", 2048))
    inputs = results([first])
    assert inputs[0].num_alt_fragments == 1 and inputs[1].num_alt_fragments == 0
    assert_placement_phasing(inputs, 1, set())


def test_path_does_not_cross_read_groups():
    first, second = link(record(), record(200, "30H30M", 2048, group="b"))
    assert_placement_phasing(results([first, second]), 1, set())


def test_mates_do_not_need_or_share_a_supplementary_path():
    first, second = record(flag=65), record(200, "30H30M", 129)
    inputs = results([first, second])
    assert_placement_phasing(inputs, 1, {"split"})
    reads = inputs[0].alt_reads + inputs[1].alt_reads
    assert count_reads(reads) == 2 and len(fragment_ids(reads)) == 1


@pytest.mark.parametrize("soft_clips", [False, True])
def test_overlapping_pieces_phase_only_outside_ambiguous_query_overlap(soft_clips):
    first, second = link(record(cigar="35M25S"), record(200, "30H30M", 2048))
    assert_placement_phasing(results([first, second], use_soft_clipped_bases=soft_clips), 1, {"split"})
    assert_placement_phasing(results([first, second], (("1", 133), ("1", 215)),
                                    use_soft_clipped_bases=soft_clips), 1, set())


def test_minimap2_span_summary_does_not_replace_actual_indel_cigar():
    first, second = link(record(cigar="30M31S"), record(200, "30H10M1I20M", 2048))
    first.set_tag("SA", "1,201,+,30S30M1I,60,1;")
    inputs = results([first, second])
    assert_placement_phasing(inputs, 1, {"split"})
    assert inputs[1].alt_reads[0].source_alignments[0][1][2] == "30H10M1I20M"


def test_approximate_sa_cannot_choose_between_observed_alternative_cigars():
    first, second = link(record(cigar="30M31S"), record(200, "30H10M1I20M", 2048))
    first.set_tag("SA", "1,201,+,30S30M1I,60,1;")
    alternative = record(200, "30H20M1I10M", 256)
    assert_placement_phasing(results([first, second, alternative]), 1, set())
    # An exact SA link can identify the intended piece among alternatives.
    first.set_tag("SA", sa_entry(second))
    assert_placement_phasing(results([first, second, alternative]), 1, {"split"})


def test_summary_equal_to_one_actual_cigar_still_cannot_choose_an_alternative():
    first, second = link(record(cigar="30M31S"), record(200, "30H30M1I", 2048))
    alternative = record(200, "30H10M1I20M", 256)
    assert_placement_phasing(results([first, second, alternative]), 1, set())


def test_paths_survive_compaction_classification_and_primary_mate_merging():
    first, supplementary = link(record(cigar="40M20S", flag=65), record(200, "40H20M", 2048 | 65))
    mate = record(start=110, cigar="30M", flag=129)
    collector = ReadCollector()
    locus = collector.locus_read_from_pysam_aligned_segment(first, 114, 115)
    compact = _CompactLocusRead.from_locus_read(locus)
    assert compact.source_alignment_paths == locus.source_alignment_paths
    allele = AlleleRead.from_locus_read(compact)
    assert allele.with_compatible_transcript_ids({"tx"}).source_alignment_paths == allele.source_alignment_paths
    inputs = results([first, supplementary, mate], (("1", 115), ("1", 210)))
    assert len(inputs[0].alt_reads) == 1 and inputs[0].num_alt_reads == 2
    assert inputs[0].alt_reads[0].source_alignment_paths
    assert_placement_phasing(inputs, 1, {"split"})


def test_protein_phasing_uses_its_own_actual_path_observations():
    first, second = link(record(), record(200, "30H30M", 2048))
    alternative = record(205, "30H30M", 256)
    inputs = results([first, second, alternative])
    inputs[1].top_protein_sequence.supporting_reads = {
        r for r in inputs[1].alt_reads if r.source_alignments[0][1][1] == 205}
    for result in annotate_phased_variants(inputs, 1):
        assert result.phased_variants_in_supporting_reads
        assert not result.phased_variants_in_protein_sequence


def test_cigar_normalization_preserves_interior_mapping_and_original_orientation():
    forward = _alignment("1", 100, False, "5H5S10=1X9M30S")
    assert forward.cigar == ((10, "S"), (20, "M"), (30, "S"))
    assert forward.query_start == 10 and forward.query_end == 30 and forward.query_length == 60
    reverse = _alignment("1", 100, True, "5H5S20M30H")
    assert reverse.query_start == 30 and reverse.query_end == 50
    assert _matches_reported(forward, _alignment("1", 100, False, "10S20M30S"))


def test_source_metadata_keeps_segment_flags_separate_from_alignment_flags():
    read, _ = link(record(flag=65 | 256), record(200, "30H30M", 65 | 256 | 2048))
    source = source_alignments_from_pysam(read, read.query_name)
    ((segment, path),) = source_alignment_paths_from_pysam(read, source)
    assert segment == ("a", "split", 64)
    assert path.secondary and not path.supplementary
    assert source == ((segment, (0, 100, "30M30S", False)),)


@pytest.mark.parametrize("kind", ["ont", "short"])
def test_original_osteosarc_rna_records_phase_only_with_observed_reciprocal_links(kind):
    with pysam.AlignmentFile(Path(__file__).parent / "data/chimeric" / f"osteosarc-{kind}.sam") as sam:
        records = [r for r in sam if r.has_tag("SA")]
        references = sam.references
    assert len(records) == 2
    # Probe one actual mapped base on each piece. These deliberately artificial
    # alleles test alignment co-occurrence, not a claim of patient mutations.
    def collect():
        values = []
        for record in records:
            pairs = record.get_aligned_pairs(matches_only=True)
            query, position = pairs[len(pairs) // 2]
            alt = record.query_sequence[query]
            ref = "A" if alt != "A" else "C"
            variant = Variant(record.reference_name.removeprefix("chr"), position + 1, ref, alt, genome="GRCh38")
            bam = MockAlignmentFile(references, [r for r in records if r.reference_id == record.reference_id])
            evidence = ReadCollector().read_evidence_for_variant(variant, bam)
            assert len(evidence.alt_reads) == 1
            protein = DummyProteinSequence({record.query_name}, "K", (), (), ())
            protein.supporting_reads = set(evidence.alt_reads)
            values.append(IsovarResult(variant, evidence, None, sorted_protein_sequences=[protein]))
        return values
    inputs = collect()
    assert_placement_phasing(inputs, 1, {records[0].query_name})
    assert_placement_phasing(inputs, 2, set())
    assert count_reads(inputs[0].alt_reads + inputs[1].alt_reads) == 1
    records[0].set_tag("SA", None)
    assert_placement_phasing(collect(), 1, set())


def test_sa_entries_parse_every_field_and_reject_malformed_tags():
    from isovar.chimeric_alignment import sa_entries
    assert sa_entries("chr2,101,-,20S30M,60,1;chr3,5,+,30M20S,0,0;") == [
        ("chr2", 100, "-", "20S30M", 60, 1), ("chr3", 4, "+", "30M20S", 0, 0)]
    for tag in ("chr2,101", "chr2,x,-,30M,60,0", 7):
        with pytest.raises(ValueError):
            sa_entries(tag)
