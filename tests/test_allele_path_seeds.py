"""Allele identity and original placement geometry must agree before assembly."""

from types import SimpleNamespace

import pysam
import pytest

from isovar.allele_path_seeds import allele_path_seeds
from isovar.read_collector import ReadCollector
from isovar.rna_reconstruction import _assemble
from isovar.variant_helpers import LiteralVariantLocus


HEADER = pysam.AlignmentHeader.from_dict({
    "SQ": [{"SN": "chr1", "LN": 1000}],
    "RG": [{"ID": "A", "SM": "sample"}, {"ID": "B", "SM": "sample"}],
})


def record(sequence="CCCTAACCC", cigar="9M", *, reverse=False, start=100, name="original", group="A", md=None):
    read = pysam.AlignedSegment(HEADER)
    read.query_name = name
    read.query_sequence = sequence
    read.query_qualities = [30] * len(sequence)
    read.reference_id = 0
    read.reference_start = start
    read.cigarstring = cigar
    read.flag = 16 if reverse else 0
    read.mapping_quality = 60
    read.set_tag("RG", group)
    if md is not None:
        read.set_tag("MD", md)
    return read


def collect(templates, locus=LiteralVariantLocus(104, "A", "T"), **options):
    return allele_path_seeds(templates, locus, "chr1", read_collector=ReadCollector(**options))


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("cigar", ["9M", "5H9M7H"])
def test_orientation_and_hard_clips_preserve_original_allele(reverse, cigar):
    read = record(reverse=reverse, cigar=cigar)
    original = read.to_string()
    result = collect([[read]])
    assert len(result.seeds) == 2 and result.fragments("ALT") == {("A", "original")}
    for seed in result.seeds:
        assert seed.state == "ALT" and seed.allele == "T"
        start, end = seed.query_interval
        strand = seed.observation.positions[start][2]
        assert seed.observation.sequence[start:end] == ("T" if strand == "+" else "A")
        assert seed.canonical_query_interval == seed.query_interval
    assert read.to_string() == original
    assert set(dict(result.original_records).values()) == {original}


@pytest.mark.parametrize("reverse", [False, True])
def test_repeat_normalization_keeps_actual_and_canonical_boundaries(reverse):
    read = record("CCCAAACCC", "6M1D3M", reverse=reverse, md="6^A3")
    result = collect([[read]], LiteralVariantLocus(104, "A", ""))
    assert len(result.seeds) == 2 and len(result.fragments("ALT")) == 1
    for seed in result.seeds:
        strand = next(position[2] for position in seed.observation.positions if position)
        assert seed.query_interval == ((6, 6) if strand == "+" else (3, 3))
        assert seed.canonical_query_interval == ((3, 3) if strand == "+" else (6, 6))
        assert seed.state == "ALT"


def test_alternative_placements_can_disagree_without_adding_molecules():
    alt = record()
    ref = record(start=99)
    ref.flag |= 256
    other_group = record(group="B")
    result = collect([[alt, ref, alt], [other_group]], use_secondary_alignments=True)
    assert result.fragments("ALT") == {("A", "original"), ("B", "original")}
    assert result.fragments("REF") == {("A", "original")}
    assert len(result.segments("ALT")) == 2
    assert not result.fragments("OTHER")
    assert len(result.original_records) == 3


def test_skipped_intron_is_not_an_observed_deletion():
    spliced = record("CCCTAACCC", "3M20N6M")
    retained = record("CCCTAACCC", "9M", name="retained")
    result = collect([[spliced], [retained]], LiteralVariantLocus(104, "A", ""))
    assert not result.fragments("ALT")
    assert result.fragments("OTHER") == {("A", "retained")}


def test_trimmed_terminal_deletion_does_not_borrow_missing_flanks():
    read = record("AACCC", "2S1D3M", md="0^A3")
    result = collect([[read]], LiteralVariantLocus(101, "A", ""))
    assert not result.seeds


@pytest.mark.parametrize("start,ref,alt,expected", [
    (104, "A", "T", (103, 104, "A", "T")),
    (104, "GA", "G", (104, 105, "A", "")),
    (104, "G", "GAA", (104, 104, "", "AA")),
])
def test_normalized_locus_conventions_match_collector(start, ref, alt, expected):
    locus = LiteralVariantLocus.from_variant(SimpleNamespace(start=start, ref=ref, alt=alt))
    assert (locus.start, locus.end, locus.ref, locus.alt) == expected


def test_explicit_allele_seeds_use_shared_extension_without_unrelated_paths():
    seed_read = record(sequence="ACGT", cigar="4M")
    background_read = record(sequence="GTAA", cigar="4M", start=102, name="background")
    target = collect([[seed_read]], LiteralVariantLocus(102, "A", "C"))
    background = collect([[background_read]], LiteralVariantLocus(104, "A", "T"))
    seeds = tuple(seed.observation for seed in target.seeds
                  if seed.observation.positions[0][2] == "+")
    pool = {seed.observation.key: seed.observation for seed in background.seeds
            if seed.observation.positions[0][2] == "+"}
    paths, _ = _assemble(pool, 10, 2, seed_observations=seeds)
    assert len(paths) == 1 and paths[0][0] == "ACGTAA"
    assert _assemble(pool, 10, 2, seed_observations=())[0] == []


@pytest.mark.parametrize("start,ref,alt", [(True, "A", "T"), (0, "A", "T"), (-1, "", "T"),
                                         (1, "A", "<DEL>"), (1, None, "T")])
def test_invalid_literal_locus_is_explicit(start, ref, alt):
    with pytest.raises(ValueError):
        LiteralVariantLocus(start, ref, alt)


def test_allele_count_survives_unusable_ambiguous_flanking_sequence():
    read = record(sequence='NCCTAACCC')
    result = collect([[read]])
    assert result.fragments('ALT') == {('A', 'original')}
    assert not result.seed_segments('ALT')
    assert not result.seeds
    assert set(dict(result.original_records).values()) == {read.to_string()}
    call, = result.calls
    assert call.state == 'ALT' and call.allele == 'T'
    assert dict(result.notes)['ambiguous_bases'] >= 1
