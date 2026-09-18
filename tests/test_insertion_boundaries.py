"""One-sided reads are retained at a locus, but do not prove reference alleles."""

import pysam
import pytest
from varcode import Variant

from isovar import ReadCollector
from isovar.allele_read import AlleleRead
from .genomes_for_testing import grch38
from .mock_objects import MockAlignmentFile, make_pysam_read


@pytest.mark.parametrize("left", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("clipped", [False, True])
@pytest.mark.parametrize("use_clips", [False, True])
@pytest.mark.parametrize("compact", [False, True])
def test_one_sided_insertion_locus_is_not_reference(tmp_path, left, reverse, clipped, use_clips, compact):
    insert = "GT" * 64
    clip = insert if clipped else ""
    sequence = "ACCA" + clip if left else clip + "CGGT"
    cigar = ("4M" + ("128S" if clipped else "") if left else
             ("128S" if clipped else "") + "4M")
    read = make_pysam_read(sequence, cigar, reference_start=100 if left else 104, mapq=60)
    read.flag = 16 if reverse else 0
    path = tmp_path / "original.bam"
    with pysam.AlignmentFile(path, "wb", header={"SQ": [{"SN": "1", "LN": 1000}]}) as bam:
        bam.write(read)
    pysam.index(str(path))
    collector = ReadCollector(use_soft_clipped_bases=use_clips)
    with pysam.AlignmentFile(path) as bam:
        locus, = collector.collect_locus_reads_with_optional_compaction(bam, "1", 104, 104, compact=compact)
        assert locus.sequence == sequence if use_clips else locus.sequence == ("ACCA" if left else "CGGT")
        assert AlleleRead.from_locus_read(locus) is None
        evidence = collector.read_evidence_for_variant(Variant("1", 104, "A", "A"+insert, ensembl=grch38), bam)
        assert not evidence.ref_reads and not evidence.alt_reads and not evidence.other_reads
    assert read.query_sequence == sequence  # No realignment/rewrite of input.


@pytest.mark.parametrize("left", [False, True])
@pytest.mark.parametrize("use_clips", [False, True])
@pytest.mark.parametrize("compact", [False, True])
def test_terminal_cigar_insertion_remains_an_observed_allele(left, use_clips, compact):
    insert = "GT" * 64
    seq = "TT" + ("ACCA"+insert if left else insert+"CGGT") + "AA"
    cigar = "2S" + ("4M128I" if left else "128I4M") + "2S"
    read = make_pysam_read(seq, cigar, reference_start=100 if left else 104, mapq=60)
    collector = ReadCollector(use_soft_clipped_bases=use_clips)
    locus, = collector.collect_locus_reads_with_optional_compaction(
        MockAlignmentFile(["1"], [read]), "1", 104, 104, compact=compact)
    assert AlleleRead.from_locus_read(locus).allele == insert


@pytest.mark.parametrize("cigar,seq,start,expected", [
    ("4M2D4M", "ACCAGTGT", 100, None),
    ("4M2D4M", "GTGTACCA", 98, None),
    ("4M2I2I2D", "ACCAGTGT", 100, "GTGT"),
    ("2D2I2I4M", "GTGTACCA", 102, "GTGT"),
    ("4M2D2I2I", "ACCAGTGT", 100, None),
    ("2I2I2D4M", "GTGTACCA", 102, None),
])
@pytest.mark.parametrize("use_clips", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
def test_deletion_at_an_anchor_never_turns_aligned_bases_into_an_insertion(cigar, seq, start, expected, use_clips, reverse):
    read = make_pysam_read(seq, cigar, reference_start=start, mapq=60)
    read.flag = 16 if reverse else 0
    locus = ReadCollector(use_soft_clipped_bases=use_clips).locus_read_from_pysam_aligned_segment(read, 104, 104)
    assert locus is not None  # Keep unassigned locus evidence, not a false allele.
    allele = AlleleRead.from_locus_read(locus)
    assert (allele.allele if allele else None) == expected


@pytest.mark.parametrize("merge", [False, True])
@pytest.mark.parametrize("alternative", [False, True])
@pytest.mark.parametrize("use_clips", [False, True])
def test_unassigned_end_does_not_cancel_a_spanning_mate_or_placement(merge, alternative, use_clips):
    first = make_pysam_read("ACCA", "4M", reference_start=100, mapq=60)
    second = make_pysam_read("ACCATGCA", "8M", reference_start=100, mapq=60)
    first.flag, second.flag = (0, 256) if alternative else (65, 129)
    collector = ReadCollector(use_soft_clipped_bases=use_clips, merge_overlapping_fragments=merge)
    variant = Variant("1", 104, "A", "AGT", ensembl=grch38)
    evidence = collector.read_evidence_for_variant(variant, MockAlignmentFile(["1"], [first, second]))
    assert len(evidence.ref_reads) == 1
    assert not evidence.alt_reads and not evidence.other_reads
    assert evidence.ref_reads[0].source_read_count == (2 if merge and not alternative else 1)


def test_original_sid_ktn1_endpoints_and_false_alternative_conflict():
    import gzip
    import hashlib
    import json
    from pathlib import Path
    from isovar.read_identity import fragment_ids, source_read_ids

    corpus = Path(__file__).parent / "data/osteosarc/figure_comparisons/corpus"
    manifest = json.loads((corpus / "insertion-boundaries-manifest.json").read_text())
    raw = (corpus / manifest["file"]).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == manifest["sha256"]
    data = json.loads(gzip.decompress(raw))
    variant = Variant(ensembl=grch38, **data["variant"])
    for product in data["products"]:
        header = pysam.AlignmentHeader.from_dict(product["header"])
        reads = [pysam.AlignedSegment.fromstring(s, header) for s in product["original_sam"]]
        evidence = ReadCollector().read_evidence_for_variant(variant, MockAlignmentFile(header.references, reads))
        actual = {key: sorted({tuple(s) for r in getattr(evidence, key+"_reads") for s in source_read_ids(r)})
                  for key in ("ref", "alt", "other")}
        assert actual == {k: [tuple(s) for s in v] for k, v in product["expected_source_ids"].items()}
        assert actual != {k: [tuple(s) for s in v] for k, v in product["before_source_ids"].items()}
        if product["source"] == "T2-short":
            assert len(fragment_ids(evidence.alt_reads)) == 1  # Two mates, one template ID.
            assert not evidence.other_reads
