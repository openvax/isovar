"""Independent original-read replays of the HCC1395 fusion protein candidate."""

from pathlib import Path
from hashlib import sha256

import pysam
import pytest

from isovar.fusion import FusionBreakpoint
from isovar.read_collector import ReadCollector
from isovar.sv_rna import reconstruct_sv_rna
from reports.hcc1395_colo829.build_hcc1395_eif3k_fusion import load_references, report


DATA = Path(__file__).parent / "data/hcc1395_eif3k_cyp39a1"
PROTEIN = "FICHVVGITYQHIDRWLLAEMLGDLSASIPKNVFLALHEKLYIMLKGKMGTVNLH"


def test_original_illumina_reads_assemble_an_observed_partial_fusion_protein():
    value = report(DATA / "original_unmapped.bam", DATA / "models.json.gz",
                   DATA / "ont_candidate_sequences.fasta.gz")
    fusion, result = value["fusion"], value["result"]
    assert len(fusion.sequence) == 167
    assert fusion.junction_start == fusion.junction_end == 79
    assert result["status"] == "ambiguous"  # Noncoding donor models remain plausible.
    assert result["paths"][0]["junctions"][0]["direct_support"]["fragments"] == 149
    assert value["proteins"] == [PROTEIN]
    translation, = result["paths"][0]["translations"]
    assert not translation["complete_5prime"]
    assert translation["acceptor_frames"] == [dict(transcript_id="ENST00000275016", in_frame=True)]
    assert value["ont"] == dict(screened_original_reads=21, exact_junction_reads=6, flank_bases=20)


def test_coding_tiling_yields_a_complete_candidate_protein():
    value = report(DATA / "original_unmapped.bam", DATA / "models.json.gz",
                   DATA / "ont_candidate_sequences.fasta.gz", DATA / "coding_tiling.bam",
                   DATA / "ont_candidate_alignments.bam", DATA / "ont_exon_models.json")
    fusion, result = value["fusion"], value["result"]
    assert fusion.sequence[:3] == "ATG" and fusion.sequence[-3:] == "TGA"
    assert len(fusion.sequence) == 1596
    assert result["status"] == "translated"
    assert result["paths"][0]["junctions"][0]["direct_support"]["fragments"] == 149
    assert fusion.provenance["parameters"]["mapped_tiling_records"] == 19
    translation, = result["paths"][0]["translations"]
    assert translation["complete_5prime"] and translation["ends_with_stop_codon"]
    assert translation["acceptor_frames"] == [dict(transcript_id="ENST00000275016", in_frame=True)]
    assert len(translation["amino_acids"]) == 531
    assert sha256(translation["amino_acids"].encode()).hexdigest() == (
        "6d2fbb1872c08a507f56a977fb97328471b90a0f077211892d0a5d5ef3053573")
    assert PROTEIN in translation["amino_acids"]
    assert value["ont"] == dict(
        screened_original_reads=21, exact_junction_reads=6, flank_bases=20,
        full_coding_span_reads=9, full_coding_span_exact_junction_reads=4,
        full_coding_span_exact_junction_aligned_bases=[1513, 1539],
        coding_path_bases=1596)


def test_coding_replay_refuses_to_fill_a_missing_read_interval(tmp_path):
    source = DATA / "coding_tiling.bam"
    incomplete = tmp_path / "missing-one.bam"
    with pysam.AlignmentFile(str(source)) as bam:
        with pysam.AlignmentFile(str(incomplete), "wb", header=bam.header) as subset:
            for index, record in enumerate(bam.fetch(until_eof=True)):
                if index != 9:
                    subset.write(record)
    with pytest.raises(ValueError, match="read-overlap chain|unobserved base"):
        report(DATA / "original_unmapped.bam", DATA / "models.json.gz",
               mapped_tiling_bam=incomplete)


def test_coding_replay_rejects_a_misplaced_source_alignment(tmp_path):
    source = DATA / "coding_tiling.bam"
    misplaced = tmp_path / "misplaced.bam"
    with pysam.AlignmentFile(str(source)) as bam:
        with pysam.AlignmentFile(str(misplaced), "wb", header=bam.header) as subset:
            for index, record in enumerate(bam.fetch(until_eof=True)):
                if index == 9:
                    record.reference_start += 1
                subset.write(record)
    with pytest.raises(ValueError, match="disagrees with the pinned exon path"):
        report(DATA / "original_unmapped.bam", DATA / "models.json.gz",
               mapped_tiling_bam=misplaced)


def test_exon_anchored_sv_replay_keeps_original_split_support_and_frame():
    references = load_references(DATA / "models.json.gz")
    with pysam.AlignmentFile(str(DATA / "remapped_grch37_lite.bam")) as bam:
        result = reconstruct_sv_rna(
            bam, event_id="EIF3K--CYP39A1-RNA-adjacency", reference_name="GRCh37-lite",
            sample_id="HCC1395", source="INTEGRATE-2014-remapped-unmapped.bam",
            # These are observed RNA exon boundaries, not claimed DNA breakpoints.
            donor=FusionBreakpoint("19", 39123318, "+"),
            acceptor=FusionBreakpoint("6", 46607405, "-"), regions=[],
            references=references, event_provenance={"type": "observed_RNA_exon_adjacency",
                                                     "source": "INTEGRATE 2014 HCC1395 Illumina"},
            read_collector=ReadCollector(min_mapping_quality=0, use_soft_clipped_bases=True),
            max_records=2000, max_queries=300, max_paths=10,
            include_regional_candidates=False, exon_anchored_clips=True)
    assert result["status"] == "event_linked_candidates"
    assert not result["limitations"]
    assert any(seed["relation"] == "exon_anchored_clip_partner_unplaced"
               for seed in result["seeds"])
    supported = [path for path in result["paths"] if any(
        junction["relation"] == "breakpoint_junction" and
        junction["direct_support"]["fragments"] >= 120
        for junction in path["junctions"])]
    assert supported
    assert any(translation["amino_acids"] in PROTEIN and
               translation["amino_acids"].endswith("GTVNLH")
               for path in supported for translation in path["translations"])
