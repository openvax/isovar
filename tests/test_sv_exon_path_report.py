"""Source-derived exon and frame checks for selected HCC1395/COLO829 SVs."""

import gzip
import json
from pathlib import Path

from reports.hcc1395_colo829.build_sv_exon_path_evidence import (
    _positions, build, build_source_extension)


DATA = Path(__file__).parent / "data/colo829_whole_run"


def test_cnih3_reference_models_are_pinned_without_live_pyensembl():
    models = json.loads(gzip.open(DATA / "ensembl87_cnih3_models.json.gz", "rt").read())
    assert models["schema"] == "isovar.colo829_sv3_reference_models.v1"
    assert models["assembly"] == "GRCh38" and models["release"] == 87
    assert len(models["models"]) == 14
    assert {model["transcript_id"] for model in models["models"]} >= {"ENST00000272133"}
    assert not any(start <= 224595092 < end or start <= 224598331 < end
                   for model in models["models"] for start, end in model["exons"])


def test_path_block_restoration_keeps_reverse_strand_orientation():
    path = {"sequence": "ACGTA", "blocks": [{"query_start": 1, "query_end": 4,
            "contig": "chr9", "reference_start": 100, "reference_end": 103,
            "strand": "-"}]}
    assert _positions(path) == (None, ("chr9", 102, "-"), ("chr9", 101, "-"),
                                ("chr9", 100, "-"), None)


def test_original_sv3_join_has_no_exon_link_or_transferred_frame():
    rows = build()
    assert len(rows) == 6
    sv3, = [row for row in rows if row["call_id"] == "COLO829_SOMATIC_SV3_DEL"]
    assert sv3["path_id"] == "88ab60e1f0f441a13cbecfe6"
    assert sv3["direct_join_fragments"] == sv3["local_exact_fragments"] == 3
    assert sv3["exact_exon_anchor_bases"] == 0
    assert sv3["anchored_transcripts"] == ""
    assert sv3["exon_path_status"] == "no_exact_exon_anchor_on_joined_assembly"
    assert sv3["frame_status"] == "no_gene_anchor"
    assert sv3["event_linked_translation_candidates"] == 0
    assert sv3["protein_status"] == "no_event_linked_protein_call"
    assert all(row["exact_exon_anchor_bases"] >= 20 for row in rows if row["sample"] == "HCC1395")
    assert not any(row["event_linked_translation_candidates"] for row in rows)


def test_sv3_full_sra_window_spots_have_no_direct_exon_anchor():
    lines = gzip.open(DATA / "SRR8615617.sv3_exact_window_spots.fastq.gz", "rt").read().splitlines()
    assert len(lines) == 40
    reads = [(lines[i].split()[0][1:], lines[i + 1]) for i in range(0, len(lines), 4)]
    spots = {name.rsplit(".", 1)[0] for name, _ in reads}
    assert spots == {"SRR8615617.7312123", "SRR8615617.7313798",
                     "SRR8615617.97504216", "SRR8615617.101713327",
                     "SRR8615617.110737497"}
    assert all({name.rsplit(".", 1)[1] for name, _ in reads
                if name.startswith(spot + ".")} == {"1", "2"} for spot in spots)
    window = "CTGTTGTCAGGCCTATAGTAGAACTTCTTTTGGAGTTGCTTCT"
    reverse = window.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    assert {name.rsplit(".", 1)[0] for name, sequence in reads
            if window in sequence or reverse in sequence} == spots
    models = json.loads(gzip.open(DATA / "ensembl87_cnih3_models.json.gz", "rt").read())
    transcripts = [row["sequence"] for row in models["models"]]
    for _, sequence in reads:
        reverse_sequence = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        assert not any(query[start:start + 20] in transcript
                       for query in (sequence, reverse_sequence)
                       for start in range(len(query) - 19)
                       for transcript in transcripts)


def test_full_sra_target_screen_extends_sv3_without_inventing_a_frame():
    row = build_source_extension()
    assert row["path_id"] == "c2162971d05bbc0d5de89d46"
    assert row["sra_spots_streamed"] == 115269780
    assert row["genome_remapped_cnih3_exon_spots"] == 6291
    assert row["path_bases"] == 2435
    assert row["placed_join_fragments"] == 3
    assert row["exact_43mer_spots"] == 5
    assert row["exact_exon_anchor_bases"] == 0
    assert row["frame_status"] == "no_gene_anchor"
    assert row["event_linked_translation_candidates"] == 0
    assert "path_limit" in row["limitations"]
