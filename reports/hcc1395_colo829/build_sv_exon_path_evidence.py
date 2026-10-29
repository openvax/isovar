"""Replay source SV joins with pinned Ensembl 87 exon and coding models."""

import csv
import gzip
import json
from pathlib import Path

import pysam

from isovar.fusion import FusionBreakpoint, FusionReference
from isovar.read_collector import ReadCollector
from isovar.sv_rna import reconstruct_sv_rna
from isovar.transcript_path_frames import TranscriptPathModel, exact_anchor_runs
from reports.hcc1395_colo829.build_sv_junction_audit import COLO_BAM, EVENTS, HCC_BAM


HERE = Path(__file__).parent
ROOT = HERE.parents[1]
HCC_MODELS = ROOT / "tests/data/hcc1395_sv_junction/ensembl87_models.json.gz"
COLO_MODELS = ROOT / "tests/data/colo829_whole_run/ensembl87_cnih3_models.json.gz"
COLO_SOURCE_BAM = ROOT / "tests/data/colo829_whole_run/SRR8615617.cnih3_source_candidates.grch38.bam"
COLO_SOURCE_READS = ROOT / "tests/data/colo829_whole_run/SRR8615617.sv3_exact_window_spots.fastq.gz"
SOURCE_AUDIT = HERE / "sv3_full_sra_screen.json"
SOURCE_FIELDS = ("sample", "call_id", "source_scope", "sra_spots_streamed",
                 "genome_remapped_cnih3_exon_spots", "path_id", "path_bases",
                 "placed_join_fragments", "exact_43mer_spots", "exact_exon_anchor_bases",
                 "frame_status", "event_linked_translation_candidates", "limitations")
WINDOW43 = "CTGTTGTCAGGCCTATAGTAGAACTTCTTTTGGAGTTGCTTCT"
PREFIX = HERE / "sv_exon_path_evidence"
GENES = ("SUGT1", "PKNOX1", "HNRNPK", "CNIH3")
FIELDS = ("sample", "call_id", "gene", "path_id", "path_bases", "direct_join_fragments",
          "local_exact_fragments", "exact_exon_anchor_bases", "anchored_transcripts",
          "exon_path_status", "frame_status", "event_linked_translation_candidates",
          "protein_status", "event_exon_clip_seeds", "limitations")


def _positions(path):
    """Restore per-base placements from the public path block representation."""
    positions = [None] * len(path["sequence"])
    for block in path["blocks"]:
        coordinates = range(block["reference_start"], block["reference_end"])
        if block["strand"] == "-":
            coordinates = reversed(coordinates)
        interval = range(block["query_start"], block["query_end"])
        if len(interval) != block["reference_end"] - block["reference_start"]:
            raise ValueError("SV path block has unequal query and reference lengths")
        for query, coordinate in zip(interval, coordinates):
            positions[query] = (block["contig"], coordinate, block["strand"])
    return tuple(positions)


def _exon_anchors(path, references):
    positions = _positions(path)
    anchors = [(reference.transcript_id, end - start)
               for reference in references
               for start, end, _ in exact_anchor_runs(path["sequence"], positions, TranscriptPathModel(reference))]
    maximum = max((length for _, length in anchors), default=0)
    transcripts = sorted({transcript for transcript, length in anchors if length >= 20})
    return maximum, transcripts


def build(hcc_bam=HCC_BAM, colo_bam=COLO_BAM, hcc_models=HCC_MODELS,
          colo_models=COLO_MODELS):
    """Replay direct joins with exon clip discovery and exact model anchors."""
    hcc = json.loads(gzip.open(hcc_models, "rt").read())
    colo = json.loads(gzip.open(colo_models, "rt").read())
    if (hcc["schema"] != "isovar.hcc1395_sv_reference_models.v1"
            or colo["schema"] != "isovar.colo829_sv3_reference_models.v1"):
        raise ValueError("Unexpected versioned SV reference model schema")
    rows = []
    for event, gene in zip(EVENTS[:4], GENES):
        sample, call_id, _, donor, acceptor, *_ = event
        bam = hcc_bam if sample == "HCC1395" else colo_bam
        models = hcc["models"][gene] if sample == "HCC1395" else colo["models"]
        references = [FusionReference(**model) for model in models]
        with pysam.AlignmentFile(str(bam), "rb") as alignment:
            result = reconstruct_sv_rna(
                alignment, event_id=call_id, reference_name="GRCh38", sample_id=sample,
                source=str(bam), donor=FusionBreakpoint(*donor),
                acceptor=FusionBreakpoint(*acceptor), regions=[], references=references,
                event_provenance={"dna_call_id": call_id},
                read_collector=ReadCollector(min_mapping_quality=0, use_soft_clipped_bases=True),
                max_records=5000, max_queries=100, max_paths=10, assemble=True,
                include_regional_candidates=False, dense_support=True, exon_anchored_clips=True)
        clip_seeds = sum(seed["relation"] == "exon_anchored_clip_partner_unplaced"
                         for seed in result["seeds"])
        for path in result["paths"]:
            joins = [junction for junction in path["junctions"]
                     if junction["relation"] == "breakpoint_junction"]
            if not joins:
                continue
            anchor_bases, transcripts = _exon_anchors(path, references)
            linked = [translation for translation in path["translations"]
                      if set(translation["departure_relations"]) &
                      {"breakpoint_junction", "event_compatible_junction"}]
            rows.append(dict(
                sample=sample, call_id=call_id, gene=gene, path_id=path["path_id"],
                path_bases=len(path["sequence"]),
                direct_join_fragments=max(j["direct_support"]["fragments"] for j in joins),
                local_exact_fragments=path["event_local_support"]["fragments"],
                exact_exon_anchor_bases=anchor_bases,
                anchored_transcripts=";".join(transcripts),
                exon_path_status=("exact_exon_anchor_on_joined_assembly" if transcripts
                                  else "no_exact_exon_anchor_on_joined_assembly"),
                frame_status=path["frame_status"],
                event_linked_translation_candidates=len(linked),
                protein_status=("reference_anchored_candidate" if linked
                                else "no_event_linked_protein_call"),
                event_exon_clip_seeds=clip_seeds,
                limitations=";".join(result["limitations"]) or "none"))
    return rows



def _exact_window_spots(fastq=COLO_SOURCE_READS):
    """Distinct original SRA spots carrying an exact strand-independent 43-mer."""
    lines = gzip.open(fastq, "rt").read().splitlines()
    if len(lines) % 4:
        raise ValueError("Incomplete FASTQ records in SV3 source window fixture")
    reverse = WINDOW43.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return sorted({lines[i].split()[0][1:].rsplit(".", 1)[0]
                   for i in range(0, len(lines), 4)
                   if WINDOW43 in lines[i + 1] or reverse in lines[i + 1]})


def build_source_extension(bam=COLO_SOURCE_BAM, models=COLO_MODELS,
                           fastq=COLO_SOURCE_READS, audit_path=SOURCE_AUDIT):
    """Replay the target-screened full SRA candidates beside the old fixture."""
    event = EVENTS[3]
    sample, call_id, _, donor, acceptor, *_ = event
    model_rows = json.loads(gzip.open(models, "rt").read())
    if model_rows["schema"] != "isovar.colo829_sv3_reference_models.v1":
        raise ValueError("Unexpected CNIH3 model schema")
    references = [FusionReference(**row) for row in model_rows["models"]]
    audit = json.loads(Path(audit_path).read_text())
    spots = _exact_window_spots(fastq)
    if (audit["schema"] != "isovar.sv3_full_sra_target_screen.v1"
            or audit["junction_window"]["sequence"] != WINDOW43
            or sorted(audit["junction_window"]["matching_spots_in_target_screen"]) != spots):
        raise ValueError("SV3 source screen audit differs from pinned original reads")
    with pysam.AlignmentFile(str(bam), "rb") as alignment:
        result = reconstruct_sv_rna(
            alignment, event_id=call_id, reference_name="GRCh38", sample_id=sample,
            source=str(bam), donor=FusionBreakpoint(*donor),
            acceptor=FusionBreakpoint(*acceptor), regions=[], references=references,
            event_provenance={"dna_call_id": call_id, "rna_accession": "SRR8615617"},
            read_collector=ReadCollector(min_mapping_quality=0, use_soft_clipped_bases=True),
            max_records=20000, max_queries=200, max_paths=20,
            max_extension_segments=1200, assemble=True,
            include_regional_candidates=False, dense_support=False,
            exon_anchored_clips=True)
    joined = [(path, junction) for path in result["paths"]
              for junction in path["junctions"] if junction["relation"] == "breakpoint_junction"]
    if not joined:
        raise ValueError("The full-source candidate fixture lost its observed SV3 join")
    path, junction = max(joined, key=lambda pair: (
        pair[1]["direct_support"]["fragments"], len(pair[0]["sequence"])))
    anchor_bases, _ = _exon_anchors(path, references)
    linked = [translation for translation in path["translations"]
              if set(translation["departure_relations"]) &
              {"breakpoint_junction", "event_compatible_junction"}]
    return dict(sample=sample, call_id=call_id,
                source_scope="all_SRR8615617_spots_target_screen_then_GRCh38_remap",
                sra_spots_streamed=audit["sra_spots_streamed"],
                genome_remapped_cnih3_exon_spots=audit["genome_remap"]["cnih3_exon_spots"],
                path_id=path["path_id"], path_bases=len(path["sequence"]),
                placed_join_fragments=junction["direct_support"]["fragments"],
                exact_43mer_spots=len(spots),
                exact_exon_anchor_bases=anchor_bases,
                frame_status=path["frame_status"],
                event_linked_translation_candidates=len(linked),
                limitations=";".join(result["limitations"]))


def write(prefix=PREFIX, **kwargs):
    """Export source exon path rows without asserting unobserved mature RNA."""
    rows = build(**kwargs)
    source = build_source_extension()
    prefix = Path(prefix)
    with prefix.with_suffix(".tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with prefix.with_name("sv3_full_sra_extension.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, SOURCE_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(source)
    lines = [
        "# Versioned exon path audit for selected SV joins", "",
        "Original HCC1395 and COLO829 selected-site BAMs were replayed with",
        "Ensembl 87 GRCh38 models, exon anchored clips, sequence overlaps,",
        "dense full-input support, and all alignment MAPQs. These are bounded",
        "**assembled path hypotheses**. The exon anchor column is the longest",
        "exact collinear match between a joined path and a versioned transcript.",
        "It does not alone prove mature splicing. A frame candidate is an RNA",
        "sequence prediction, not protein expression or somatic causation.",
        "Counts here belong to each retained path. Enabling exon clip seeds",
        "changes the bounded path set, so direct counts here can differ from",
        "the original join table and must not be summed across paths.", "",
        "| Event | Path | Joined RNA bases | Direct join | Local exact | Longest exon anchor | Frame | Event-linked translations | Limits |",
        "|---|---|---:|---:|---:|---:|---|---:|---|",
    ]
    for row in rows:
        lines.append("| %s %s | `%s` | %d | %d | %d | %d | %s | %d | %s |" % (
            row["sample"], row["call_id"], row["path_id"][:8], row["path_bases"],
            row["direct_join_fragments"], row["local_exact_fragments"],
            row["exact_exon_anchor_bases"], row["frame_status"],
            row["event_linked_translation_candidates"],
            row["limitations"].replace(";", ", ") or "none"))
    lines += [
        "", "## SV3 full SRA target screen", "",
        f"All {source['sra_spots_streamed']:,} SRR8615617 spots were streamed against the",
        "GRCh38 CNIH3 locus plus the prior 1,213-base SV3 path. Candidate",
        "reads were remapped against whole GRCh38. Of those,",
        f"{source['genome_remapped_cnih3_exon_spots']:,} distinct spots mapped to Ensembl 87",
        "CNIH3 exons. The committed 4 MB",
        "CNIH3 candidate BAM retains actual original selected-site alignments",
        "and full-source reads outside those windows, including alternative",
        "placements. This is an alignment-based target screen, so completely",
        "unaligned connectors may be missed. The exact-window count is a",
        "lower bound on distinct SRA spots, not a count of independent",
        "molecules. These reads are from a different culture or passage",
        "than the nominated DNA call.", "",
        "| Screened source | Joined path | RNA bases | Placed join fragments | Exact 43-mer spots | Longest exact CNIH3 exon anchor | Frame | Limits |",
        "|---|---|---:|---:|---:|---:|---|---|",
        "| SRR8615617 CNIH3 target screen | `%s` | %d | %d | %d | %d | %s | %s |" % (
            source["path_id"][:8], source["path_bases"], source["placed_join_fragments"],
            source["exact_43mer_spots"], source["exact_exon_anchor_bases"],
            source["frame_status"], source["limitations"].replace(";", ", ") or "none"),
        "", "The five exact 43-mer spots include the three placed split",
        "fragments and two with clipped or partial placements. The latter",
        "add sequence evidence without adding placed breakpoint support.",
        "Their ten original 101-base mates and qualities are pinned in",
        "[the source FASTQ](../../tests/data/colo829_whole_run/SRR8615617.sv3_exact_window_spots.fastq.gz).",
        "None of those ten reads has an exact 20-base CNIH3 transcript",
        "match. The expanded joined path also has no exact versioned exon",
        "anchor or transferred frame. The 37-aa and 162-aa ATG-to-stop",
        "sequences remain exploratory ORFs sharing the same nucleotide",
        "witnesses; their starts and frames are unestablished. No single",
        "read must span the entire assembled path. The old 1,213-base",
        "full-interval count of zero is only a strict one-read measure.",
        "Low MAPQ marks uncertain placement without erasing join sequence.",
        "The 20-path search hit `path_limit`, so unretained alternatives",
        "remain unevaluated. See the [screen provenance](sv3_full_sra_screen.json),",
        "[machine-readable source row](sv3_full_sra_extension.tsv) and",
        "[candidate BAM](../../tests/data/colo829_whole_run/SRR8615617.cnih3_source_candidates.grch38.bam).",
        "The [full SRA screen methods](SV3_FULL_SRA_SCREEN.md) document",
        "source retrieval, global remapping and candidate selection.",
        "The missing mature path is tracked in",
        "[Isovar #563](https://github.com/openvax/isovar/issues/563);",
        "dense paired-template recount scaling is tracked in",
        "[Isovar #562](https://github.com/openvax/isovar/issues/562).", "",
        "Regenerate both report tables with",
        "`python -m reports.hcc1395_colo829.build_sv_exon_path_evidence`.",
        "The model fixture is generated by",
        "`python -m reports.hcc1395_colo829.build_colo829_cnih3_models`.", "",
    ]
    prefix.with_suffix(".md").write_text("\n".join(lines))
    return rows


if __name__ == "__main__":
    write()
