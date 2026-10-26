"""Audit every pinned COLO829 CNIH3 candidate for an SV3 frame anchor."""

import argparse
from collections import defaultdict
import gzip
from hashlib import sha256
import json
from pathlib import Path

import pysam

from isovar.sv_rna_path_support import exact_path_matches, reverse_complement
from reports.hcc1395_colo829.build_sv3_path_fusion import _raw_ends
from reports.hcc1395_colo829.build_sv_exon_path_evidence import (
    COLO_MODELS, COLO_SOURCE_BAM, COLO_SOURCE_READS, WINDOW43)


HERE = Path(__file__).parent
PATH_REPORT = HERE / "sv3_path_fusion.json"
POLISHED_REPORT = HERE / "colo829_sv3_polished_check.json"
OUTPUT = HERE / "sv3_frame_anchor_audit"
MIN_PATH_MATCH = 20
MIN_EXON_OVERLAP = 15


def _sha256(path):
    digest = sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _coding_intervals(model):
    """Return genomic CDS blocks from a plus-strand transcript's cDNA offsets."""
    if model["cds_start"] is None:
        return []
    if model["strand"] != "+":
        raise ValueError("Pinned CNIH3 models must be on the plus strand")
    blocks, cdna_offset = [], 0
    for start, end in model["exons"]:
        left = max(model["cds_start"], cdna_offset)
        right = min(model["cds_end"], cdna_offset + end - start)
        if left < right:
            blocks.append((start + left - cdna_offset, start + right - cdna_offset))
        cdna_offset += end - start
    return blocks


def _catalog(models):
    exons = defaultdict(lambda: dict(transcripts=set(), coding_intervals=[]))
    coding_starts = []
    for model in models:
        coding = _coding_intervals(model)
        if coding:
            coding_starts.append(dict(transcript_id=model["transcript_id"],
                                      start_grch38_1based=coding[0][0] + 1))
        for exon in model["exons"]:
            item = exons[model["contig"], *exon]
            item["transcripts"].add(model["transcript_id"])
            item["coding_intervals"].extend(
                (a, b, model["transcript_id"]) for a, b in coding
                if min(b, exon[1]) > max(a, exon[0]))
    return exons, sorted(coding_starts, key=lambda row: row["transcript_id"])


def _clip_qualities(read):
    qualities = read.query_qualities
    cigar = read.cigartuples or []
    left = cigar[0][1] if cigar and cigar[0][0] == 4 else 0
    right = cigar[-1][1] if cigar and cigar[-1][0] == 4 else 0
    return dict(
        left_clip_bases=left, right_clip_bases=right,
        left_clip_min_q=min(qualities[:left]) if qualities is not None and left else None,
        right_clip_min_q=min(qualities[-right:]) if qualities is not None and right else None)


def _exon_placements(read, exons):
    if read.is_unmapped:
        return []
    placements = []
    for block_start, block_end in read.get_blocks():
        for (contig, exon_start, exon_end), info in exons.items():
            if contig != read.reference_name:
                continue
            overlap = min(block_end, exon_end) - max(block_start, exon_start)
            if overlap < MIN_EXON_OVERLAP:
                continue
            coding_ids = sorted({transcript for a, b, transcript in info["coding_intervals"]
                                 if min(block_end, b) - max(block_start, a) >= MIN_EXON_OVERLAP})
            placements.append(dict(
                exon_start_0based=exon_start, exon_end_0based=exon_end,
                aligned_exon_bases=overlap,
                transcript_ids=sorted(info["transcripts"]),
                coding_transcript_ids=coding_ids))
    return placements


def _path_genome_interval(match, blocks):
    for block in blocks:
        if block["query_start"] <= match.start and match.end <= block["query_end"]:
            if block["strand"] == "+":
                start = block["reference_start"] + match.start - block["query_start"]
                end = block["reference_start"] + match.end - block["query_start"]
            else:
                start = block["reference_end"] - (match.end - block["query_start"])
                end = block["reference_end"] - (match.start - block["query_start"])
            return [block["contig"], start, end]
    return None


def _relative_to_join(match, left, right):
    if match.end <= left:
        return "upstream_reference"
    if match.start >= right:
        return "downstream_reference"
    return "overlaps_unplaced_join"


def _aligned_exon_clip_contains(record, segment):
    """Require the exon placement's own query to carry a clipped match."""
    query = record["query_sequence"]
    if not query:
        return False
    start, end = record["query_aligned_interval"]
    for match_sequence in (segment, reverse_complement(segment)):
        offset = query.find(match_sequence)
        while offset >= 0:
            if offset < start or offset + len(match_sequence) > end:
                return True
            offset = query.find(match_sequence, offset + 1)
    return False


def build(bam=COLO_SOURCE_BAM, models=COLO_MODELS, fastq=COLO_SOURCE_READS,
          path_report=PATH_REPORT, polished_report=POLISHED_REPORT):
    """Count exon and exact-junction links without merging alternate RNA paths."""
    model_data = json.loads(gzip.open(models, "rt").read())
    path_data = json.loads(Path(path_report).read_text())
    polished = json.loads(Path(polished_report).read_text())
    if (model_data["schema"] != "isovar.colo829_sv3_reference_models.v1"
            or path_data["schema"] != "isovar.sv3_path_fusion.v1"
            or polished["schema"] != "isovar.colo829_sv3_polished_check.v1"
            or model_data["assembly"] != path_data["assembly"]
            or model_data["gene_id"] != "ENSG00000143786"):
        raise ValueError("Incompatible SV3 frame-anchor inputs")
    sequence = path_data["joined_path_sequence"]
    if len(sequence) != path_data["path_bases"] or not path_data.get("path_blocks"):
        raise ValueError("The joined path needs its exact placed blocks")
    left, right = path_data["junction_path_interval"]
    exons, coding_starts = _catalog(model_data["models"])
    reads, exon_by_spot, sequences_by_spot = [], defaultdict(list), defaultdict(set)
    records_scanned, bam_spots, coding_exon_spots = 0, set(), set()
    with pysam.AlignmentFile(str(bam), "rb") as alignment:
        for read in alignment.fetch(until_eof=True):
            records_scanned += 1
            spot = read.query_name
            bam_spots.add(spot)
            mate = 1 if read.is_read1 else 2 if read.is_read2 else 0
            if read.query_sequence:
                reads.append((spot, mate, read.query_sequence,
                              None if read.is_unmapped else read.query_alignment_start,
                              None if read.is_unmapped else read.query_alignment_end))
                sequences_by_spot[spot].add(read.query_sequence.upper())
            placements = _exon_placements(read, exons)
            if placements:
                if any(exon["coding_transcript_ids"] for exon in placements):
                    coding_exon_spots.add(spot)
                exon_by_spot[spot].append(dict(
                    mate=mate, reference_start_0based=read.reference_start,
                    reference_end_0based=read.reference_end,
                    cigar=read.cigarstring, mapq=read.mapping_quality,
                    secondary=read.is_secondary, supplementary=read.is_supplementary,
                    query_sequence=read.query_sequence,
                    query_aligned_interval=[read.query_alignment_start,
                                            read.query_alignment_end],
                    exons=placements, **_clip_qualities(read)))
    bam_exact_join_spots = sorted(spot for spot, sequences in sequences_by_spot.items()
                                  if any(WINDOW43 in seq or reverse_complement(WINDOW43) in seq
                                         for seq in sequences))
    raw = _raw_ends(fastq)
    for (spot, mate), original in raw.items():
        reads.append((spot, mate, original, None, None))
        sequences_by_spot[spot].add(original)
    reverse_window = reverse_complement(WINDOW43)
    exact_join_spots = sorted(spot for spot, sequences in sequences_by_spot.items()
                              if any(WINDOW43 in seq or reverse_window in seq
                                     for seq in sequences))
    if exact_join_spots != path_data["exact_junction_window_spots"]:
        raise ValueError("Exact SV3 junction spots differ from the pinned path audit")
    matches = exact_path_matches(sequence, reads, min_match=MIN_PATH_MATCH, seed_size=15)
    by_spot = defaultdict(list)
    path_match_spots = {match.spot for match in matches}
    for match in matches:
        if match.spot in exon_by_spot:
            by_spot[match.spot].append(match)
    links = []
    for spot in sorted(by_spot):
        placements = exon_by_spot[spot]
        candidates = []
        for match in by_spot[spot]:
            same_end = [record for record in placements if record["mate"] == match.mate]
            mate_end = [record for record in placements if record["mate"] != match.mate]
            segment = sequence[match.start:match.end]
            candidates.append(dict(
                mate=match.mate, path_interval=[match.start, match.end],
                exact_match_bases=match.end - match.start,
                genomic_interval_0based_half_open=_path_genome_interval(
                    match, path_data["path_blocks"]),
                path_region=_relative_to_join(match, left, right),
                distance_after_join_bases=max(0, match.start - right),
                uses_clipped_bases=match.clipped,
                same_end_exon_alignment=bool(same_end),
                same_end_exon_clip_contains_match=any(
                    _aligned_exon_clip_contains(record, segment) for record in same_end),
                opposite_mate_exon_alignment=bool(mate_end)))
        links.append(dict(
            spot=spot, exact_sv3_43mer=spot in exact_join_spots,
            matches=sorted(candidates, key=lambda row: (
                row["path_interval"], row["mate"])),
            exon_alignments=sorted(placements, key=lambda row: (
                row["mate"], row["reference_start_0based"], row["cigar"]))))
    linked_join_spots = sorted(set(by_spot) & set(exact_join_spots))
    joined_exon_spots = sorted(set(exon_by_spot) & set(exact_join_spots))
    polished_exact = polished["exact_polished_join_status"]
    has_frame = path_data["event_linked_translation_candidates"] > 0
    return dict(
        schema="isovar.sv3_frame_anchor_audit.v1",
        event_id=path_data["event_id"], source="SRA:SRR8615617",
        assembly=path_data["assembly"], annotation=path_data["annotation"],
        input_scope="alignment-based CNIH3 target-screen candidate BAM plus ten exact-junction raw ends",
        input_sha256=dict(bam=_sha256(bam), models=_sha256(models),
                          fastq=_sha256(fastq), path_report=_sha256(path_report),
                          polished_report=_sha256(polished_report)),
        parameters=dict(min_path_exact_bases=MIN_PATH_MATCH,
                        min_exon_aligned_bases=MIN_EXON_OVERLAP,
                        min_mapping_quality=0,
                        path_search_limit=path_data["parameters"]["max_paths"]),
        bam_alignment_records_scanned=records_scanned,
        bam_spots_scanned=len(bam_spots),
        bam_spots_with_versioned_exon_alignment=len(exon_by_spot),
        bam_spots_with_coding_exon_alignment=len(coding_exon_spots),
        spots_with_20_base_path_match=len(path_match_spots),
        raw_exact_junction_spots_missing_from_bam=sorted(set(exact_join_spots) - bam_spots),
        raw_exact_junction_spots_without_exact_bam_window=sorted(
            set(exact_join_spots) - set(bam_exact_join_spots)),
        joined_path_id=path_data["path_id"],
        joined_path_bases=path_data["path_bases"],
        joined_path_blocks=path_data["path_blocks"],
        junction_path_interval=[left, right],
        exact_sv3_43mer_spots=exact_join_spots,
        placed_sv3_split_spots=path_data["placed_split_join_spots"],
        exact_exon_anchor_bases_on_joined_path=path_data["exact_exon_anchor_bases"],
        joined_path_frame_status=path_data["frame_status"],
        event_linked_translation_candidates=path_data["event_linked_translation_candidates"],
        coding_reference_starts=coding_starts,
        exon_and_path_spot_count=len(links),
        exon_and_path_spots_with_30_base_match=sum(
            any(row["exact_match_bases"] >= 30 for row in link["matches"]) for link in links),
        exon_and_path_spots_with_cds_alignment=sum(
            any(exon["coding_transcript_ids"] for record in link["exon_alignments"]
                for exon in record["exons"]) for link in links),
        exact_sv3_join_and_exon_same_spot=linked_join_spots,
        exact_sv3_join_spots_with_exon_alignment=joined_exon_spots,
        exon_path_links=links,
        exploratory_orfs=[dict(amino_acids=row["amino_acids"],
                               start_path_offset=row["start_path_offset"],
                               exact_start_spots=len(row["exact_start_spots"]),
                               same_end_start_to_core_join_spots=len(
                                   row["same_end_start_to_core_join_spots"]),
                               same_end_start_to_junction_spots=len(
                                   row["same_end_start_to_junction_spots"]),
                               full_orf_sequence_overlap_chain_spots=len(
                                   row["full_orf_sequence_overlap_chain"]),
                               annotated_start=row["annotated_start"],
                               initiation_observed=row["initiation_observed"],
                               frame_status=row["frame_status"])
                          for row in path_data["candidates"]],
        polished_isoform_exact_join_status=polished_exact,
        polished_isoform_representatives=polished["representative_count"],
        polished_motif_positive_representatives=[
            row["isoform_id"] for row in polished["motif_positive_representatives"]],
        limitations=sorted(set(path_data["limitations"] + [
            "target_screen_not_all_original_reads", "no_matched_normal_rna",
            "original_pacbio_reads_not_ingested"])),
        conclusion=("event_linked_frame_candidate_in_path_report" if has_frame else
                    "no_event_linked_frame_anchor_in_scanned_inputs"))


def write(prefix=OUTPUT, **kwargs):
    """Write the machine-readable regional inventory and its interpretation."""
    report = build(**kwargs)
    prefix = Path(prefix)
    prefix.with_suffix(".json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    no_frame = report["conclusion"] == "no_event_linked_frame_anchor_in_scanned_inputs"
    all_downstream = all(match["path_region"] == "downstream_reference"
                         for link in report["exon_path_links"] for match in link["matches"])
    distances = {match["distance_after_join_bases"]
                 for link in report["exon_path_links"] for match in link["matches"]}
    lines = [
        "# COLO829 SV3: regional RNA evidence for a reading-frame anchor", "",
        ("The result is **no event-linked frame anchor in the scanned inputs**." if no_frame
         else "The path report contains an event-linked frame candidate requiring review."),
        "The original-read search found exon-linked sequences in the same",
        "downstream reference block as the SV3 RNA join, but none of those",
        "spots contains the SV3-specific 43-base junction. A shared reference",
        "segment does not phase two upstream RNA branches into one transcript.",
        "This advances [#563](https://github.com/openvax/isovar/issues/563).",
        "The exploratory ORFs remain separate from assigned proteins.", "",
        f"The [machine-readable audit]({prefix.name}.json) scanned",
        f"**{report['bam_alignment_records_scanned']:,} alignments** from",
        f"**{report['bam_spots_scanned']:,} SRA spots** in the [pinned CNIH3",
        "candidate BAM](../../tests/data/colo829_whole_run/SRR8615617.cnih3_source_candidates.grch38.bam),",
        "without a MAPQ threshold. Exact matches require at",
        "least 20 path bases; exon alignments require 15 aligned bases. The",
        "ten raw ends of the five known SV3 junction spots were also checked.",
        f"Across both sources, {report['spots_with_20_base_path_match']:,} spots",
        "have at least 20 exact bases of the selected joined path; most matches",
        "are to ordinary reference sequence away from the SV3 junction.",
        f"The BAM contains {report['bam_spots_with_versioned_exon_alignment']:,}",
        "spots with a qualifying Ensembl 87 exon alignment, including",
        f"{report['bam_spots_with_coding_exon_alignment']:,} with a coding-exon",
        "alignment. Coding-region expression alone does not link an RNA",
        "molecule to the SV3 join.", "",
        "## Exon-linked candidates", "",
        f"**{report['exon_and_path_spot_count']} spots** have both a versioned",
        "CNIH3 exon alignment and at least 20 exact bases of the selected SV3",
        f"path; **{report['exon_and_path_spots_with_30_base_match']}** reach 30 bases.",
        ("Every matched segment is in the ordinary downstream reference block,"
         if all_downstream else "Matched segments occupy multiple path regions; see the JSON."),
        (f"starting {next(iter(distances))} bases after the 19 unplaced SV3 join bases."
         if all_downstream and len(distances) == 1 else
         "Distances from the SV3 join vary across the observed matches."),
        "These are",
        "candidate alternative exon-to-region RNA links, not reads spanning",
        "the SV3-specific join.", "",
        "| SRA spot | Best exact path match, GRCh38 1-based | Exon end(s), 1-based | Link | Clip quality | SV3 43-mer |",
        "|---|---:|---|---|---|---|",
    ]
    for link in report["exon_path_links"]:
        best = max(link["matches"], key=lambda row: (
            row["same_end_exon_clip_contains_match"],
            row["exact_match_bases"]))
        same = [a for a in link["exon_alignments"] if a["mate"] == best["mate"]]
        ends = sorted({exon["exon_end_0based"] for a in same for exon in a["exons"]})
        qualities = sorted({a["right_clip_min_q"] for a in same
                            if a["right_clip_min_q"] is not None})
        kind = "same-end clipped" if best["same_end_exon_clip_contains_match"] else (
            "same-end" if best["same_end_exon_alignment"] else "mate only")
        interval = best["genomic_interval_0based_half_open"]
        location = (f"{interval[0]}:{interval[1] + 1}-{interval[2]}" if interval
                    else "unplaced")
        lines.append(
            f"| `{link['spot']}` | {best['exact_match_bases']} nt at {location} "
            f"(path `{best['path_interval'][0]}:{best['path_interval'][1]}`) | "
            f"{', '.join(map(str, ends)) or 'none'} | {kind} | "
            f"right-clip min Q {', '.join(map(str, qualities)) or 'unavailable'} | "
            f"{'yes' if link['exact_sv3_43mer'] else 'no'} |")
    starts = ", ".join(f"`{row['transcript_id']}` at chr1:{row['start_grch38_1based']:,}"
                       for row in report["coding_reference_starts"])
    lines += [
        "", "The exon coordinates above are **exon ends**, not the start of",
        "the downstream clipped segment. The JSON retains each original",
        "alignment, mate, CIGAR, MAPQ, exon model, clip quality, and exact",
        "path interval. Multiple alignments and overlapping transcript exons",
        "do not add independent spot support.",
        f"**{report['exon_and_path_spots_with_cds_alignment']}** exon-linked spots align to",
        "a coding block in Ensembl 87. A soft-clipped sequence match",
        "does not by itself establish a mature splice or reading frame.", "",
        "## Event-specific and frame evidence", "",
        f"The exact SV3 43-mer has **{len(report['exact_sv3_43mer_spots'])}**",
        f"distinct SRA spots; **{len(report['placed_sv3_split_spots'])}** have",
        f"placed split-junction support. **{len(report['exact_sv3_join_spots_with_exon_alignment'])}**",
        "have any qualifying exon alignment in the BAM, and none is among",
        "the exon-linked path matches.",
        f"{len(report['raw_exact_junction_spots_without_exact_bam_window'])} exact 43-mer",
        "spot requires the pinned raw FASTQ end because its sequence is not",
        "an exact 43-mer in the candidate BAM record.",
        f"The retained joined assembly has **{report['exact_exon_anchor_bases_on_joined_path']}**",
        "exact versioned exon-anchor bases and",
        f"**{report['event_linked_translation_candidates']}** event-linked",
        f"frame translations. Ensembl 87 coding CNIH3 start: {starts}.",
        "", "| Exploratory ORF | Original spots matching its ATG | Single ends covering ATG to join core / complete 43-mer | Full ORF overlap chain | Annotated or observed initiation |",
        "|---|---:|---:|---:|---|",
    ]
    for orf in report["exploratory_orfs"]:
        lines.append(f"| {orf['amino_acids']} aa | {orf['exact_start_spots']} | "
                     f"{orf['same_end_start_to_core_join_spots']} / "
                     f"{orf['same_end_start_to_junction_spots']} | "
                     f"{orf['full_orf_sequence_overlap_chain_spots']} spots | "
                     f"{'yes' if orf['annotated_start'] or orf['initiation_observed'] else 'no'} |")
    lines += [
        "", "An ATG sequence match is not evidence of translation initiation.",
        "An overlap chain combines different spots and does not establish",
        "one full-length molecule or a selected reading frame. The exact",
        "start-to-join spot IDs are in the [path-fusion audit](sv3_path_fusion.json).",
        "The [polished-isoform cross-check](colo829_sv3_polished_check.md)",
        "found no exact SV3 junction among",
        f"{report['polished_isoform_representatives']:,} representatives; its",
        "reverse-motif representatives are",
        ", ".join(f"`{name}`" for name in report["polished_motif_positive_representatives"]) + ".",
        "counts are a separate assay, not original-read support for this join.", "",
        "## Scope and next evidence", "",
        "This BAM retains candidates from an alignment-based screen of the",
        "full SRR8615617 source. Completely unaligned connectors can be missed.",
        f"The Isovar path search reached its {report['parameters']['path_search_limit']}-path limit.",
        "Matched-normal RNA",
        "and original COLO829 PacBio reads were not assessed. The DNA call and",
        "RNA source also came from different COLO829 cultures or passages.",
        "A supported path must connect the exact 19-base SV3 join to a",
        "versioned exon or other established ORF anchor, with orientation and",
        "frame evidence. Overlapping reads can establish that path; one",
        "full-length read is not required. The alternative exon-to-downstream",
        "links here are specific candidates for comparison, not a substituted",
        "SV3 proteoform.", "",
        "Regenerate offline with",
        "`python -m reports.hcc1395_colo829.build_sv3_frame_anchor_audit`.", "",
    ]
    prefix.with_suffix(".md").write_text("\n".join(lines))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, default=COLO_SOURCE_BAM)
    parser.add_argument("--models", type=Path, default=COLO_MODELS)
    parser.add_argument("--fastq", type=Path, default=COLO_SOURCE_READS)
    parser.add_argument("--path-report", type=Path, default=PATH_REPORT)
    parser.add_argument("--polished-report", type=Path, default=POLISHED_REPORT)
    parser.add_argument("--output-prefix", type=Path, default=OUTPUT)
    args = parser.parse_args()
    write(args.output_prefix, bam=args.bam, models=args.models,
          fastq=args.fastq, path_report=args.path_report,
          polished_report=args.polished_report)


if __name__ == "__main__":
    main()
