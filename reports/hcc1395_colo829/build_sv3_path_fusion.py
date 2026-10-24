"""Audit original read, mate, exon, and ORF evidence for COLO829 SV3."""

import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path

import pysam

from isovar.fusion import FusionBreakpoint, FusionReference
from isovar.read_collector import ReadCollector
from isovar.sv_rna import reconstruct_sv_rna
from isovar.sv_rna_path_support import (
    exact_path_matches, interval_spots, overlap_chain, reverse_complement)
from reports.hcc1395_colo829.build_sv_exon_path_evidence import (
    COLO_MODELS, COLO_SOURCE_BAM, COLO_SOURCE_READS, WINDOW43, _exon_anchors)
from reports.hcc1395_colo829.build_sv_junction_audit import EVENTS


HERE = Path(__file__).parent
OUTPUT = HERE / "sv3_path_fusion"
EVENT = EVENTS[3]


def _raw_ends(fastq):
    """Read the ten pinned original ends, preserving the SRA spot identity."""
    ends = {}
    with gzip.open(fastq, "rt") as stream:
        while name := stream.readline().strip():
            sequence = stream.readline().strip()
            plus = stream.readline().strip()
            quality = stream.readline().strip()
            if not name.startswith("@") or not plus.startswith("+") or len(sequence) != len(quality):
                raise ValueError("Malformed SV3 source FASTQ")
            spot, mate = name[1:].split()[0].rsplit(".", 1)
            if mate not in ("1", "2") or (spot, int(mate)) in ends:
                raise ValueError("Unexpected SV3 source spot or mate")
            ends[spot, int(mate)] = sequence.upper()
    return ends


def _bam_reads(bam):
    """Yield full original query sequences with available clip boundaries."""
    with pysam.AlignmentFile(str(bam), "rb") as alignment:
        for read in alignment.fetch(until_eof=True):
            if not read.query_sequence:
                continue
            mate = 1 if read.is_read1 else 2 if read.is_read2 else 0
            interval = (None, None) if read.is_unmapped else (
                read.query_alignment_start, read.query_alignment_end)
            yield read.query_name, mate, read.query_sequence, *interval


def _exact_exon_match(sequence, kmers):
    reverse = reverse_complement(sequence)
    return any(query[i:i + 20] in kmers for query in (sequence, reverse)
               for i in range(len(query) - 19))


def _alignment_records(bam, spots):
    """Retain every available placement for cited original spots."""
    records = defaultdict(list)
    with pysam.AlignmentFile(str(bam), "rb") as alignment:
        for read in alignment.fetch(until_eof=True):
            if read.query_name not in spots:
                continue
            records[read.query_name].append(dict(
                mate=1 if read.is_read1 else 2 if read.is_read2 else 0,
                contig=None if read.is_unmapped else read.reference_name,
                start=None if read.is_unmapped else read.reference_start,
                mapq=read.mapping_quality, cigar=read.cigarstring,
                secondary=read.is_secondary, supplementary=read.is_supplementary))
    return {spot: sorted(records[spot], key=lambda row: (
        row["mate"], row["contig"] or "", row["start"] or -1,
        row["secondary"], row["supplementary"])) for spot in sorted(spots)}


def build(bam=COLO_SOURCE_BAM, models=COLO_MODELS, fastq=COLO_SOURCE_READS,
          max_paths=80):
    """Replay the source fixture and assess two protein hypotheses separately."""
    model_data = json.loads(gzip.open(models, "rt").read())
    if model_data["schema"] != "isovar.colo829_sv3_reference_models.v1":
        raise ValueError("Unexpected versioned CNIH3 models")
    references = [FusionReference(**model) for model in model_data["models"]]
    sample, call_id, _, donor, acceptor, *_ = EVENT
    with pysam.AlignmentFile(str(bam), "rb") as alignment:
        result = reconstruct_sv_rna(
            alignment, event_id=call_id, reference_name="GRCh38", sample_id=sample,
            source="SRA:SRR8615617", donor=FusionBreakpoint(*donor),
            acceptor=FusionBreakpoint(*acceptor), regions=[], references=references,
            event_provenance={"dna_call_id": call_id},
            read_collector=ReadCollector(min_mapping_quality=0, use_soft_clipped_bases=True),
            max_records=20000, max_queries=200, max_paths=max_paths,
            max_extension_segments=1200, assemble=True,
            include_regional_candidates=False, dense_support=True,
            exon_anchored_clips=True)
    joined = [(path, junction) for path in result["paths"] for junction in path["junctions"]
              if junction["relation"] == "breakpoint_junction"]
    if not joined:
        raise ValueError("Pinned source BAM lost the original SV3 join")
    path, junction = max(joined, key=lambda pair: (
        pair[1]["direct_support"]["fragments"], len(pair[0]["sequence"])))
    raw = _raw_ends(fastq)
    reverse_window = reverse_complement(WINDOW43)
    window_ends = {key for key, sequence in raw.items()
                   if WINDOW43 in sequence or reverse_window in sequence}
    window_spots = sorted({spot for spot, _ in window_ends})
    placed_spots = sorted({result["observations"][key]["identity"][1]
                          for key in junction["direct_observations"]})
    if not set(placed_spots) <= set(window_spots):
        raise ValueError("Placed join has no exact raw junction-window sequence")

    reads = list(_bam_reads(bam))
    reads += [(spot, mate, sequence, None, None)
              for (spot, mate), sequence in raw.items()]
    matches = exact_path_matches(path["sequence"], reads)
    match_by_spot = defaultdict(list)
    for match in matches:
        if match.spot in window_spots:
            match_by_spot[match.spot].append(match)
    exon_kmers = {reference.sequence[i:i + 20] for reference in references
                  for i in range(len(reference.sequence) - 19)}
    source_spots = []
    for spot in window_spots:
        junction_mates = {mate for name, mate in window_ends if name == spot}
        opposite = {mate for mate in (1, 2) if mate not in junction_mates}
        source_spots.append(dict(
            spot=spot, exact_junction_mates=sorted(junction_mates),
            placed_split_join=spot in placed_spots,
            junction_end_has_clipped_alignment=any(
                match.mate in junction_mates and match.clipped
                for match in match_by_spot[spot]),
            opposite_mate_exact_path=any(
                match.mate in opposite and match.end - match.start >= 30
                for match in match_by_spot[spot]),
            opposite_mate_exact_cnih3_exon=any(
                _exact_exon_match(raw[spot, mate], exon_kmers) for mate in opposite)))

    anchor_bases, transcripts = _exon_anchors(path, references)
    left, right = junction["query_interval"]
    window_start, window_end = left - 12, right + 12
    candidates = []
    for orf in path["exploratory_orfs"]["candidates"]:
        if len(orf["amino_acids"]) not in (37, 162):
            continue
        start = orf["query_interval"][0]
        candidates.append(dict(
            amino_acids=len(orf["amino_acids"]),
            sequence=orf["amino_acids"], start_path_offset=start,
            start_codon=path["sequence"][start:start + 3],
            exact_start_spots=interval_spots(matches, start, start + 3),
            same_end_start_to_core_join_spots=interval_spots(matches, start, right),
            same_end_start_to_junction_spots=interval_spots(
                matches, start, window_end),
            start_to_junction_overlap_chain=overlap_chain(
                matches, start, window_end),
            full_orf_sequence_overlap_chain=overlap_chain(
                matches, start, orf["query_interval"][1], max_reads=12),
            annotated_start=orf["annotated_start"],
            initiation_observed=orf["initiation_observed"],
            frame_status=path["frame_status"]))
    if {row["amino_acids"] for row in candidates} != {37, 162}:
        raise ValueError("The two pinned SV3 exploratory ORFs changed")
    cited_spots = set(window_spots)
    for candidate in candidates:
        cited_spots.update(part["spot"] for key in (
            "start_to_junction_overlap_chain", "full_orf_sequence_overlap_chain")
                           for part in candidate[key])
    alignment_records = _alignment_records(bam, cited_spots)
    for spot in source_spots:
        spot["mapping_qualities"] = sorted({row["mapq"] for row in
                                            alignment_records[spot["spot"]]})
    return dict(
        schema="isovar.sv3_path_fusion.v1", event_id=call_id,
        source="SRA:SRR8615617", assembly="GRCh38", annotation="Ensembl 87",
        input_scope="CNIH3-target-screen candidate BAM plus ten original junction-spot reads",
        pacbio_original_long_reads="published_source_not_ingested",
        matched_normal_rna="not_evaluated",
        parameters=dict(max_records=20000, max_queries=200, max_paths=max_paths,
                        max_extension_segments=1200, min_exact_match=30,
                        min_overlap=30, max_orf_chain_reads=12,
                        exon_kmer=20, min_mapping_quality=0),
        input_records_scanned=result["support_acquisition"]["records_scanned"],
        path_count=len(result["paths"]), joined_path_count=len(joined),
        path_id=path["path_id"], path_bases=len(path["sequence"]),
        joined_path_sequence=path["sequence"],
        junction_path_interval=[left, right], junction_window_interval=[window_start, window_end],
        unplaced_junction_bases=junction["unplaced_bases"],
        exact_junction_window_spots=window_spots, placed_split_join_spots=placed_spots,
        source_spots=source_spots, cited_alignment_records=alignment_records,
        exact_exon_anchor_bases=anchor_bases, anchored_transcripts=transcripts,
        frame_status=path["frame_status"],
        event_linked_translation_candidates=sum(
            bool(set(translation["departure_relations"]) &
                 {"breakpoint_junction", "event_compatible_junction"})
            for translation in path["translations"]),
        candidates=sorted(candidates, key=lambda row: row["amino_acids"]),
        limitations=result["limitations"],
        conclusion="mature_path_and_frame_unresolved")


def write(prefix=OUTPUT, **kwargs):
    """Write a repeatable source-backed JSON assessment and short report."""
    report = build(**kwargs)
    prefix = Path(prefix)
    prefix.with_suffix(".json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    lines = [
        "# COLO829 SV3: fused original RNA evidence and protein hypotheses", "",
        "The source candidate BAM was replayed with dense paired-template support,",
        "soft-clipped bases, exon-anchored discovery, and no MAPQ threshold.",
        "The ten pinned original mates of the five exact junction spots were",
        "checked separately. The [machine-readable report](sv3_path_fusion.json)",
        "contains the assembled sequence, spot IDs, path intervals and parameters.", "",
        "Source: [SRR8615617](https://www.ncbi.nlm.nih.gov/sra/SRR8615617).",
        "Alignment and clipping terms follow the",
        "[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).", "",
        f"The retained joined path has **{report['path_bases']:,} bases**,",
        f"**{len(report['placed_split_join_spots'])} placed split-junction spots** and",
        f"**{len(report['exact_junction_window_spots'])} exact 43-base junction spots**.",
        f"It has **{report['exact_exon_anchor_bases']} exact CNIH3 exon-anchor bases**",
        f"and **{report['event_linked_translation_candidates']} event-linked frame translations**.",
        f"The bounded search retained {report['path_count']} paths and reported",
        f"{', '.join(report['limitations']) or 'no limits'}.", "",
        "| Spot | Exact-window end | Placed split | Junction-end match includes soft clip | Opposite mate on joined path | Opposite mate on CNIH3 exon | MAPQs |",
        "|---|---:|---|---|---|---|---|",
    ]
    for spot in report["source_spots"]:
        yn = lambda key: "yes" if spot[key] else "no"
        lines.append(f"| `{spot['spot']}` | {','.join(map(str, spot['exact_junction_mates']))} | "
                     f"{yn('placed_split_join')} | {yn('junction_end_has_clipped_alignment')} | "
                     f"{yn('opposite_mate_exact_path')} | {yn('opposite_mate_exact_cnih3_exon')} | "
                     f"{','.join(map(str, spot['mapping_qualities'])) or 'none'} |")
    lines += [
        "", "Mate sequence agreement with the joined path is a paired-fragment",
        "link. It does not fill the unsequenced insert or establish a splice",
        "path. Alternative placements and low MAPQ remain placement uncertainty;",
        "the exact 19-base unplaced join is still directly observed.", "",
        "| Exploratory ORF | Start-codon spot matches | One end start-to-core join | Start-to-43-mer overlap chain | Full ORF sequence tiling | Start/frame |",
        "|---|---:|---:|---|---|---|",
    ]
    for row in report["candidates"]:
        chain = row["start_to_junction_overlap_chain"]
        chain_text = " → ".join(f"`{part['spot']}`" for part in chain) if chain else "none"
        orf_chain = row["full_orf_sequence_overlap_chain"]
        orf_chain_text = f"{len(orf_chain)} distinct spots" if orf_chain else "none found"
        lines.append(f"| {row['amino_acids']} aa | {len(row['exact_start_spots'])} | "
                     f"{len(row['same_end_start_to_core_join_spots'])} | {chain_text} | "
                     f"{orf_chain_text} | unassigned |")
    lines += [
        "", "The two amino-acid sequences come from different ATGs and frames",
        "on the **same assembled nucleotide path**. The tilings use exact",
        "30-base-or-longer overlaps among distinct original spots; the JSON",
        "gives each interval. They support local RNA sequence content, not one",
        "observed full-length molecule or translation initiation.",
        "Start-codon spot matches each have at least 30 exact path bases",
        "covering the ATG. A single read spanning an entire ORF is not required.",
        "Neither start is annotated or experimentally observed. Their shared",
        "junction fragments cannot rank the proteins. There is no observed",
        "mature RNA connection to a versioned CNIH3 exon and no established",
        "frame, so neither exploratory ORF is an SV3 proteoform call.", "",
        "The input is an alignment-based target-screen fixture, not every RNA",
        "molecule or all possible splice paths. Its search hit the path limit.",
        "The RNA and nominated DNA were obtained from different COLO829",
        "cultures/passages. Matched-normal RNA was not evaluated.",
        "The [prior polished Iso-Seq check](colo829_sv3_polished_check.md)",
        "found no exact SV3 junction in representative sequences. Published",
        "[PacBio original COLO829 subreads](https://downloads-ap.pacbcloud.com/public/dataset/Melanoma2019_IsoSeq/subreads/COLO829T/)",
        "were not ingested: the direct PacBio host returned HTTP 403 and the",
        "Asia-Pacific host timed out on 2026-10-24. Original long-read RNA or an independent splice",
        "witness remains needed to resolve [#563](https://github.com/openvax/isovar/issues/563).", "",
        "Regenerate with `python -m reports.hcc1395_colo829.build_sv3_path_fusion`.",
    ]
    prefix.with_suffix(".md").write_text("\n".join(lines) + "\n")
    return report


def main():
    """Build the pinned SV3 path evidence report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, default=COLO_SOURCE_BAM)
    parser.add_argument("--models", type=Path, default=COLO_MODELS)
    parser.add_argument("--fastq", type=Path, default=COLO_SOURCE_READS)
    parser.add_argument("--output-prefix", type=Path, default=OUTPUT)
    parser.add_argument("--max-paths", type=int, default=80)
    args = parser.parse_args()
    write(prefix=args.output_prefix, bam=args.bam, models=args.models,
          fastq=args.fastq, max_paths=args.max_paths)


if __name__ == "__main__":
    main()
