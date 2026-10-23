"""Regenerate the pinned HCC1395 EIF3K--CYP39A1 source-read fixture.

Requires the extracted public INTEGRATE test-data.tar.gz archive. Optionally
accepts the locally realigned HCC1395 ONT BAM for its independent read subset.
"""

import argparse
from dataclasses import asdict
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

import pysam
from pyensembl import EnsemblRelease

from isovar.fusion import FusionReference
from build_hcc1395_eif3k_fusion import (
    ACCEPTOR_EDGE, ACCEPTOR_ID, DONOR_EDGE, DONOR_ID,
)


def build_models(source_dir):
    """Use GRCh37 exons and only CDS offsets whose cDNA matches Ensembl 87."""
    fasta = pysam.FastaFile(str(source_dir / "reference.fasta"))
    release = EnsemblRelease(87)
    models = []
    for line in (source_dir / "annot.ensembl.test.txt").read_text().splitlines():
        if line.startswith("#"):
            continue
        transcript_id, chromosome, strand, _, _, _, starts, ends, gene = line.split("\t")
        if gene not in ("EIF3K", "CYP39A1"):
            continue
        exons = tuple(zip((int(x) for x in starts.strip(",").split(",")),
                           (int(x) for x in ends.strip(",").split(","))))
        sequence = "".join(fasta.fetch(chromosome[3:], start, end) for start, end in exons)
        if strand == "-":
            sequence = pysam.reverse_complement(sequence)
        fields = (transcript_id, "GRCh37-lite", "INTEGRATE 2014 hg19 Ensembl exons",
                  chromosome[3:], strand, exons, sequence)
        try:
            transcript = release.transcript_by_id(transcript_id)
        except ValueError:
            transcript = None
        if transcript is not None and transcript.sequence.upper() != sequence:
            raise ValueError("Ensembl cDNA changed for " + transcript_id)
        complete = (transcript is not None and transcript.contains_start_codon and
                    transcript.contains_stop_codon and transcript.start_codon_complete and
                    transcript.stop_codon_complete)
        start = transcript.first_start_codon_spliced_offset if complete else None
        end = transcript.last_stop_codon_spliced_offset + 1 if complete else None
        try:
            model = FusionReference(*fields, start, end)
        except ValueError:
            model = FusionReference(*fields)
        models.append(model)
    if not {DONOR_ID, ACCEPTOR_ID} <= {model.transcript_id for model in models}:
        raise ValueError("Missing canonical fusion models")
    return models


def _one_placement(sequence, full):
    """Find a unique exact forward or reverse-complement transcript placement."""
    hits = []
    for oriented in (sequence, pysam.reverse_complement(sequence)):
        start = full.find(oriented)
        while start >= 0:
            hits.append((start, oriented))
            start = full.find(oriented, start + 1)
    return hits[0] if len(hits) == 1 else None


def select_coding_tiling(source_dir, output, models):
    """Keep original mapped records that connect the coding path by 20-base overlaps."""
    donor = next(model for model in models if model.transcript_id == DONOR_ID)
    acceptor = next(model for model in models if model.transcript_id == ACCEPTOR_ID)
    junction = donor.spliced_offset(DONOR_EDGE - 1) + 1
    acceptor_offset = acceptor.spliced_offset(ACCEPTOR_EDGE - 1)
    full = donor.sequence[:junction] + acceptor.sequence[acceptor_offset:]
    lo, hi = donor.cds_start, junction + acceptor.cds_end - acceptor_offset
    candidates = []
    with pysam.AlignmentFile(str(source_dir / "unmapped.bam"), check_sq=False) as bam:
        for record in bam.fetch(until_eof=True):
            placement = _one_placement(record.query_sequence, full)
            if placement is None:
                continue
            start, oriented = placement
            if min(junction - start, start + len(oriented) - junction) >= 10:
                candidates.append((start, start + len(oriented), None, record.query_name))
    with pysam.AlignmentFile(str(source_dir / "accepted_hits.bam")) as bam:
        for chrom, start, end in (("19", 39109736, DONOR_EDGE),
                                  ("6", 46517540, ACCEPTOR_EDGE)):
            for record in bam.fetch(chrom, start, end):
                if record.is_unmapped or record.mapping_quality < 20 or not record.query_sequence:
                    continue
                placement = _one_placement(record.query_sequence, full)
                if placement is None:
                    continue
                offset, oriented = placement
                if not ((offset + len(oriented) <= junction and chrom == "19") or
                        (offset >= junction and chrom == "6")):
                    continue
                if offset >= hi or offset + len(oriented) <= lo:
                    continue
                key = (record.query_name, record.flag, record.reference_start)
                candidates.append((offset, offset + len(oriented), key, record.query_name))
    reach, chosen = lo, set()
    while reach < hi:
        options = (row for row in candidates
                   if row[0] <= max(lo, reach - 20) and row[1] > reach)
        best = max(options, key=lambda row: (row[1], row[2] is None, row[3]), default=None)
        if best is None:
            raise ValueError("The original reads do not tile the complete fusion CDS")
        if best[2] is not None:
            chosen.add(best[2])
        reach = best[1]
    with pysam.AlignmentFile(str(source_dir / "accepted_hits.bam")) as bam:
        with pysam.AlignmentFile(str(output), "wb", header=bam.header) as subset:
            for record in bam.fetch(until_eof=True):
                if (record.query_name, record.flag, record.reference_start) in chosen:
                    subset.write(record)
    return len(chosen)


def remap_unmapped(source_dir, output):
    """Recreate the small, exon-anchored short-read BAM with pinned commands."""
    with tempfile.TemporaryDirectory() as scratch:
        fastq = Path(scratch) / "original.fastq"
        sam = Path(scratch) / "remapped.sam"
        with pysam.AlignmentFile(str(source_dir / "unmapped.bam"), check_sq=False) as bam:
            with fastq.open("w") as handle:
                for record in bam.fetch(until_eof=True):
                    handle.write("@%s\n%s\n+\n%s\n" % (
                        record.query_name, record.query_sequence,
                        record.qual or "I" * len(record.query_sequence)))
        with sam.open("w") as handle:
            subprocess.run(["minimap2", "-t", "4", "-ax", "sr", "-k11", "-w5", "-Y",
                            str(source_dir / "reference.fasta"), str(fastq)],
                           stdout=handle, check=True)
        subprocess.run(["samtools", "sort", "-o", str(output), str(sam)], check=True)
        subprocess.run(["samtools", "index", str(output)], check=True)


def select_ont_candidates(bam_path, output):
    """Pin original SRA query sequences with local EIF3K and CYP39A1 placements."""
    names = {}
    with pysam.AlignmentFile(str(bam_path)) as bam:
        for chrom, start, end in (("chr19", 38617000, 38639000),
                                  ("chr6", 46548000, 46655000)):
            for record in bam.fetch(chrom, start, end):
                if record.is_unmapped:
                    continue
                parts = names.setdefault(record.query_name, {})
                old = parts.get(chrom)
                if old is None or len(record.query_sequence or "") > len(old):
                    parts[chrom] = record.query_sequence or ""
    rows = [">%s\n%s\n" % (name, max(parts.values(), key=len))
            for name, parts in sorted(names.items()) if len(parts) == 2]
    with output.open("wb") as raw, gzip.GzipFile(filename="", fileobj=raw,
                                                   mode="wb", mtime=0) as compressed:
        compressed.write("".join(rows).encode())
    return len(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--ont-bam", type=Path)
    args = parser.parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    models = build_models(args.source_dir)
    data = dict(schema="isovar.hcc1395_eif3k_cyp39a1_models.v1",
                reference_name="GRCh37-lite",
                exon_source="INTEGRATE 2014 annot.ensembl.test.txt and reference.fasta",
                cds_source="Ensembl 87 GRCh38 for exact-sequence-matched transcripts",
                models=[asdict(model) for model in models])
    with (args.output_dir / "models.json.gz").open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            compressed.write(json.dumps(data, sort_keys=True, separators=(",", ":")).encode())
    shutil.copyfile(args.source_dir / "unmapped.bam", args.output_dir / "original_unmapped.bam")
    count = select_coding_tiling(args.source_dir, args.output_dir / "coding_tiling.bam", models)
    remap_unmapped(args.source_dir, args.output_dir / "remapped_grch37_lite.bam")
    print("Selected %d original mapped reads for the coding tile" % count)
    if args.ont_bam is not None:
        print("Selected %d original ONT queries" % select_ont_candidates(
            args.ont_bam, args.output_dir / "ont_candidate_sequences.fasta.gz"))


if __name__ == "__main__":
    main()
