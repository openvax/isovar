"""Rebuild original long-RNA fixtures for three unresolved osteosarc events.

Use --local-corpus for a regional read-corpus directory containing
source_inventory.json and alignments/<source_id>/reads.bam, or --snapshot
for bounded acquisition through osteosarc. Ensembl 115 must already be cached.
See three-fusions/README.md for coordinate, selection and interpretation rules.
"""
import argparse
from collections import defaultdict
import gzip
from hashlib import sha256
import json
from pathlib import Path
import tempfile

import pysam
from pyensembl import EnsemblRelease

from isovar.sid_data import extract_regions, minimal_header, open_dataset, sam_digest

HERE = Path(__file__).parent
BUCKET = "https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/"
WINDOW = 2000
RELEASE = 115
EVENTS = {
    "PARD3B--CDKN2B-AS1-CDKN2B": {
        "donor": {"contig": "chr2", "position": 205310103, "strand": "+"},
        "acceptor": {"contig": "chr9", "position": 22007647, "strand": "+"},
        "genes": ["ENSG00000116117", "ENSG00000147883", "ENSG00000240498", "ENSG00000264545"],
    },
    "GABBR1--SLC29A1": {
        "donor": {"contig": "chr6", "position": 29612905, "strand": "-"},
        "acceptor": {"contig": "chr6", "position": 44218907, "strand": "+"},
        "genes": ["ENSG00000204681", "ENSG00000112759", "ENSG00000262179"],
    },
    "OTUD7A--FMN1": {
        "donor": {"contig": "chr15", "position": 31743669, "strand": "-"},
        "acceptor": {"contig": "chr15", "position": 33043478, "strand": "+"},
        "genes": ["ENSG00000169918", "ENSG00000248905"],
    },
}
SOURCES = {
    "ONT-T1": "53f498a544883d51", "ONT-T2": "1120a096937e29e1",
    "ONT-T3": "58c4d68e5f7e55b5", "PacBio-T1": "0066232879babe83",
}


def source_url(label):
    if label == "PacBio-T1":
        return BUCKET + "pacbio/IPISRC044_T1_sclrs_live_pbmm2_mapped.bam"
    sample = "IPISRC044_%s_sclrs_ONT" % label.split("-")[-1]
    return BUCKET + "ONT/IPISRC044_ONT_upload/IPISRC044_ONT/processed/%s/%s/%s.tagged.bam" % (sample, sample, sample)


def write_gzip(path, value):
    # GzipFile avoids the platform-specific OS byte of gzip.compress in Python 3.11/12.
    with path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, filename="", mode="wb", mtime=0) as handle:
            handle.write(json.dumps(value, sort_keys=True, separators=(",", ":")).encode())
    return sha256(path.read_bytes()).hexdigest()


def file_digest(path):
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def windows(event):
    return [(event[k]["contig"], event[k]["position"] - WINDOW,
             event[k]["position"] + WINDOW) for k in ("donor", "acceptor")]


def segment(read):
    return (read.get_tag("RG") if read.has_tag("RG") else "", read.query_name,
            read.flag & 0xC0 if read.is_paired else 0)


def select(bam, event):
    """Keep unedited local records for segments overlapping both end windows.

    Identical SAM lines retrieved in both queries are represented once. No
    sequence, quality, ORF, native-orientation or CB/UB selection is applied.
    This is a bounded path fixture, not an exhaustive event-support count.
    """
    sides, lines = defaultdict(set), {}
    for side, (contig, start, end) in enumerate(windows(event)):
        for read in bam.fetch(contig, start, end):
            if read.is_unmapped:
                continue
            key = segment(read)
            sides[key].add(side)
            lines[read.to_string()] = key
    selected = sorted(line for line, key in lines.items() if sides[key] == {0, 1})
    return selected, sum(value == {0, 1} for value in sides.values())


def references(genome, event):
    models = []
    for gene in sorted(event["genes"]):
        for tid in sorted(genome.transcript_ids_of_gene_id(gene)):
            t = genome.transcript_by_id(tid)
            sequence = t.sequence
            if sequence is None or set(sequence) - set("ACGT"):
                continue
            models.append(dict(
                transcript_id=t.id, reference_name="GRCh38", annotation="Ensembl %d" % RELEASE,
                contig="chr" + t.contig, strand=t.strand,
                exons=sorted((a - 1, b) for a, b in t.exon_intervals), sequence=sequence,
                cds_start=min(t.start_codon_spliced_offsets) if t.complete else None,
                cds_end=max(t.stop_codon_spliced_offsets) + 1 if t.complete else None))
    return models


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--local-corpus", type=Path)
    mode.add_argument("--snapshot")
    parser.add_argument("--cache")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    genome = EnsemblRelease(RELEASE)
    models = {name: references(genome, event) for name, event in EVENTS.items()}
    model_hash = write_gzip(args.output / "references.json.gz", models)
    manifest = dict(
        reference_name="GRCh38", annotation="Ensembl 115", events=EVENTS,
        annotation_sources=[genome.gtf_url] + genome.transcript_fasta_urls,
        references={"file": "references.json.gz", "sha256": model_hash},
        event_source="https://osteosarc.com/fusions/", window=WINDOW,
        coordinate_contract="0-based interbase retained donor/acceptor boundaries",
        selection="Original local SAM lines for RG/QNAME/mate segments touching both +/-2000 windows; identical lines once",
        software={"pysam": pysam.__version__, "samtools": pysam.__samtools_version__},
        sources={}, fixtures=[])
    if args.local_corpus:
        inventory_path = args.local_corpus / "source_inventory.json"
        inventory = {s["source_id"]: s for s in json.loads(inventory_path.read_text())}
        manifest["local_inventory_sha256"] = file_digest(inventory_path)
        inventory_file = "source-inventory.json.gz"
        manifest["source_inventory"] = dict(file=inventory_file, sha256=write_gzip(
            args.output / inventory_file, {sid: inventory[sid] for sid in SOURCES.values()}))
    dataset = open_dataset(args.snapshot, args.cache) if args.snapshot else None
    with tempfile.TemporaryDirectory() as scratch:
        for label, source_id in SOURCES.items():
            url = source_url(label)
            if args.local_corpus:
                entry = inventory[source_id]
                if BUCKET + entry["key"] != url:
                    raise ValueError("Unexpected source URL: " + source_id)
                path = args.local_corpus / "alignments" / source_id / "reads.bam"
                provenance = dict(source_id=source_id, url=url, local_bam_sha256=file_digest(path))
            else:
                path = Path(scratch) / (label + ".bam")
                regions = sorted({"%s:%d-%d" % (c, a + 1, b) for e in EVENTS.values() for c, a, b in windows(e)})
                subset = extract_regions(url, regions, "GRCh38", path, dataset=dataset)
                provenance = dict(source_id=source_id, url=url, acquisition=subset.receipt)
            manifest["sources"][label] = provenance
            with pysam.AlignmentFile(str(path)) as bam:
                for name, event in EVENTS.items():
                    lines, count = select(bam, event)
                    header = minimal_header(bam.header.to_dict(), lines)
                    if not lines:
                        contigs = {c for c, _, _ in windows(event)}
                        header["SQ"] = [row for row in bam.header.to_dict()["SQ"] if row["SN"] in contigs]
                    payload = dict(header=header,
                                   records={sam_digest(line): line for line in lines})
                    filename = name + "." + label + ".json.gz"
                    digest = write_gzip(args.output / filename, payload)
                    manifest["fixtures"].append(dict(
                        event=name, source=label, file=filename, sha256=digest,
                        segments=count, records=len(lines), regions=windows(event),
                        source_header_sha256=sha256(str(bam.header).encode()).hexdigest()))
                    print(name, label, count, "segments", len(lines), "records", flush=True)
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
