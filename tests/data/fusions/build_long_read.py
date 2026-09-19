"""Pin original Sid long-read records at nominated fusion breakpoints.

For each event, the breakpoint windows (+-1000) and every supplied reference
exon are fetched from the public, indexed BAM (``samtools view -M``). Segments
with records overlapping both breakpoint windows are kept unchanged: split and
spliced-through junction reads, never reads re-aligned or edited here.

Usage: python -m tests.data.fusions.build_long_read --output tests/data/fusions/long-read
"""
import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tempfile

import pysam

HERE = Path(__file__).parent
BUCKET = "https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/"
SOURCES = {
    "ONT-T1-tagged": "ONT/IPISRC044_ONT_upload/IPISRC044_ONT/processed/IPISRC044_T1_sclrs_ONT/"
                     "IPISRC044_T1_sclrs_ONT/IPISRC044_T1_sclrs_ONT.tagged.bam",
    "PacBio-T1": "pacbio/IPISRC044_T1_sclrs_live_pbmm2_mapped.bam",
}
# (event, source, input carrying breakpoints and reference models)
FIXTURES = [
    ("TPST1--CRCP", "PacBio-T1", "corpus/TPST1--CRCP-T1.input.json.gz"),
    ("FOXO3--STRADA-CCDC47", "ONT-T1-tagged", "corpus/FOXO3--STRADA-CCDC47-T1.input.json.gz"),
    ("ATP5MG--KMT2A", "PacBio-T1", "coding-corpus/ATP5MG--KMT2A.input.json.gz"),
]
WINDOW = 1000


def event_input(path):
    data = json.loads(gzip.decompress((HERE / path).read_bytes()))
    return data["fusion"]["donor"], data["fusion"]["acceptor"], data["references"]


def intervals(donor, acceptor, references):
    spans = [(b["contig"], max(0, b["position"] - WINDOW), b["position"] + WINDOW) for b in (donor, acceptor)]
    spans += [(r["contig"], a, b) for r in references for a, b in r["exons"]]
    merged = []
    for contig, start, end in sorted(spans):
        if merged and merged[-1][0] == contig and start <= merged[-1][2]:
            merged[-1][2] = max(end, merged[-1][2])
        else:
            merged.append([contig, start, end])
    return merged


def segment(read):
    return (read.get_tag("RG") if read.has_tag("RG") else "", read.query_name,
            read.flag & 0xC0 if read.is_paired else 0)


def build(event, source, input_path, output, scratch):
    donor, acceptor, references = event_input(input_path)
    url = BUCKET + SOURCES[source]
    regions = ["%s:%d-%d" % (c, a + 1, b) for c, a, b in intervals(donor, acceptor, references)]
    regional = scratch / ("%s.%s.bam" % (event, source))
    command = ["samtools", "view", "--no-PG", "-b", "-M", "-o", str(regional), "-X", url, url + ".bai"]
    subprocess.run(command + regions, check=True, cwd=scratch)  # htslib saves remote indexes in cwd.
    windows = [(b["contig"], b["position"] - WINDOW, b["position"] + WINDOW) for b in (donor, acceptor)]
    touches, records = {}, {}
    with pysam.AlignmentFile(str(regional)) as bam:
        header = str(bam.header)
        for read in bam:
            if read.is_unmapped:
                continue
            key = segment(read)
            records.setdefault(key, []).append(read.to_string())
            for k, (contig, start, end) in enumerate(windows):
                if read.reference_name == contig and read.reference_start < end and start < read.reference_end:
                    touches.setdefault(key, set()).add(k)
    kept = sorted(key for key, sides in touches.items() if sides == {0, 1})
    lines = sorted(line for key in kept for line in records[key])
    path = output / ("%s.%s.sam.gz" % (event, source))
    path.write_bytes(gzip.compress((header + "".join(line + "\n" for line in lines)).encode(), mtime=0))
    return dict(event=event, source=source, url=url, input=input_path, window=WINDOW, regions=len(regions),
                regional_records=sum(len(v) for v in records.values()), segments=len(kept), records=len(lines),
                command=" ".join(command[:5] + ["-o", "<output>", "-X", "<url>", "<url>.bai", "<regions>"]), file=path.name,
                sha256=sha256(path.read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as scratch:
        manifest = [build(*fixture, args.output, Path(scratch)) for fixture in FIXTURES]
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
