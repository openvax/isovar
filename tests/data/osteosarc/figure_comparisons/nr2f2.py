"""Reproducible matched CeGaT DNA/RNA audit of the nonfocal NR2F2 deletion.

Acquisition is explicit, never performed by tests. Supply complete regional
BAMs named normal.bam, tumor.bam and rna.bam plus hg19-reference.json from
the URLs below; original records and hashes are retained in the offline pin.
"""

import argparse
from collections import Counter, defaultdict
import gzip
from hashlib import sha256
import json
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[4]
CORPUS = Path(__file__).parent / "corpus"
BUCKET = "https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/"
REFERENCE_URL = "https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr15;start=96875510;end=96875600"
REGION = "chr15:96875200-96875850"
SOURCES = {
    "normal": ("Blood DNA", "vendor/cegat/P116686_1_S000048/P116686_1.bam"),
    "tumor": ("Tumor DNA", "vendor/cegat/P116686_2_S000048/P116686_2.bam"),
    "rna": ("Tumor RNA", "vendor/cegat/P116686_2_S000048/P116686_3.bam"),
}
FOCAL = 96875527
DELETION = (96875576, 96875579)
EXCLUDED = 4 | 256 | 512 | 1024 | 2048


def recount(source, reference, min_quality=20, min_mapq=20, flank=8):
    """Exact anchored alleles, fragment counts and same-segment phase evidence."""
    header = pysam.AlignmentHeader.from_text(source["header"])
    lo, hi = DELETION[0] - flank, DELETION[1] + flank
    ref = reference["dna"][lo-reference["start"]:hi-reference["start"]].upper()
    assert ref[flank:flank+3] == "GTG"
    hypotheses = {ref: "reference", ref[:flank] + ref[flank+3:]: "deletion"}
    states = {name: defaultdict(set) for name in ("focal", "deletion", "joint")}
    deletion_records = []
    for sam in source["records"]:
        read = pysam.AlignedSegment.fromstring(sam, header)
        if read.flag & EXCLUDED or read.mapping_quality < min_mapq or read.query_qualities is None:
            continue
        coordinates = {g: q for q, g in read.get_aligned_pairs(matches_only=True)}
        key = (read.get_tag("RG") if read.has_tag("RG") else "", read.query_name)
        q = coordinates.get(FOCAL)
        focal = None
        if q is not None and read.query_qualities[q] >= min_quality:
            base = read.query_sequence[q]
            focal = {"G": "reference", "T": "alternate"}.get(base, "other")
            states["focal"][key].add(focal)
        left, right = coordinates.get(lo), coordinates.get(hi - 1)
        if left is None or right is None or min(read.query_qualities[left:right+1]) < min_quality:
            continue
        observed = read.query_sequence[left:right+1]
        deletion = hypotheses.get(observed, "other")
        states["deletion"][key].add(deletion)
        if focal is not None:
            states["joint"][key].add(focal + "/" + deletion)
        if deletion == "deletion":
            deletion_records.append(dict(
                name=read.query_name, flag=read.flag, mapq=read.mapping_quality,
                start=read.reference_start, end=read.reference_end, cigar=read.cigarstring,
                reverse=read.is_reverse, mate_start=read.next_reference_start,
                query_start=left, query_end=right+1, sequence=observed,
                minimum_quality=min(read.query_qualities[left:right+1]), focal=focal,
                sam_sha256=sha256(sam.encode()).hexdigest()))
    counts = {name: dict(sorted(Counter(
        next(iter(values)) if len(values) == 1 else "conflicting_template"
        for values in observations.values()).items())) for name, observations in states.items()}
    return dict(counts=counts, reference_window=ref,
                deletion_window=ref[:flank] + ref[flank+3:], window=[lo, hi],
                deletion_records=deletion_records,
                deletion_alignment_paths=len({(r["start"], r["cigar"], r["reverse"]) for r in deletion_records}),
                deletion_fragment_endpoints=len({(r["start"], r["end"], r["mate_start"]) for r in deletion_records}))


def build(inputs, output):
    reference = json.loads((inputs / "hg19-reference.json").read_text())
    assert (reference["genome"], reference["chrom"], reference["start"], reference["end"]) == (
        "hg19", "chr15", 96875510, 96875600)
    data = dict(reference=reference, reference_url=REFERENCE_URL, region=REGION,
                focal=dict(position=FOCAL+1, ref="G", alt="T"),
                deletion=dict(interval=list(DELETION), ref="GTG", alt=""),
                metadata_url="https://osteosarc.com/bams/bams.json",
                sample="CeGaT T0, tissue collected 2022-12-16; sequenced 2023-03-22",
                excluded_flags=EXCLUDED, sources={})
    for key, (label, path) in SOURCES.items():
        bam_path = inputs / (key + ".bam")
        with pysam.AlignmentFile(bam_path) as bam:
            assert bam.get_reference_length("chr15") == 102531392
            source = dict(label=label, url=BUCKET+path, index_url=BUCKET+path[:-4]+".bai",
                          regional_bam_sha256=sha256(bam_path.read_bytes()).hexdigest(),
                          header=str(bam.header), records=[r.to_string() for r in bam])
        source["audits"] = {
            "q20": recount(source, reference),
            "q30": recount(source, reference, min_quality=30),
            "mapq1": recount(source, reference, min_mapq=1),
        }
        data["sources"][key] = source
        print(key, {name: audit["counts"] for name, audit in source["audits"].items()})
    # Independent pinned transcript reference must agree with the entire
    # hg19 genomic window, not just the three bases nominated for deletion.
    from tests.data.osteosarc.expansion.references import load_reference
    from tests.osteosarc_protein_helpers import transcript_offset
    directory = ROOT / "tests/data/osteosarc/expansion/corpus/references/GRCh37-cegat"
    _, models = load_reference(directory)
    model = models["ENST00000394166"]
    start = transcript_offset(model["exons"], model["strand"], reference["start"]+1)
    assert model["cdna"][start:start+len(reference["dna"])] == reference["dna"].upper()
    data["reference_transcript"] = model["transcript_version"]
    raw = gzip.compress(json.dumps(data, sort_keys=True).encode(), mtime=0)
    with output.open("xb") as handle:
        handle.write(raw)
    manifest = dict(file=output.name, sha256=sha256(raw).hexdigest())
    with output.with_name("nr2f2-manifest.json").open("x") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")


def load():
    manifest = json.loads((CORPUS / "nr2f2-manifest.json").read_text())
    raw = (CORPUS / manifest["file"]).read_bytes()
    if sha256(raw).hexdigest() != manifest["sha256"]:
        raise ValueError("NR2F2 evidence checksum mismatch")
    return json.loads(gzip.decompress(raw))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build(args.inputs, args.output)
