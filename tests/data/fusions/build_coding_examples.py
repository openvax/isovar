"""Pin real coding hypotheses without discarding incompatible CDS alternatives."""
import argparse
from collections import defaultdict
import gzip
from hashlib import sha256
import json
from pathlib import Path

import pysam
from pyensembl import EnsemblRelease

from .build_osteosarc import blocks, extract, references

SID_URL = ("https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/"
           "kamil/oncoanalyser/IPISRC044_T0_personalis/alignments/rna/IPISRC044_tumor_T0_personalis_rna.md.bam")
K562_URL = "https://jbrowse.org/demos/cancer_sv/K562_isoseq.bam"


def sid_window(read):
    """The intrachromosomal RNA join is represented by a CIGAR N, not SA."""
    if read.flag & (4 | 256 | 512 | 1024 | 2048) or read.mapping_quality < 20:
        return None
    mapping = {q: (read.reference_name, g) for q, g in read.get_aligned_pairs(matches_only=True)}
    inverse = {p: q for q, p in mapping.items()}
    if ("chr11", 118401716) not in inverse or ("chr11", 118468774) not in inverse:
        return None
    a, b = inverse["chr11", 118401716] + 1, inverse["chr11", 118468774]
    lo, hi = a - 60, b + 45
    if (a != b or lo < 0 or hi > len(read.query_sequence) or read.query_qualities is None
            or min(read.query_qualities[lo:hi]) < 20 or any(q not in mapping for q in range(lo, hi))):
        return None
    return dict(sequence=read.query_sequence[lo:hi], junction_start=60, junction_end=60,
                blocks=blocks(mapping, lo, hi), source_query_start=lo,
                minimum_base_quality=min(read.query_qualities[lo:hi]), read_id=read.query_name,
                fragment_id=read.query_name, original_sam=read.to_string())


def build(sid_bam, k562_bam, output):
    output.mkdir(parents=True, exist_ok=True)
    genome = EnsemblRelease(87)
    manifest = []
    for event, sample, source, path, genes, donor, acceptor in [
        ("ATP5MG--KMT2A", "Sid-T0-Personalis", SID_URL, sid_bam, ["ATP5L", "KMT2A"],
         ("chr11", 118401717), ("chr11", 118468774)),
        ("BCR--ABL1", "K562-external-control", K562_URL, k562_bam, ["BCR", "ABL1"],
         ("chr22", 23290413), ("chr9", 130854063)),
    ]:
        groups = defaultdict(list)
        with pysam.AlignmentFile(path) as bam:
            for read in bam:
                row = (sid_window(read) if event.startswith("ATP5MG") else
                       extract(read, donor[0], donor[1], acceptor[0], acceptor[1] + 1, 120))
                if row is None:
                    continue
                key = json.dumps({k: row[k] for k in ("sequence", "junction_start", "junction_end", "blocks")}, sort_keys=True)
                if row["read_id"] not in {r["read_id"] for r in groups[key]}:
                    groups[key].append(row)
        rows = sorted(groups.values(), key=lambda r: (-len({v["fragment_id"] for v in r}), r[0]["sequence"]))[0]
        assert len({r["fragment_id"] for r in rows}) >= 2
        row = rows[0]
        models = references(genome, genes)  # All available models, including unavailable CDSs.
        data = dict(fusion=dict(event_id=event, reference_name="GRCh38",
            **{k: row[k] for k in ("sequence", "junction_start", "junction_end", "blocks")},
            donor=dict(contig=donor[0], position=donor[1], strand="+"),
            acceptor=dict(contig=acceptor[0], position=acceptor[1], strand="+"),
            provenance=dict(sample_id=sample, method="identical original aligned RNA windows", version="1",
                parameters=dict(left_flank=row["junction_start"], right_flank=len(row["sequence"]) - row["junction_end"],
                                min_mapq=20, min_base_quality=20 if event.startswith("ATP5MG") else 10),
                source=source, contig_id=event + "-observed-window", source_bam_sha256=sha256(path.read_bytes()).hexdigest(),
                caller=("Isofox Id_14364; Personalis" if event.startswith("ATP5MG") else "BCR exon 14 / ABL1 exon 2"))),
            references=models, reference_names={r["transcript_id"]: genome.transcript_by_id(r["transcript_id"]).name for r in models},
            reads=[dict(sample_id=sample, library_id=sample, fragment_id=r["fragment_id"], read_id=r["read_id"],
                        source=source, source_query_start=r["source_query_start"], cdna_start=0,
                        sequence=r["sequence"], blocks=r["blocks"]) for r in rows],
            original_records=[dict(sam=r["original_sam"], minimum_base_quality=r["minimum_base_quality"]) for r in rows],
            note="Observed local RNA, not a full-length transcript or confirmed tumor-specific protein. All compatible annotation alternatives retained.")
        path = output / (event + ".input.json.gz")
        path.write_bytes(gzip.compress((json.dumps(data, sort_keys=True) + "\n").encode(), mtime=0))
        manifest.append(dict(event=event, sample=sample, input=path.name, sha256=sha256(path.read_bytes()).hexdigest(),
                             fragments=len(rows), expected_status="ambiguous"))
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sid-bam", type=Path, required=True)
    parser.add_argument("--k562-bam", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build(args.sid_bam, args.k562_bam, args.output)
