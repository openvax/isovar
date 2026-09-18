"""DLG5 DNA calls and bounded original tumor/normal alignment evidence."""

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import gzip
import json
from pathlib import Path
import re

import pysam

from tests.data.osteosarc.expansion.inventory import digest, fetch_snapshot
from tests.data.osteosarc.expansion.liftover import genomic_sequence
from tests.data.osteosarc.figure_comparisons.footprints import ROOT, acquire, deletion_footprint, split_deletion_paths
from tests.data.fusions.build_osteosarc import extract, segment_key

PREFIXES = {
    "T0": "personalis/normal-personalis_tumor-personalis-T0_ds.766d6a422f244f43a9c055851dd75eab/normal-personalis_tumor-personalis-T0",
    "T1": "ucla_T1/normal-ucla_tumor-ucla-T1_ds.6bf94c171a8f4af68274fcbb2b9b3e01/normal-ucla_tumor-ucla-T1",
    "T2": "ucla_T2/normal-ucla_tumor-ucla-T2_ds.0babb746c09042c79869258feb9391fc/normal-ucla_tumor-ucla-T2",
    "T1-organoid": "ucla_T1_organoid/normal-ucla_tumor-ucla-T1-organoid_ds.b3ae0836c914464aa9390aea60f1778a/normal-ucla_tumor-ucla-T1-organoid",
}
INTERVAL = (77850921, 77930452)
REGIONS = ["chr10:77845922-77935452"]
START_CODON = (77926517, 77926520)
EXCLUDE = 4 | 256 | 512 | 1024 | 2048


def local_evidence(records, signatures, start_codon=START_CODON):
    """Exact sequence evidence is separate from splice/deletion CIGAR evidence.

    A primary sequence can retain the junction even when its supplementary
    record was dropped. Count that sequence explicitly without fabricating
    a chimeric path. QNAME/RG templates are not independent molecules.
    """
    supported = {name: defaultdict(set) for name in signatures}
    start_states, original, quality_unresolved = defaultdict(set), {}, set()
    for r in records:
        if r.flag & EXCLUDE or r.mapping_quality < 20 or r.query_sequence is None:
            continue
        key = segment_key(r)[:2]
        sequence, quality = r.query_sequence, r.query_qualities
        for name, signature in signatures.items():
            positions = [i for i in range(len(sequence) - len(signature) + 1) if sequence.startswith(signature, i)]
            for i in positions:
                supported[name]["sequence_only"].add(key)
                original[r.to_string()] = None
                for q in (10, 20, 30):
                    if quality is not None and min(quality[i:i + len(signature)]) >= q:
                        supported[name]["q%d" % q].add(key)
        mapping = {g: q for q, g in r.get_aligned_pairs(matches_only=True)}
        positions = [mapping.get(g) for g in range(*start_codon)]
        if any(q is None for q in positions) or positions != list(range(positions[0], positions[0] + 3)):
            continue
        original[r.to_string()] = None
        if quality is None or min(quality[q] for q in positions) < 20:
            quality_unresolved.add(key)
        else:
            # DLG5 is minus-strand; genomic CAT corresponds to coding ATG.
            state = "reference_start" if sequence[positions[0]:positions[-1] + 1] == "CAT" else "other"
            start_states[key].add(state)
    start_counts = dict(Counter(next(iter(v)) if len(v) == 1 else "conflicting_template"
                               for v in start_states.values()))
    if quality_unresolved.difference(start_states):
        start_counts["quality_unresolved"] = len(quality_unresolved.difference(start_states))
    return dict(signatures={name: {label: len(states[label]) for label in ("sequence_only", "q10", "q20", "q30")}
                            for name, states in supported.items()},
                start_codon=start_counts,
                original_records=sorted(original))


def validate_context(context, calls):
    """Preserve caller disagreement; never repair VCF fields to fit the reads."""
    references = [json.loads(p.read_text()) for p in sorted(context.glob("hg38-*.json")) if ".receipt." not in p.name]
    def sequence(start, end):
        ref, = [r for r in references if r["start"] <= start and r["end"] >= end]
        return ref["dna"][start - ref["start"]:end - ref["start"]].upper()

    dragen = []
    for call in calls:
        for raw in call["records"]:
            fields = raw.split("\t")
            info = dict(field.split("=", 1) for field in fields[7].split(";") if "=" in field)
            start, end = int(fields[1]), int(info["END"])
            if sequence(start - 1, start) != fields[3]:
                raise ValueError("DLG5 DRAGEN reference anchor disagrees")
            signature = sequence(start - 20, start) + info["SVINSSEQ"] + sequence(end, end + 20)
            dragen.append(dict(sample=call["sample"], interval=[start, end], insertion=info["SVINSSEQ"],
                               signature=signature, contig=info["CONTIG"], signature_in_contig=signature in info["CONTIG"]))
    purple = []
    for path in sorted(context.glob("*.purple.sv.vcf.gz")):
        with gzip.open(path, "rt") as handle:
            lines = list(handle)
        selected = []
        for raw in lines:
            if raw.startswith("#"):
                continue
            fields = raw.strip().split("\t")
            if fields[0] != "chr10" or not abs(int(fields[1]) - INTERVAL[0]) <= 50:
                continue
            match = re.fullmatch(r"([ACGT]+)\[chr10:(\d+)\[", fields[4])
            if not match or abs(int(match[2]) - INTERVAL[1]) > 50:
                continue
            start, end = int(fields[1]), int(match[2]) - 1
            if not match[1].startswith(fields[3]) or sequence(start - 1, start) != fields[3]:
                raise ValueError("DLG5 BND reference anchor disagrees")
            insertion = match[1][len(fields[3]):]
            signature = sequence(start - 20, start) + insertion + sequence(end, end + 20)
            selected.append(dict(record=raw.strip(), interval=[start, end], insertion=insertion, signature=signature))
        purple.append(dict(sample=path.name.split(".")[0], records=selected,
                           receipt=json.loads(path.with_name(path.name + ".receipt.json").read_text())))
    signature, = {r["signature"] for p in purple for r in p["records"]}
    return dict(references=references, dragen=dragen, purple=purple,
                signatures={"dragen_fields": dragen[0]["signature"], "purple_bnd": signature},
                purple_signature_in_all_dragen_contigs=all(signature in d["contig"] for d in dragen))


def audit(dna, rna, context, footprint_audit, output):
    output.mkdir(parents=True, exist_ok=False)
    calls = json.loads((dna / "calls.json").read_text())
    evidence = validate_context(context, calls)
    evidence["calls"] = calls
    base = json.loads((footprint_audit / "evidence.json").read_text())
    models = next(e["models"] for e in base["deletions"] if e["gene"] == "DLG5")
    evidence["models"], evidence["products"] = models, []
    intervals = {"dragen_fields": INTERVAL, "purple_bnd": (77850914, 77930460), "observed_cigar": (77850921, 77930460)}
    for directory in (dna, rna):
        sources = json.loads((directory / "sources.json").read_text())
        for source in sources:
            receipt = json.loads((directory / source["id"] / "receipt.json").read_text())
            path = directory / source["id"] / "regions.bam"
            if digest(path) != receipt["bam_sha256"]:
                raise ValueError("DLG5 source drift: " + source["id"])
            with pysam.AlignmentFile(path) as bam:
                records = list(bam.fetch("chr10", 77845921, 77935452))
                header = str(bam.header)
            result = dict(source=source, receipt=receipt, header=header,
                          **local_evidence(records, evidence["signatures"]))
            result["footprints"] = {label: deletion_footprint(records, dict(start=a, end=b), models)
                                    for label, (a, b) in intervals.items()}
            result["split_paths"] = {label: split_deletion_paths(records, a, b) for label, (a, b) in intervals.items()}
            groups = defaultdict(list)
            for r in records:
                if r.has_tag("SA"):
                    groups[segment_key(r)].append(r)
            windows = []
            for group in groups.values():
                for r in group:
                    # The exact BND's retained bases, not rounded gene labels.
                    window = extract(r, "chr10", 77850914, "chr10", 77930461, 20, group)
                    if window:
                        windows.append(window)
            result["observed_windows"] = windows
            evidence["products"].append(result)
            print(source["id"], result["signatures"], result["start_codon"], len(windows), flush=True)
    (output / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    packed = output / "dlg5.json.gz"
    packed.write_bytes(gzip.compress(json.dumps(evidence, sort_keys=True).encode(), mtime=0))
    (output / "dlg5-manifest.json").write_text(json.dumps(dict(file=packed.name, sha256=digest(packed)), indent=2) + "\n")
    return evidence


def acquire_context(output):
    output.mkdir(parents=True, exist_ok=False)
    for start, end in [(77850721, 77851121), (77930252, 77930652), (77926497, 77926540)]:
        genomic_sequence(output, "hg38", "chr10", start, end)
    for sample in ("T0_personalis", "T1_ucla", "T2_ucla", "T1_organoid_ucla"):
        prefix = "IPISRC044_tumor_" + sample
        url = ROOT + "kamil/oncoanalyser/IPISRC044_" + sample + "/purple/" + prefix + ".purple.sv.vcf.gz"
        fetch_snapshot(url, output / (sample + ".purple.sv.vcf.gz"))


def acquire_dna(output):
    output.mkdir(parents=True, exist_ok=False)
    calls, sources = [], []
    for sample, prefix in PREFIXES.items():
        url = ROOT + "kamil/basespace/results/" + prefix
        vcf = output / (sample + ".sv.vcf.gz")
        receipt = fetch_snapshot(url + ".sv.vcf.gz", vcf)
        with pysam.VariantFile(vcf) as handle:
            selected = [str(r).strip() for r in handle
                        if r.contig == "chr10" and
                        abs(r.pos - INTERVAL[0]) <= 50 and abs(r.stop - INTERVAL[1]) <= 50]
            calls.append(dict(sample=sample, receipt=receipt, header=str(handle.header), records=selected))
        sources.append(dict(id=sample + "-DNA", sample=sample, product="DNA", tissue="tumor",
                            url=url + "_tumor.bam", processing_family=sample + "-DNA"))
        if sample in ("T0", "T1"):
            sources.append(dict(id=sample + "-normal-DNA", sample=sample, product="DNA", tissue="normal",
                                url=url + ".bam", processing_family=sample + "-normal-DNA"))
    (output / "calls.json").write_text(json.dumps(calls, indent=2) + "\n")

    def fetch(source):
        index = output / (source["id"] + ".source.bai")
        source["index_receipt"] = fetch_snapshot(source["url"] + ".bai", index)
        return acquire(source, output, query_regions=REGIONS, index_path=index.resolve())

    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(fetch, sources))
    (output / "sources.json").write_text(json.dumps(sources, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--acquire", type=Path)
    group.add_argument("--context", type=Path)
    group.add_argument("--audit", type=Path)
    parser.add_argument("--dna", type=Path)
    parser.add_argument("--rna", type=Path)
    parser.add_argument("--context-inputs", type=Path)
    parser.add_argument("--footprints", type=Path)
    args = parser.parse_args()
    if args.acquire:
        acquire_dna(args.acquire)
    elif args.context:
        acquire_context(args.context)
    else:
        if any(x is None for x in (args.dna, args.rna, args.context_inputs, args.footprints)):
            parser.error("--audit requires --dna, --rna, --context-inputs and --footprints")
        audit(args.dna, args.rna, args.context_inputs, args.footprints, args.audit)
