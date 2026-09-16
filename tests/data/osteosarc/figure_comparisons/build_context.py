"""Pin audited original RNA/DNA records for length, phase and haplotype figures."""
import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path

import pysam


def records(path, chrom, lo, hi):
    with pysam.AlignmentFile(path) as bam:
        name = chrom if chrom in bam.references else chrom.removeprefix("chr")
        return [r.to_string() for r in bam.fetch(name,lo,hi)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostics",type=Path,required=True)
    parser.add_argument("--reference-models",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    read = lambda name: json.loads((args.diagnostics/name).read_text())
    znf = [r for r in read("reconstruction.json")["results"] if r["variant_id"] == "ZNF436-chr1-23362762"]
    phase = read("cd109-phase.json")
    map2 = read("map2.json")
    for r in znf:
        path = Path(r.pop("input_path"))
        assert sha256(path.read_bytes()).hexdigest() == r["input_sha256"]
        r["original_records"] = records(path,"chr1",23362761,23362762)
    for r in phase:
        path = Path(r.pop("path"))
        r["source_sha256"] = sha256(path.read_bytes()).hexdigest()
        r["original_records"] = records(path,"chr6",73811000-1,73820000)
    selected = {"d9200d085cb7774d","707928e142bdb21a","261f8be4ce5aa8c8","95a1e72e7ec10990","2bc1fc291308debb"}
    for r in map2["sources"]:
        path = Path(r.pop("path"))
        receipt = json.loads(path.with_suffix(".json").read_text())
        r["source_sha256"] = sha256(path.read_bytes()).hexdigest()
        assert r["source_sha256"] == receipt["bam_sha256"]
        if r["source_id"] in selected:
            r["original_records"] = records(path,"chr2",209694763,209694806)
    models = json.loads(gzip.decompress(args.reference_models.read_bytes()))
    payload = dict(znf436=znf,cd109=phase,map2=map2,verification=read("verification.json"),
                   models={k:models[k] for k in ("ENST00000287097","ENST00000360351","ENST00000447185")},
                   reference_models_sha256=sha256(args.reference_models.read_bytes()).hexdigest())
    args.output.mkdir(parents=True,exist_ok=True)
    output = args.output/"context-evidence.json.gz"
    output.write_bytes(gzip.compress(json.dumps(payload,sort_keys=True).encode(),mtime=0))
    (args.output/"context-manifest.json").write_text(json.dumps(dict(
        file=output.name,sha256=sha256(output.read_bytes()).hexdigest(),
        description="Original source records plus audited exact local sequence/quality and independent translation checks"),indent=2)+"\n")
    print(output, output.stat().st_size)


if __name__ == "__main__":
    main()
