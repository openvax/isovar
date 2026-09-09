"""Validate a coordinate-shifted ND5 annotation for CeGaT's hg19 MT.

Only ungapped forward mappings are accepted. Original Ensembl cDNA/protein
expectations remain unchanged: historical hg19 reference substitutions are
recorded explicitly, not silently injected into the biological expectation.
This is not a general annotation liftover or sequence substitution.
"""

import argparse
from copy import deepcopy
import gzip
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from tests.data.osteosarc.expansion.inventory import digest, write_json  # noqa: E402
from tests.data.osteosarc.expansion.liftover import chain_blocks, genomic_sequence, map_interval  # noqa: E402
from tests.data.osteosarc.expansion.references import (  # noqa: E402
    apply_variant, gzip_bytes, load_reference, raw_features, reverse_complement,
)


def shifted_interval(blocks, start, end):
    result = map_interval(blocks, "chrM", start - 1, end)
    if result["chrom"] != "chrM" or result["strand"] != "+":
        raise ValueError("Only same-contig forward mitochondrial annotation shifts are supported")
    return result["start"] + 1, result["end"]


def build(destination, reference, output):
    destination, reference, output = map(Path, (destination, reference, output))
    original_manifest, original_models = load_reference(reference)
    if original_manifest["assembly"] != "GRCh37":
        raise ValueError("CeGaT requires the validated GRCh37 nuclear reference")
    manifest, models = deepcopy(original_manifest), deepcopy(original_models)
    chain = destination / "genomic-references/hg38ToHg19.over.chain.gz"
    with gzip.open(chain, "rt") as handle:
        blocks = list(chain_blocks(handle, {"chrM"}))
    validations = {}
    for tid, model in models.items():
        if model["contig"] != "MT":
            continue
        model["exons"] = [shifted_interval(blocks, start, end) for start, end in model["exons"]]
        sequences, receipts = [], []
        for start, end in sorted(model["exons"]):
            sequence, receipt = genomic_sequence(destination / "genomic-references", "hg19", "chrM", start - 1, end)
            sequences.append(sequence)
            receipts.append(receipt)
        sequence = "".join(sequences)
        if model["strand"] == "-":
            sequence = reverse_complement(sequence)
        if len(sequence) != len(model["cdna"]) or set(sequence) - set("ACGT"):
            raise ValueError("Mitochondrial reference length/sequence cannot validate this coding frame")
        differences = [dict(cdna_offset=i, ensembl_base=a, hg19_base=b)
                       for i, (a, b) in enumerate(zip(model["cdna"], sequence)) if a != b]
        validations[tid] = dict(receipts=receipts, reference_substitutions=differences,
                                expectation="original Ensembl cDNA/protein, not historical hg19 haplotype")
    if not validations:
        raise ValueError("No independently validated mitochondrial transcript")
    lines = []
    for line, fields, _ in raw_features(reference / "reference.gtf.gz"):
        if fields is not None and fields[0] == "MT":
            start, end = shifted_interval(blocks, int(fields[3]), int(fields[4]))
            fields[3:5] = [str(start), str(end)]
            line = "\t".join(fields) + "\n"
        lines.append(line)
    inventory = json.loads((destination / "inventory-GRCh37-hg19MT-validated.json").read_text())
    for record in inventory["variants"]:
        for tid in manifest["variant_transcripts"][record["variant_id"]]:
            apply_variant(record, models[tid])
    output.mkdir(parents=True, exist_ok=False)
    gzip_bytes(output / "reference.gtf.gz", "".join(lines).encode())
    gzip_bytes(output / "models.json.gz", json.dumps(models, sort_keys=True).encode())
    for name in ("reference.cdna.fa.gz", "reference.pep.fa.gz"):
        shutil.copyfile(reference / name, output / name)
    manifest["files"] = {name: digest(output / name) for name in manifest["files"]}
    manifest["dataset_identity"] += "-hg19MT-validated-coordinate-shift"
    manifest["mitochondrial_annotation_derivation"] = dict(
        original_manifest_sha256=digest(reference / "manifest.json"), chain_sha256=digest(chain),
        sequence_validations=validations,
        method="forward ungapped coordinate mapping; full-transcript comparison with explicit reference substitutions; no cDNA edits",
    )
    write_json(output / "manifest.json", manifest)
    print("Validated CeGaT mitochondrial coordinate-shifted annotation", sorted(validations), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    build(args.destination, args.reference, args.output)
