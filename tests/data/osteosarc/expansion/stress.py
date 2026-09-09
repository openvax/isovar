"""Separately labelled original-RNA insertion/compound-allele stress probes."""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import json
import logging
from pathlib import Path
import shutil
import sys

import pysam
from varcode import Variant

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from tests.data.osteosarc.expansion.acquire import acquire_regions  # noqa: E402
from tests.data.osteosarc.expansion.inventory import digest, write_json  # noqa: E402
from tests.data.osteosarc.expansion.liftover import verify_allele  # noqa: E402
from tests.data.osteosarc.expansion.references import apply_variant, load_reference, minimal_edit, reference_genome  # noqa: E402
from tests.data.osteosarc.expansion.runner import audit_mode, PRIMARY_EXCLUDE_FLAGS  # noqa: E402
from tests.real_rna_helpers import record_digest  # noqa: E402


SOURCES = ("f30f618fb76a0e49", "7010dd389c50467d", "0066232879babe83", "1b66c15da594a3ef")


def exact_allele(read, record):
    """Walk CIGAR with external anchors; N is never a genomic deletion."""
    pos, ref, _ = minimal_edit(record)
    if read.is_unmapped or read.query_sequence is None or read.reference_name != record["chrom"]:
        return None
    left, right = pos - 2, pos - 1 + len(ref)
    positions, skips = {}, []
    rpos, qpos = read.reference_start, 0
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            for target in (left, right):
                if rpos <= target < rpos + length:
                    positions[target] = qpos + target - rpos
        elif op == 3:
            skips.append((rpos, rpos + length))
        if op in (0, 1, 4, 7, 8):
            qpos += length
        if op in (0, 2, 3, 7, 8):
            rpos += length
    if left not in positions or right not in positions or any(a <= right and b > left for a, b in skips):
        return None
    return read.query_sequence[positions[left] + 1:positions[right]]


def records(destination):
    destination = Path(destination)
    inventory = json.loads((destination / "inventory-GRCh38-validated.json").read_text())
    vaccine_ids = {v["variant_id"] for v in inventory["variants"]}
    selected = {}
    for row in csv.DictReader((destination / "vafs.tsv").open(), delimiter="\t"):
        if row["variant_id"] not in ("ACSL6-chr5-131988563", "KTN1-chr14-55627965", "EPPK1-chr8-143870409"):
            continue
        variant = {k: row[k] for k in ("variant_id", "gene", "chrom", "ref", "alt")}
        variant.update(pos=int(row["pos"]), assembly="GRCh38", cohort="non_vaccine_stress", vaccine_count=0)
        if variant["variant_id"] in vaccine_ids:
            raise ValueError("Stress cohort unexpectedly overlaps vaccine inventory")
        if row["variant_id"] in selected and selected[row["variant_id"]] != variant:
            raise ValueError("Ambiguous stress allele")
        selected[row["variant_id"]] = variant
    ntf3 = next(v for v in inventory["variants"] if v["variant_id"] == "NTF3-chr12-5494381")
    metadata = json.loads((destination / "source-metadata.json").read_text())
    details = next(v for v in metadata["variants"] if v["variant_id"] == ntf3["variant_id"])
    definition = details["source_sections"]["Genomic context"]
    if "cDNA change AG>GT Refseq strand fwd" not in definition:
        raise ValueError("NTF3 source compound definition changed")
    selected["NTF3-compound"] = dict(ntf3, variant_id=ntf3["variant_id"] + "-compound",
                                     ref="AG", alt="GT", cohort="source_linked_compound_stress", vaccine_count=None,
                                     membership_scope="original vaccine row is the SNV; compound allele is separately tested",
                                     source_definition_receipt=details["receipt"])
    validated = [verify_allele(v, "hg38", destination / "genomic-references") for v in selected.values()]
    return dict(inventory, variants=validated)


def prepare(destination):
    destination = Path(destination)
    inventory = records(destination)
    output = destination / "inventory-stress.json"
    if output.exists():
        if json.loads(output.read_text()) != inventory:
            raise ValueError("Stress definition drift")
    else:
        write_json(output, inventory)
    sources = [s for s in inventory["alignments"] if s["source_id"] in SOURCES]
    with ThreadPoolExecutor(max_workers=4) as executor:
        for result in executor.map(lambda s: acquire_regions(s, inventory["variants"], destination, "GRCh38", "GRCh38-stress"), sources):
            print("STRESS SOURCE", result["source_id"], result["status"], flush=True)


def observations(destination):
    destination = Path(destination)
    inventory = json.loads((destination / "inventory-stress.json").read_text())
    rows = []
    for sid in SOURCES:
        directory = destination / "alignments" / sid
        receipt = json.loads((directory / "regions-GRCh38-stress.json").read_text())
        if receipt["status"] != "ok":
            raise ValueError("Stress acquisition not complete")
        path = directory / "regions-GRCh38-stress.bam"
        if digest(path) != receipt["bam_sha256"]:
            raise ValueError("Stress source checksum mismatch")
        with pysam.AlignmentFile(path) as bam:
            for record in inventory["variants"]:
                counts, witnesses, total = Counter(), [], 0
                _, ref, alt = minimal_edit(record)
                for read in bam.fetch(record["chrom"], record["pos"] - 2, record["pos"] + len(record["ref"])):
                    total += 1
                    if read.flag & PRIMARY_EXCLUDE_FLAGS or read.mapping_quality < 1:
                        continue
                    observed = exact_allele(read, record)
                    group = "uncallable" if observed is None else "ref" if observed == ref else "alt" if observed == alt else "other"
                    counts[group] += 1
                    if group == "alt":
                        witnesses.append(dict(sam_sha256=record_digest(read), name=read.query_name))
                rows.append(dict(source_id=sid, variant_id=record["variant_id"], region_alignments=total,
                                 counts={g: counts[g] for g in ("ref", "alt", "other", "uncallable")}, witnesses=witnesses,
                                 source_bam_sha256=receipt["bam_sha256"]))
                print("STRESS", sid, record["variant_id"], counts, flush=True)
    write_json(destination / "stress-observations.json", rows)


def build_corpus(destination, output, cache):
    """Copy complete original stress regions; no base/tag/quality rewriting."""
    destination, output, cache = map(Path, (destination, output, cache))
    logging.disable(logging.CRITICAL)
    inventory = json.loads((destination / "inventory-stress.json").read_text())
    observations = json.loads((destination / "stress-observations.json").read_text())
    sources = {s["source_id"]: s for s in inventory["alignments"]}
    output.mkdir(parents=True, exist_ok=False)
    shutil.copytree(destination / "reference-stress-GRCh38", output / "reference")
    manifest, models = load_reference(output / "reference")
    genome = reference_genome(output / "reference", cache)
    cases = []
    for sid in SOURCES:
        directory = destination / "alignments" / sid
        source = directory / "regions-GRCh38-stress.bam"
        receipt = json.loads(source.with_suffix(".json").read_text())
        if receipt["status"] != "ok" or digest(source) != receipt["bam_sha256"]:
            raise ValueError("Stress original source failed integrity validation")
        for record in inventory["variants"]:
            name = f"{record['variant_id']}-{sid}.bam"
            path = output / name
            selected = []
            with pysam.AlignmentFile(source) as bam:
                with pysam.AlignmentFile(path, "wb", header=bam.header) as target:
                    for ordinal, read in enumerate(bam.fetch(record["chrom"], record["pos"] - 2,
                                                              record["pos"] + len(record["ref"]))):
                        target.write(read)
                        selected.append(dict(regional_ordinal=ordinal, sam_sha256=record_digest(read)))
            pysam.index(str(path))
            variant = Variant(record["chrom"].removeprefix("chr"), record["pos"], record["ref"], record["alt"], ensembl=genome)
            expected = {tid: apply_variant(record, models[tid]) for tid in manifest["variant_transcripts"][record["variant_id"]]}
            result = audit_mode(path, variant, expected)
            if result["status"] != "ok" or result["outcome"] == "protein_validation_error":
                raise ValueError(f"Stress pipeline/validation failure: {record['variant_id']}: {result}")
            observation = next(r for r in observations if r["source_id"] == sid and r["variant_id"] == record["variant_id"])
            cases.append(dict(case_id=name, bam=name, variant=record, source_id=sid, source_url=sources[sid]["url"],
                              source_bam_sha256=receipt["bam_sha256"], selected_records=selected,
                              files={p.name: digest(p) for p in (path, Path(str(path) + ".bai"))},
                              independent_primary=observation, defaults=result))
            print("STRESS FIXTURE", name, result["outcome"], flush=True)
    write_json(output / "manifest.json", dict(cases=cases, selection="complete original expanded region, not vaccine membership",
                                             inventory_sha256=digest(destination / "inventory-stress.json")))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--observe", action="store_true")
    parser.add_argument("--corpus", type=Path)
    parser.add_argument("--cache", type=Path)
    args = parser.parse_args()
    if args.corpus:
        if args.cache is None or args.observe:
            parser.error("Corpus creation requires --cache and cannot combine --observe")
        build_corpus(args.destination, args.corpus, args.cache)
    else:
        (observations if args.observe else prepare)(args.destination)
