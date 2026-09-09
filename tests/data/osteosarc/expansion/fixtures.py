"""Build small, deliberately selected original-read regression fixtures.

Selection is allele-balanced and may favor names supporting a returned
protein. It is NOT a random sample or an estimator of full-source counts.
Records are copied unchanged, including flags, tags and multiplicity.
"""

import argparse
from collections import defaultdict
import json
import logging
from pathlib import Path
import shutil
import sys

import pysam
from varcode import Variant

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from tests.data.osteosarc.expansion.inventory import digest, write_json  # noqa: E402
from tests.data.osteosarc.expansion.references import (  # noqa: E402
    apply_variant, load_reference, minimal_edit, reference_genome,
)
from tests.data.osteosarc.expansion.runner import audit_mode, independent_counts, PRIMARY_EXCLUDE_FLAGS  # noqa: E402
from tests.real_rna_helpers import cigar_observation, record_digest  # noqa: E402


def select_names(bam, record, preferred_names=(), limit=16):
    contig = record["chrom"].removeprefix("chr")
    names = set(bam.references) & ({"MT", "M", "chrM", "chrMT"} if contig in ("M", "MT") else {contig, "chr" + contig})
    if len(names) != 1:
        raise ValueError("Fixture requires an unambiguous native contig")
    native = dict(record, chrom=next(iter(names)))
    region = (native["chrom"], max(0, record["pos"] - 2), record["pos"] + len(record["ref"]))
    _, ref, alt = minimal_edit(record)
    groups = defaultdict(dict)
    for read in bam.fetch(*region):
        observation = cigar_observation(read, native)
        group = ("uncallable" if observation is None else "ref" if observation["allele"] == ref else
                 "alt" if observation["allele"] == alt else "other")
        groups[group][read.query_name] = min(groups[group].get(read.query_name, "f" * 64), record_digest(read))
    selected = set(preferred_names[:32])
    for values in groups.values():
        selected.update(sorted(values, key=lambda name: (values[name], name))[:limit])
    return selected, region, native


def original_fixture(source_bam, output, record, preferred_names=(), complete=False):
    with pysam.AlignmentFile(source_bam) as source:
        selected, region, native = select_names(source, record, preferred_names)
        records = []
        with pysam.AlignmentFile(output, "wb", header=source.header) as target:
            for ordinal, read in enumerate(source.fetch(*region)):
                if complete or read.query_name in selected:
                    target.write(read)
                    records.append(dict(regional_ordinal=ordinal, sam_sha256=record_digest(read), name=read.query_name))
    pysam.index(str(output))
    return native, records


def choose_cases(destination, baseline):
    inventory = json.loads((destination / "inventory-GRCh38-validated.json").read_text())
    candidates = defaultdict(list)
    for path in baseline.glob("*/*-chr*.json"):
        row = json.loads(path.read_text())
        result = row.get("defaults", {})
        if result.get("status") == "ok":
            candidates[row["variant_id"]].append(row)
    cases = []
    for record in inventory["variants"]:
        rows = candidates[record["variant_id"]]
        if not rows:
            raise ValueError("No baseline selection record for vaccine variant")
        if record["chrom"] == "chrM":
            chosen = next(r for r in rows if r["source_id"] == "f30f618fb76a0e49")
        else:
            chosen = max(rows, key=lambda r: (bool(r["defaults"]["proteins"]), r["defaults"]["counts"]["reads"]["alt"],
                                              r["source_id"] == "f30f618fb76a0e49", r["source_id"]))
        preferred = chosen["defaults"]["proteins"][0]["supporting_read_names"] if chosen["defaults"]["proteins"] else []
        cases.append(dict(record=record, source_id=chosen["source_id"], reference="GRCh38",
                          preferred_names=preferred, complete_region=record["chrom"] == "chrM",
                          selection_row_sha256=digest(baseline / chosen["source_id"] / (record["variant_id"] + ".json"))))
    # Independent native-assembly and long-read mitochondrial examples. Keep
    # full deep-source counts in the report, not in these bounded fixtures.
    for reference, inventory_name, sid in (
        ("GRCh37", "inventory-GRCh37-validated.json", "238f917202dd9997"),
        ("GRCh37-cegat", "inventory-GRCh37-hg19MT-validated.json", "fbeb5d8a415e4d18"),
        ("GRCh38", "inventory-GRCh38-validated.json", "7010dd389c50467d"),
        ("GRCh38", "inventory-GRCh38-validated.json", "0066232879babe83"),
    ):
        native = json.loads((destination / inventory_name).read_text())
        record = next(r for r in native["variants"] if r["chrom"] == "chrM")
        cases.append(dict(record=record, source_id=sid, reference=reference, preferred_names=[], complete_region=False))
    # This native GRCh37 source has genuine callable NR2F2 alternate evidence
    # although none of the checked GRCh38 products did. Its additional RNA
    # differences must not be relabelled as an absent alternate or bad code.
    native = json.loads((destination / "inventory-GRCh37-hg19MT-validated.json").read_text())
    record = next(r for r in native["variants"] if r["variant_id"] == "NR2F2-chr15-96332299")
    cases.append(dict(record=record, source_id="fbeb5d8a415e4d18", reference="GRCh37-cegat",
                      preferred_names=[], complete_region=False))
    return inventory, cases


def build(destination, baseline, output, cache):
    destination, baseline, output, cache = map(Path, (destination, baseline, output, cache))
    logging.disable(logging.CRITICAL)
    output.mkdir(parents=True, exist_ok=False)
    reference_sources = {"GRCh38": "reference-GRCh38", "GRCh37": "reference-GRCh37-v2",
                         "GRCh37-cegat": "reference-GRCh37-cegat"}
    references = {}
    for name, source in reference_sources.items():
        shutil.copytree(destination / source, output / "references" / name)
        manifest, models = load_reference(output / "references" / name)
        references[name] = (manifest, models, reference_genome(output / "references" / name, cache / name))
    inventory, cases = choose_cases(destination, baseline)
    sources = {s["source_id"]: s for s in inventory["alignments"]}
    results = []
    for index, case in enumerate(cases):
        record, sid = case["record"], case["source_id"]
        manifest, models, genome = references[case["reference"]]
        assembly = manifest["assembly"]
        directory = destination / "alignments" / sid
        receipt = json.loads((directory / f"regions-{assembly}.json").read_text())
        source_bam = directory / f"regions-{assembly}.bam"
        if receipt["status"] != "ok" or digest(source_bam) != receipt["bam_sha256"]:
            raise ValueError("Fixture source integrity failure")
        name = f"{index:02d}-{record['variant_id']}-{sid}"
        bam = output / (name + ".bam")
        native, selected = original_fixture(source_bam, bam, record, case["preferred_names"], case["complete_region"])
        primary = output / (name + ".primary.bam")
        pysam.view("--no-PG", "-b", "-F", str(PRIMARY_EXCLUDE_FLAGS), "-o", str(primary), str(bam), catch_stdout=False)
        pysam.index(str(primary))
        contig = "MT" if record["chrom"] == "chrM" else record["chrom"].removeprefix("chr")
        variant = Variant(contig, record["pos"], record["ref"], record["alt"], ensembl=genome)
        expectations = {tid: apply_variant(record, models[tid]) for tid in manifest["variant_transcripts"][record["variant_id"]]}
        modes = {"defaults": audit_mode(bam, variant, expectations), "primary_only": audit_mode(primary, variant, expectations)}
        for mode in modes.values():
            if mode["status"] != "ok" or mode["outcome"] == "protein_validation_error":
                raise ValueError(f"Fixture exposed a pipeline/validation failure: {name}: {mode}")
        with pysam.AlignmentFile(bam) as handle:
            observations = independent_counts(handle, native)
        result = dict(
            case_id=name, variant=record, source_id=sid, reference=case["reference"],
            source_url=sources[sid]["url"], source_region_sha256=receipt["bam_sha256"],
            source_receipt_sha256=digest(directory / f"regions-{assembly}.json"),
            selection_row_sha256=case.get("selection_row_sha256"),
            bam=bam.name, primary_bam=primary.name, selected_records=selected,
            selection="complete expanded locus" if case["complete_region"] else "allele-balanced hash-selected names plus up to 32 protein-supporting names",
            independent_primary=observations, **modes,
        )
        result["files"] = {p.name: digest(p) for p in (bam, primary, Path(str(bam) + ".bai"), Path(str(primary) + ".bai"))}
        results.append(result)
        print("FIXTURE", name, len(selected), modes["defaults"]["outcome"], flush=True)
    write_json(output / "manifest.json", dict(
        schema_version=1, selection_baseline="pre-fix full-region matrix; selection only, not expected translation truth",
        source_inventory_sha256=digest(destination / "inventory-GRCh38-validated.json"), cases=results,
        scope="44 vaccine variants plus four mitochondrial platform/reference cases and native GRCh37 NR2F2; fixtures are not full-source counts",
    ))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("destination", "baseline", "output", "cache"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    build(args.destination, args.baseline, args.output, args.cache)
