"""Coordinate-table working sets at every nuclear locus, measured on Isovar 1.8.3 (#234).

Run in a fresh process with PYTHONPATH pointing at a clean v1.8.3 checkout. The
scan instruments that release's per-call LRU coordinate cache to record lookups
and hits, and records the exact size of the dict table Isovar 1.8.4 builds for
the same coordinates on this interpreter. It is a diagnostic for choosing loci
and bounding the table, not a timing benchmark.
"""

import argparse
from functools import lru_cache
import json
import logging
from pathlib import Path
import platform
import sys

import pysam
from varcode import Variant

import isovar
import isovar.read_collector as read_collector
from isovar.variant_helpers import base0_interval_for_variant_fields, trim_variant
from tests.data.osteosarc.expansion.inventory import digest
from tests.data.osteosarc.expansion.references import gzip_bytes


def table_bytes(distinct):
    """Size of a dict filled the way Isovar 1.8.4 fills it, on this interpreter."""
    table = {}
    share = table.setdefault
    for position in range(10**7, 10**7 + distinct):
        share(position, position)
    return sys.getsizeof(table)


def scan(destination):
    if not hasattr(read_collector, "lru_cache"):
        raise ValueError("Run against Isovar 1.8.3, whose collector shares coordinates through an LRU cache")
    caches = []

    def observed(**options):
        def decorate(function):
            wrapped = lru_cache(**options)(function)
            caches.append(wrapped)
            return wrapped
        return decorate

    read_collector.lru_cache = observed
    variants = json.loads((destination / "inventory-GRCh38-validated.json").read_text())["variants"]
    sources = sorted(p for p in (destination / "alignments").iterdir() if (p / "regions-GRCh38.bam").exists())
    sizes, rows, maxsizes = {}, [], set()
    for record in variants:
        if record["chrom"].removeprefix("chr") in ("M", "MT"):
            continue  # The whole mitochondrial genome is below the old bound.
        position, ref, alt = trim_variant(Variant(record["chrom"], record["pos"], record["ref"], record["alt"]))
        start, end = base0_interval_for_variant_fields(position, ref, alt)
        for source in sources:
            with pysam.AlignmentFile(source / "regions-GRCh38.bam") as bam:
                chromosome = read_collector.ReadCollector._infer_chromosome_name(record["chrom"], bam.references)
                if chromosome is None:
                    continue
                caches.clear()
                reads = read_collector.ReadCollector(merge_overlapping_fragments=False).get_locus_reads(
                    bam, chromosome, start, end, position, ref, alt)
            if not reads:
                continue
            info = caches[-1].cache_info()
            maxsizes.add(info.maxsize)
            distinct = len({p for read in reads for p in read.reference_positions if p is not None})
            if distinct not in sizes:
                sizes[distinct] = table_bytes(distinct)
            rows.append(dict(source_id=source.name, variant_id=record["variant_id"], reads=len(reads),
                             lookups=info.hits + info.misses, hits=info.hits, distinct=distinct,
                             table_bytes=sizes[distinct]))
    return rows, sorted(maxsizes)


def run(args):
    destination = Path(args.destination)
    inventory = destination / "inventory-GRCh38-validated.json"
    receipts = {p.parent.name: json.loads(p.read_text())["bam_sha256"]
                for p in sorted((destination / "alignments").glob("*/regions-GRCh38.json"))
                if (p.parent / "regions-GRCh38.bam").exists()}
    disabled = logging.root.manager.disable
    logging.disable(logging.CRITICAL)
    try:
        rows, maxsizes = scan(destination)
    finally:
        logging.disable(disabled)
    report = dict(
        identity=dict(isovar=isovar.__version__, python=platform.python_version(), platform=platform.platform(),
                      scanner_sha256=digest(__file__), inventory_sha256=digest(inventory),
                      receipt_bam_sha256=receipts, lru_maxsize=maxsizes, int_bytes=sys.getsizeof(10**7)),
        rows=rows)
    gzip_bytes(args.output, json.dumps(report, sort_keys=True, indent=1).encode() + b"\n")
    print(json.dumps(dict(rows=len(rows), max_distinct=max(r["distinct"] for r in rows)), sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", required=True)
    parser.add_argument("--output", required=True, help="a new .json.gz file")
    run(parser.parse_args())
