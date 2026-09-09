"""Full-depth #227 benchmark: separate production and independent checking.

Run as ``python -m tests.data.osteosarc.expansion.benchmark --help``.
Each invocation uses a fresh process and a new output directory. Input BAMs
are verified against acquisition receipts; neither reads nor defaults change.
"""

import argparse
from contextlib import contextmanager
import cProfile
import gzip
from hashlib import sha256
import json
import logging
import math
from pathlib import Path
import platform
import pstats
import resource
import time

import pysam
from varcode import Variant

import isovar
from isovar import __version__, run_isovar
from tests.data.osteosarc.expansion.inventory import digest, write_json
from tests.data.osteosarc.expansion.references import apply_variant, load_reference, reference_genome
from tests.data.osteosarc.expansion.runner import (
    AuditTimeout, TracedCollector, TracedCreator, counts_from_evidence,
    primary_alignment, protein_check, time_limit,
)


class BenchmarkCreator(TracedCreator):
    def variant_sequences_from_reads(self, *args, **kwargs):
        self.rna_sequences = super().variant_sequences_from_reads(*args, **kwargs)
        return self.rna_sequences


class BenchmarkCollector(TracedCollector):
    def read_evidence_for_variant(self, *args, **kwargs):
        started = time.perf_counter()
        cpu_started = time.process_time()
        result = super().read_evidence_for_variant(*args, **kwargs)
        self.trace["read_collection_seconds"] = time.perf_counter() - started
        self.trace["read_collection_cpu_seconds"] = time.process_time() - cpu_started
        self.trace["read_collection_peak_rss_bytes"] = peak_rss_bytes()
        return result


def peak_rss_bytes():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss if platform.system() == "Darwin" else rss * 1024


@contextmanager
def measure(name, output, measurements, timeout, profile):
    profiler = cProfile.Profile() if profile else None
    start = time.perf_counter()
    cpu_start = time.process_time()
    status = "ok"
    try:
        if profiler:
            profiler.enable()
        with time_limit(timeout):
            yield
    except BaseException:
        status = "interrupted"
        raise
    finally:
        if profiler:
            profiler.disable()
        measurements[name] = dict(wall_seconds=time.perf_counter() - start,
                                  process_cpu_seconds=time.process_time() - cpu_start,
                                  process_peak_rss_bytes=peak_rss_bytes(), status=status)
        if profiler:
            profiler.dump_stats(str(output / (name + ".prof")))
            with (output / (name + ".txt")).open("w") as handle:
                pstats.Stats(profiler, stream=handle).strip_dirs().sort_stats("cumulative").print_stats(60)
        write_json(output / (name + ".json"), measurements[name])
        print(json.dumps(dict(phase=name, **measurements[name]), sort_keys=True), flush=True)


def rna_fingerprint(sequences):
    """Include ordered candidates, every original read field, and coverage."""
    read_hashes = {}
    checksum = sha256()
    for sequence in sequences:
        reads = []
        for read in sequence.reads:
            # Objects remain alive in sequences throughout this pass. An
            # identity-keyed cache avoids repeatedly hashing full read fields;
            # only content hashes, never identities, enter the fingerprint.
            key = id(read)
            if key not in read_hashes:
                values = [read.prefix, read.allele, read.suffix, read.name, read.source_read_count]
                read_hashes[key] = sha256(json.dumps(values, separators=(",", ":")).encode()).hexdigest()
            reads.append(read_hashes[key])
        row = [sequence.prefix, sequence.alt, sequence.suffix,
               sorted(reads), sequence.coverage().tolist()]
        checksum.update(json.dumps(row, separators=(",", ":")).encode() + b"\n")
    return checksum.hexdigest()


def run(args):
    if not __debug__:
        raise ValueError("Independent validation requires normal Python mode (no -O)")
    if not math.isfinite(args.timeout) or not 0 < args.timeout <= 1800:
        raise ValueError("Use a finite timeout greater than zero and at most 1800 seconds")
    logging.disable(logging.CRITICAL)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    destination = Path(args.destination)
    source = destination / "alignments" / args.source_id
    bam = source / "regions-GRCh38.bam"
    receipt = source / "regions-GRCh38.json"
    checksum = digest(bam)
    if checksum != json.loads(receipt.read_text())["bam_sha256"]:
        raise ValueError("Regional BAM checksum mismatch")
    if args.mode == "primary_only":
        bam = primary_alignment(bam, checksum)
    reference = destination / "reference-GRCh38"
    manifest, models = load_reference(reference)
    inventory = destination / "inventory-GRCh38-validated.json"
    record = next(r for r in json.loads(inventory.read_text())["variants"] if r["variant_id"] == args.variant_id)
    genome = reference_genome(reference, output / "reference-cache")
    contig = record["chrom"].removeprefix("chr")
    variant = Variant("MT" if contig == "M" else contig, record["pos"], record["ref"], record["alt"], ensembl=genome)
    expectations = {tid: apply_variant(record, models[tid]) for tid in manifest["variant_transcripts"][args.variant_id]}
    root = Path(isovar.__file__).resolve().parent.parent
    identity = dict(
        source_id=args.source_id, variant_id=args.variant_id, mode=args.mode,
        source_bam_sha256=checksum, input_bam_sha256=digest(bam), input_index_sha256=digest(str(bam) + ".bai"),
        receipt_sha256=digest(receipt), inventory_sha256=digest(inventory),
        reference_manifest_sha256=digest(reference / "manifest.json"),
        benchmark_sha256=digest(__file__),
        source_files={str(p.relative_to(root)): digest(p) for directory in
                      (root / "isovar", root / "tests/data/osteosarc/expansion")
                      for p in sorted(directory.glob("*.py"))},
        isovar=__version__, python=platform.python_version(), platform=platform.platform(),
        profile=args.profile, timeout_seconds_per_phase=args.timeout,
    )
    write_json(output / "identity.json", identity)
    trace, measurements = {}, {}
    collector = BenchmarkCollector(trace, merge_overlapping_fragments=True)
    creator = BenchmarkCreator(trace, capture_all_ranked=True)
    result = dict(status="running", trace=trace)
    try:
        with measure("production", output, measurements, args.timeout, args.profile):
            with pysam.AlignmentFile(bam) as handle:
                public_result, = run_isovar([variant], handle, transcript_id_whitelist=set(expectations),
                                           read_collector=collector, protein_sequence_creator=creator)
        result["public_protein_count"] = len(public_result.sorted_protein_sequences)
        result["rna_sha256"] = rna_fingerprint(getattr(creator, "rna_sequences", []))
        checked = []
        trace["stage"] = "independent_protein_validation"
        with measure("validation", output, measurements, args.timeout, args.profile):
            for protein in creator.uncapped_ranked or []:
                checked.append(protein_check(protein, expectations, creator.protein_sequence_length,
                                             creator.protein_context_peptide_length))
                result["checked_protein_count"] = len(checked)
        result.update(status="ok", proteins_sha256=sha256(json.dumps(checked, sort_keys=True).encode()).hexdigest(),
                      validation_status="ok" if all(p["validation_status"] == "ok" for p in checked) else "error")
        with gzip.open(output / "proteins.json.gz", "wt") as handle:
            json.dump(checked, handle, sort_keys=True)
    except AuditTimeout as error:
        result.update(status="resource_limit", error=str(error))
    except Exception as error:
        result.update(status="error", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        result["counts"] = counts_from_evidence(collector.evidence)
        write_json(output / "result.json", result)
    print(json.dumps(dict(result, measurements=measurements), sort_keys=True), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", required=True)
    parser.add_argument("--source-id", required=True)
    parser.add_argument("--variant-id", default="MT_ND5-chrM-12994")
    parser.add_argument("--mode", choices=["defaults", "primary_only"], default="defaults")
    parser.add_argument("--output", required=True)
    parser.add_argument("--timeout", type=float, default=180)
    parser.add_argument("--profile", action="store_true")
    run(parser.parse_args())
