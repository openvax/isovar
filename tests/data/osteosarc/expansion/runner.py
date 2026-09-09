"""Run original regional RNA through Isovar and independent protein checks.

The output is one checkpoint per variant/source; exceptions cannot masquerade
as zero support. No source acquisition, reference padding or threshold changes.
"""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
from hashlib import sha256
from importlib.metadata import PackageNotFoundError, version
import json
import logging
import math
from pathlib import Path
import platform
import signal
import sys

import pysam
from varcode import Variant

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from isovar import __version__, run_isovar  # noqa: E402
from isovar.read_collector import ReadCollector  # noqa: E402
from isovar.protein_sequence_creator import ProteinSequenceCreator  # noqa: E402
from isovar.protein_sequence_helpers import mutant_peptide_window_count  # noqa: E402
from tests.data.osteosarc.expansion.acquire import error_record  # noqa: E402
from tests.data.osteosarc.expansion.inventory import digest, write_json  # noqa: E402
from tests.data.osteosarc.expansion.references import (  # noqa: E402
    apply_variant, check_translation, load_reference, reference_genome,
)
from tests.real_rna_helpers import cigar_observation  # noqa: E402


PRIMARY_EXCLUDE_FLAGS = 256 | 1024 | 2048


def installed_version(package):
    """Record unavailable distribution metadata without inventing a version."""
    try:
        return version(package)
    except PackageNotFoundError:
        return None


class AuditTimeout(TimeoutError):
    """The audit budget expired; this says nothing about biological support."""


@contextmanager
def time_limit(seconds):
    """Bound a Unix main-thread audit, without replacing an existing timer."""
    if seconds is None:
        yield
        return
    if not math.isfinite(seconds) or seconds <= 0:
        raise ValueError("Audit time limit must be finite and positive")
    if signal.getitimer(signal.ITIMER_REAL) != (0.0, 0.0):
        raise ValueError("An existing alarm prevents safely timing this audit")
    previous = signal.getsignal(signal.SIGALRM)

    def expire(signum, frame):
        raise AuditTimeout(f"Audit mode exceeded {seconds:g} seconds; #227 tracks high-depth matching")

    signal.signal(signal.SIGALRM, expire)
    try:
        signal.setitimer(signal.ITIMER_REAL, seconds)
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


class TracedCollector(ReadCollector):
    """Observe stage/witness information without skipping or altering reads."""

    def __init__(self, trace, **options):
        super().__init__(**options)
        self.trace = trace
        self.evidence = None

    def read_evidence_for_variant(self, *args, **kwargs):
        self.trace["stage"] = "read_collection"
        self.evidence = super().read_evidence_for_variant(*args, **kwargs)
        return self.evidence

    def locus_read_from_pysam_aligned_segment(self, pysam_aligned_segment, *args, **kwargs):
        try:
            return super().locus_read_from_pysam_aligned_segment(pysam_aligned_segment, *args, **kwargs)
        except Exception:
            self.trace["witness_sam"] = pysam_aligned_segment.to_string()
            raise


class TracedCreator(ProteinSequenceCreator):
    """Record actual intermediate counts on the public production path."""

    def __init__(self, trace, capture_all_ranked=False):
        super().__init__()
        self.trace = trace
        self.capture_all_ranked = capture_all_ranked
        self.uncapped_ranked = None

    def sorted_protein_sequences_for_variant(self, *args, **kwargs):
        if not self.capture_all_ranked:
            return super().sorted_protein_sequences_for_variant(*args, **kwargs)
        # The production method computes/sorts the full set before its final
        # slice. Observe that set once, then restore exactly the default cap.
        # This is audit instrumentation, not an extra assembly/translation run.
        cap = self.max_protein_sequences_per_variant
        try:
            self.max_protein_sequences_per_variant = None
            self.uncapped_ranked = super().sorted_protein_sequences_for_variant(*args, **kwargs)
        finally:
            self.max_protein_sequences_per_variant = cap
        return self.uncapped_ranked[:cap] if cap else self.uncapped_ranked

    def variant_sequences_from_reads(self, *args, **kwargs):
        self.trace["stage"] = "rna_sequence_creation"
        sequences = super().variant_sequences_from_reads(*args, **kwargs)
        self.trace["rna_candidate_count"] = len(sequences)
        self.trace["stage"] = "reference_matching"
        return sequences

    def all_pairs_translations(self, variant_sequences, reference_contexts):
        self.trace["reference_context_count"] = len(reference_contexts)
        self.trace["stage"] = "translation"
        result = super().all_pairs_translations(variant_sequences, reference_contexts)
        self.trace["translation_count"] = len(result)
        self.trace["stage"] = "protein_ranking"
        return result

    def translation_from_variant_sequence_and_reference_context(self, *args, **kwargs):
        result = super().translation_from_variant_sequence_and_reference_context(*args, **kwargs)
        key = ("reference_matching_rejected" if result is None else
               "translations_without_mutation" if not result.contains_mutation else "translations_with_mutation")
        self.trace[key] = self.trace.get(key, 0) + 1
        return result


def counts_from_evidence(evidence):
    if evidence is None:
        return None
    groups = {name: getattr(evidence, name + "_reads") for name in ("ref", "alt", "other")}
    return dict(
        reads={name: sum(r.source_read_count for r in reads) for name, reads in groups.items()},
        read_objects={name: len(reads) for name, reads in groups.items()},
        template_names={name: len({r.name for r in reads}) for name, reads in groups.items()},
        all_template_names=len({r.name for reads in groups.values() for r in reads}),
    )


def protein_check(protein, expectations, length, peptide_length):
    checks = []
    for translation in protein.translations:
        for transcript in translation.reference_context.transcripts:
            if transcript.id not in expectations:
                raise ValueError(f"Unpinned translated transcript: {transcript.id}")
            try:
                checks.append(dict(status="ok", **check_translation(translation, expectations[transcript.id], length)))
            except (AssertionError, ValueError) as error:
                checks.append(dict(error_record(error, "independent_protein_validation"),
                                   transcript_id=transcript.id, matches_expected=None))
    if not checks:
        raise ValueError("No independent transcript expectation for returned protein")
    return dict(
        amino_acids=protein.amino_acids, length=len(protein.amino_acids),
        mutation_interval=[protein.mutation_start_idx, protein.mutation_end_idx],
        contains_mutation=protein.contains_mutation, frameshift=protein.frameshift,
        ends_with_stop_codon=protein.ends_with_stop_codon,
        supporting_reads=protein.num_supporting_reads, supporting_template_names=protein.num_supporting_fragments,
        supporting_read_names=sorted(protein.read_names_supporting_protein_sequence),
        retained_cdna_min_coverages=sorted({int(t.untrimmed_variant_sequence.min_coverage()) for t in protein.translations}),
        mutation_containing_peptide_windows=mutant_peptide_window_count(protein, peptide_length),
        checks=checks, validation_status="ok" if all(c["status"] == "ok" for c in checks) else "error",
        matches_expected=all(c["matches_expected"] is True for c in checks),
    )


def audit_mode(path, variant, expectations, timeout=None, capture_all_ranked=False):
    """Run the real API, retaining counts even if a later stage raises."""
    trace = {}
    collector = TracedCollector(trace, merge_overlapping_fragments=True)
    creator = TracedCreator(trace, capture_all_ranked=capture_all_ranked)
    try:
        with time_limit(timeout):
            with pysam.AlignmentFile(path) as bam:
                result, = run_isovar([variant], bam, transcript_id_whitelist=set(expectations),
                                      read_collector=collector, protein_sequence_creator=creator)
            trace["stage"] = "independent_protein_validation"
            observed = creator.uncapped_ranked if capture_all_ranked else result.sorted_protein_sequences
            if observed is None:
                observed = []
            checked = [protein_check(p, expectations, creator.protein_sequence_length,
                                      creator.protein_context_peptide_length) for p in observed]
            proteins = checked[:len(result.sorted_protein_sequences)]
        counts = counts_from_evidence(collector.evidence)
        outcome = (
            "no_callable_allele" if sum(counts["reads"].values()) == 0 else
            "no_alt_reads" if not counts["reads"]["alt"] else
            "missing_independent_reference_model" if not expectations else
            "no_rna_candidate" if not trace.get("rna_candidate_count") else
            "no_reference_context" if not trace.get("reference_context_count") else
            "no_translation" if not trace.get("translation_count") else
            "no_ranked_protein" if not proteins else
            "protein_validation_error" if any(p["validation_status"] != "ok" for p in proteins) else
            "expected_top_protein" if proteins[0]["matches_expected"] else
            "different_top_protein")
        audited = dict(status="ok", outcome=outcome, counts=counts, stages=trace,
                       passes_all_filters=result.passes_all_filters,
                       filter_values=result.filter_values, proteins=proteins,
                       all_ranked_match_expected=all(p["matches_expected"] for p in proteins) if proteins else None)
        if capture_all_ranked:
            audited.update(uncapped_ranked_proteins=checked,
                           uncapped_validation_status="ok" if all(p["validation_status"] == "ok" for p in checked) else "error",
                           uncapped_all_match_expected=all(p["matches_expected"] for p in checked) if checked else None,
                           returned_protein_limit=creator.max_protein_sequences_per_variant)
        return audited
    except AuditTimeout as error:
        result = error_record(error, trace.get("stage", "run_isovar"))
        result.update(status="resource_limit", outcome="audit_timeout", time_limit_seconds=timeout,
                      counts=counts_from_evidence(collector.evidence), stages=trace)
        return result
    except Exception as error:
        return dict(error_record(error, trace.get("stage", "run_isovar")),
                    counts=counts_from_evidence(collector.evidence), stages=trace)


def independent_counts(bam, record):
    """Descriptive exact-CIGAR quality and barcode counts, not calibrated scores."""
    from tests.data.osteosarc.expansion.references import minimal_edit

    _, ref, alt = minimal_edit(record)
    counts = Counter()
    tags = {group: {"cells": set(), "umis": set()} for group in ("ref", "alt", "other")}
    quality = {str(q): 0 for q in (0, 10, 20, 30)}
    unknown_quality = 0
    total, spliced_out, filtered = 0, 0, 0
    mapqs = Counter()
    region = (record["chrom"], max(0, record["pos"] - 2), record["pos"] + len(record["ref"]))
    for read in bam.fetch(*region):
        total += 1
        if read.get_overlap(region[1], region[2]) == 0:
            spliced_out += 1
            continue
        if read.flag & PRIMARY_EXCLUDE_FLAGS or read.mapping_quality < 1:
            filtered += 1
            continue
        observation = cigar_observation(read, record)
        if observation is None:
            continue
        group = "ref" if observation["allele"] == ref else "alt" if observation["allele"] == alt else "other"
        counts[group] += 1
        rg = read.get_tag("RG") if read.has_tag("RG") else None
        if read.has_tag("CB"):
            cell = (rg, read.get_tag("CB"))
            tags[group]["cells"].add(cell)
            if read.has_tag("UB"):
                tags[group]["umis"].add((*cell, read.get_tag("UB")))
        if group == "alt":
            mapqs[str(read.mapping_quality)] += 1
            if observation["qualities"] is None:
                unknown_quality += 1
            else:
                for cutoff in quality:
                    quality[cutoff] += min(observation["qualities"]) >= int(cutoff)
    return dict(
        region_alignments=total, no_aligned_base_in_expanded_interval=spliced_out,
        primary_or_mapq_filtered=filtered, counts={g: counts[g] for g in ("ref", "alt", "other")},
        alt_min_query_quality_retention=quality, alt_missing_quality=unknown_quality,
        alt_mapq_histogram=dict(sorted(mapqs.items())),
        tag_counts={g: {"distinct_rg_cb": len(t["cells"]), "distinct_rg_cb_ub": len(t["umis"])} for g, t in tags.items()},
    )


def audit_identity(inventory_path, reference_dir, timeout=None):
    """Make stale checkpoints fail closed when code, inputs or settings change."""
    paths = [*sorted((ROOT / "isovar").rglob("*.py")), Path(__file__),
             Path(__file__).with_name("references.py"), ROOT / "tests/real_rna_helpers.py",
             ROOT / "tests/osteosarc_protein_helpers.py"]
    return dict(
        sources={str(p.relative_to(ROOT)): digest(p) for p in paths},
        inventory_sha256=digest(inventory_path), reference_manifest_sha256=digest(reference_dir / "manifest.json"),
        software=dict(isovar=__version__, python=platform.python_version(), pysam=pysam.__version__,
                      dependencies={name: installed_version(name) for name in ("varcode", "pyensembl", "numpy", "biopython")}),
        time_limit_seconds_per_mode=timeout,
        ranking_capture="observe all pre-cap ranked proteins; return exactly the configured public cap",
        settings=dict(read_collector=vars(ReadCollector(merge_overlapping_fragments=True)),
                      protein_sequence_creator={k: v for k, v in vars(ProteinSequenceCreator()).items()
                                                if not k.startswith("_")}),
        modes=dict(defaults="public API defaults, validated transcript whitelist",
                   primary_only=dict(exclude_flags=PRIMARY_EXCLUDE_FLAGS)),
    )


def primary_alignment(path, source_checksum):
    """Generate a separately checksummed flag-filtered derivative."""
    output = path.with_name("primary-v1-" + path.name)
    receipt = output.with_suffix(".json")
    request = dict(source_bam_sha256=source_checksum, exclude_flags=PRIMARY_EXCLUDE_FLAGS)
    if receipt.exists():
        saved = json.loads(receipt.read_text())
        if saved["request"] != request or saved["bam_sha256"] != digest(output):
            raise ValueError("Primary derivative checksum/request mismatch")
        if saved["index_sha256"] != digest(str(output) + ".bai"):
            raise ValueError("Primary derivative index checksum mismatch")
    else:
        if output.exists():
            raise ValueError("Unreceipted primary derivative; inspect before reuse")
        pysam.view("--no-PG", "-b", "-F", str(PRIMARY_EXCLUDE_FLAGS), "-o", str(output), str(path), catch_stdout=False)
        pysam.index(str(output))
        write_json(receipt, dict(request=request, bam_sha256=digest(output), index_sha256=digest(str(output) + ".bai")))
    return output


def run_source(destination, source_id, reference_dir, output_dir, inventory_name="inventory.json", acquisition_label=None,
               timeout=None):
    logging.disable(logging.CRITICAL)
    destination, reference_dir, output_dir = map(Path, (destination, reference_dir, output_dir))
    inventory_path = destination / inventory_name
    inventory = json.loads(inventory_path.read_text())
    source = next(s for s in inventory["alignments"] if s["source_id"] == source_id)
    directory = destination / "alignments" / source_id
    manifest, models = load_reference(reference_dir)
    assembly = manifest["assembly"]
    label = acquisition_label or assembly
    receipt_path = directory / f"regions-{label}.json"
    if not receipt_path.exists():
        return dict(source_id=source_id, status="source_not_acquired")
    receipt = json.loads(receipt_path.read_text())
    if receipt["status"] != "ok":
        return dict(source_id=source_id, status="source_acquisition_error")
    bam_path = directory / f"regions-{label}.bam"
    if digest(bam_path) != receipt["bam_sha256"]:
        raise ValueError("Regional source checksum mismatch")
    primary = primary_alignment(bam_path, receipt["bam_sha256"])
    genome = reference_genome(reference_dir, output_dir / "reference-cache")
    rows = output_dir / source_id
    rows.mkdir(parents=True, exist_ok=True)
    identity = audit_identity(inventory_path, reference_dir, timeout)
    identity["source_bam_sha256"] = receipt["bam_sha256"]
    identity["acquisition_receipt_sha256"] = digest(receipt_path)
    identity_path = rows / "run.json"
    if identity_path.exists():
        if json.loads(identity_path.read_text()) != identity:
            raise ValueError("Audit checkpoint identity changed; choose a new output directory")
    else:
        if list(rows.glob("*.json")):
            raise ValueError("Unversioned audit checkpoints; choose a new output directory")
        write_json(identity_path, identity)
    with pysam.AlignmentFile(bam_path) as bam:
        names = set(bam.references)
    requested_ids = receipt.get("variant_ids", [v["variant_id"] for v in inventory["variants"]])
    records = [v for v in inventory["variants"] if v["variant_id"] in requested_ids]
    if len(records) != len(requested_ids):
        raise ValueError("Acquisition contains unknown/duplicate variant identities")
    for record in records:
        output = rows / (record["variant_id"] + ".json")
        if output.exists():
            continue
        contig = record["chrom"].removeprefix("chr")
        available = names & ({"M", "MT", "chrM", "chrMT"} if contig in ("M", "MT") else {contig, "chr" + contig})
        row = dict(source_id=source_id, variant_id=record["variant_id"], source_bam_sha256=receipt["bam_sha256"],
                   software=dict(isovar=__version__, python=platform.python_version(), pysam=pysam.__version__))
        if len(available) != 1:
            row.update(status="missing_or_ambiguous_contig", counts=None)
        else:
            contig = "MT" if contig in ("M", "MT") else contig
            stage = "variant_reference_setup"
            try:
                variant = Variant(contig, record["pos"], record["ref"], record["alt"], ensembl=genome)
                row["native_variant"] = {k: record[k] for k in ("assembly", "chrom", "pos", "ref", "alt")}
                row["normalized_variant"] = dict(contig=variant.contig, start=variant.start, ref=variant.ref, alt=variant.alt)
                expectations = {tid: apply_variant(record, models[tid]) for tid in manifest["variant_transcripts"][record["variant_id"]]}
                row["reference_transcripts"] = sorted(expectations)
                stage = "independent_cigar_counting"
                with pysam.AlignmentFile(bam_path) as bam:
                    row["independent_primary"] = independent_counts(bam, dict(record, chrom=next(iter(available))))
                row["defaults"] = audit_mode(bam_path, variant, expectations, timeout, capture_all_ranked=True)
                row["primary_only"] = audit_mode(primary, variant, expectations, timeout, capture_all_ranked=True)
                if row["independent_primary"]["region_alignments"] == 0:
                    for mode in ("defaults", "primary_only"):
                        if row[mode]["status"] == "ok":
                            row[mode]["outcome"] = "no_overlapping_alignment"
                if contig == "MT":
                    row["mitochondrial_interpretation"] = dict(
                        translation_compartment="mitochondrial, conditional on read origin",
                        numt_origin_status="not_ruled_out", dna_heteroplasmy=None, cancer_cell_fraction=None,
                        allele_fraction_scope="RNA alignment/template observations, not diploid clonality",
                    )
                row["status"] = "audited"
            except Exception as error:
                row.update(error_record(error, stage), counts=None)
        write_json(output, row)
        print(source_id, record["variant_id"], row.get("defaults", {}).get("outcome", row["status"]), flush=True)
    return dict(source_id=source_id, status="audited", variants=len(records))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--source-id")
    parser.add_argument("--inventory", default="inventory.json")
    parser.add_argument("--acquisition-label")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--time-limit", type=float, default=120,
                        help="Seconds per source/locus/mode; expiration is a reported failure, never zero support")
    args = parser.parse_args()
    if not __debug__:
        parser.error("Independent validation assertions require normal Python mode (no -O)")
    if not 1 <= args.workers <= 4:
        parser.error("Use 1 through 4 audit workers")
    if not math.isfinite(args.time_limit) or not 0 < args.time_limit <= 1800:
        parser.error("Use a finite time limit greater than zero and at most 1800 seconds")
    args.output.mkdir(parents=True, exist_ok=True)
    logging.disable(logging.CRITICAL)
    reference_genome(args.reference, args.output / "reference-cache")
    if args.source_id:
        print(run_source(args.destination, args.source_id, args.reference, args.output, args.inventory,
                         args.acquisition_label, args.time_limit))
    else:
        inventory = json.loads((args.destination / args.inventory).read_text())
        sources = [r for r in inventory["alignments"] if r["scope"] == "rna_candidate"]
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = [executor.submit(run_source, args.destination, s["source_id"], args.reference, args.output,
                                       args.inventory, args.acquisition_label, args.time_limit) for s in sources]
            for future in as_completed(futures):
                print(future.result(), flush=True)
