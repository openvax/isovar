"""Audit full downloaded regions without changing Isovar's behavior.

Compare independent primary CIGAR counts, Isovar on primary-only BAMs, and
Isovar defaults on original BAMs. Exceptions are per-locus error results,
never successful zero counts. Protein validation uses six pinned Ensembl
transcripts and an independent sequence oracle, not vaccine constructs.

Usage: python tests/data/osteosarc/audit.py SOURCE_DIR NEW_REPORT.json
Acquire SOURCE_DIR with fetch_sources.py. Source checksums must match the
fixture manifest; this command does not download or edit original reads.
"""

import argparse
from collections import Counter
from hashlib import sha256
from importlib.metadata import version
import json
import logging
from pathlib import Path
import platform
import sys
import tempfile
from types import SimpleNamespace

import pysam
from varcode import Variant

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from isovar import __version__  # noqa: E402
from isovar.read_collector import ReadCollector  # noqa: E402
from isovar.protein_sequence_creator import ProteinSequenceCreator  # noqa: E402
from isovar.reference_context_helpers import reference_contexts_for_variant  # noqa: E402
from isovar.variant_sequence_creator import VariantSequenceCreator  # noqa: E402
from tests.real_rna_helpers import cigar_observation  # noqa: E402
from tests.osteosarc_protein_helpers import (  # noqa: E402
    REFERENCE, check_translation, expected_variant, load_references, reference_genome,
)

HERE = Path(__file__).resolve().parent
PRIMARY_EXCLUDE_FLAGS = 256 | 1024 | 2048


def error_result(error, stage):
    return {"status": "error", "stage": stage,
            "error": {"type": type(error).__name__, "message": str(error)}}


def sequence_result(variant, alt_reads, coverage):
    try:
        sequences = VariantSequenceCreator(
            min_variant_sequence_coverage=coverage).reads_to_variant_sequences(variant, alt_reads)
    except Exception as error:
        return error_result(error, "sequence_creation")
    return {"status": "ok", "min_coverage": coverage, "candidate_count": len(sequences),
            "represented_read_names": len({n for s in sequences for n in s.read_names}),
            "lengths": sorted({len(s) for s in sequences})}


def protein_result(variant, evidence, expected, coverage=2):
    """Exercise the public protein pipeline; retain failures at their stage."""
    creator = ProteinSequenceCreator(min_variant_sequence_coverage=coverage)
    whitelist = {expected["transcript_id"]}
    result = {"status": "ok", "min_coverage": coverage,
              "transcript_version": expected["transcript_version"],
              "protein_version": expected["protein_version"],
              "strand": expected["strand"],
              "reference_protein_length": len(expected["protein"]),
              "expected_mutant_protein_length": len(expected["mutant_protein"]),
              "expected_mutant_protein_sha256": sha256(expected["mutant_protein"].encode()).hexdigest()}
    stage = "protein_sequence_creation"
    try:
        sequences = creator._variant_sequence_creator.reads_to_variant_sequences(variant, evidence.alt_reads)
        result["rna_candidate_count"] = len(sequences)
        stage = "reference_matching"
        contexts = reference_contexts_for_variant(variant, creator._cdna_sequence_length, whitelist)
        result["reference_context_count"] = len(contexts)
        stage = "translation"
        translations = creator.translate_variant_reads(variant, evidence.alt_reads, whitelist)
        result["translation_count"] = len(translations)
        stage = "protein_validation"
        checks = [check_translation(t, expected, creator.protein_sequence_length) for t in translations]
        result["expected_matching_translation_count"] = sum(c["matches_expected"] for c in checks)
        stage = "protein_ranking"
        proteins = creator.sorted_protein_sequences_for_variant(variant, evidence, whitelist)
        stage = "protein_validation"
        result["top_proteins"] = []
        for protein in proteins:
            top_checks = [check_translation(t, expected, creator.protein_sequence_length)
                          for t in protein.translations]
            result["top_proteins"].append({
                "amino_acids": protein.amino_acids,
                "mutation_interval": [protein.mutation_start_idx, protein.mutation_end_idx],
                "contains_mutation": protein.contains_mutation, "frameshift": protein.frameshift,
                "ends_with_stop_codon": protein.ends_with_stop_codon,
                "supporting_read_names": protein.num_supporting_fragments,
                "matches_expected": all(c["matches_expected"] for c in top_checks),
                "expected_amino_acids": sorted({c["expected_amino_acids"] for c in top_checks}),
                "protein_starts_1based": sorted({c["protein_start_1based"] for c in top_checks}),
            })
        result["outcome"] = (
            "no_alt_reads" if not evidence.alt_reads else
            "no_rna_candidate" if not sequences else
            "no_reference_context" if not contexts else
            "no_translation" if not translations else
            "no_ranked_protein" if not proteins else
            "expected_top_protein" if all(p["matches_expected"] for p in result["top_proteins"])
            else "different_top_protein")
    except Exception as error:
        result.update(error_result(error, stage))
    return result


def isovar_result(bam, variant, expected=None):
    """Isolate locus failures; never skip an offending read within Isovar."""
    try:
        evidence = ReadCollector().read_evidence_for_variant(variant, bam)
    except Exception as error:
        return dict(error_result(error, "read_evidence"), counts=None)
    result = {
        "status": "ok",
        "counts": {"ref": len(evidence.ref_reads), "alt": len(evidence.alt_reads),
                   "other": len(evidence.other_reads)},
        "alt_read_names": len(evidence.alt_read_names),
        "sequence_creation_default": sequence_result(
            variant, evidence.alt_reads, VariantSequenceCreator().min_variant_sequence_coverage),
        # Sensitivity diagnostic, NOT a changed production default.
        "sequence_creation_coverage1": sequence_result(variant, evidence.alt_reads, 1),
    }
    if expected is not None:
        result["protein_default"] = protein_result(variant, evidence, expected)
        result["protein_coverage1"] = protein_result(variant, evidence, expected, coverage=1)
    return result


def independent_result(reads, record):
    """Direct anchored CIGAR evidence, not normalized or molecular counts."""
    min_mapq = ReadCollector().min_mapping_quality
    observations = [o for r in reads if not r.flag & PRIMARY_EXCLUDE_FLAGS
                    and r.mapping_quality >= min_mapq
                    and (o := cigar_observation(r, record)) is not None]
    ref, alt = record["ref"], record["alt"]
    if len(ref) != len(alt):
        # These six pinned cases contain anchored simple deletions only.
        assert len(alt) == 1 and ref.startswith(alt)
        ref, alt = ref[1:], alt[1:]
    counts = Counter(o["allele"] for o in observations)
    return {"counts": {"ref": counts[ref], "alt": counts[alt],
                       "other": len(observations) - counts[ref] - counts[alt]},
            "alt_quality_retention": {str(q): sum(
                o["allele"] == alt and o["qualities"] is not None and min(o["qualities"]) >= q
                for o in observations) for q in (0, 10, 20, 30)}}


def audit_alignment(name, full_bam, primary_bam, variants, genome=None, references=None):
    """Analyze an indexed source and its independently filtered primary copy."""
    rows = []
    with pysam.AlignmentFile(full_bam) as bam:
        reads = list(bam)
    for record in variants:
        record = dict(record, pos=int(record["pos"]))
        variant = SimpleNamespace(contig=record["chrom"].removeprefix("chr"),
                                  start=record["pos"], ref=record["ref"], alt=record["alt"],
                                  gene_names=(), short_description=record["variant_id"])
        expected = None
        if genome is not None:
            variant = Variant(variant.contig, variant.start, variant.ref, variant.alt, ensembl=genome)
            expected = expected_variant(record, references[record["gene"]])
        row = {"sample": name, "variant": record, "independent": independent_result(reads, record)}
        for mode, path in (("primary_only", primary_bam), ("defaults", full_bam)):
            with pysam.AlignmentFile(path) as bam:
                row[mode] = isovar_result(bam, variant, expected)
        rows.append(row)
    return rows


def build_report(source_dir):
    selection = json.loads((HERE / "selection.json").read_text())
    manifest = json.loads((HERE / "manifest.json").read_text())
    # Fail before the audit if input files drift; a source error is not a
    # biological zero or an Isovar allele-classification failure.
    for name, data in manifest["datasets"].items():
        source = source_dir / (name + ".sam")
        if sha256(source.read_bytes()).hexdigest() != data["source_sam_sha256"]:
            raise ValueError(f"Original regional source checksum mismatch: {source.name}")
    report = {
        "schema_version": 2,
        "software": {"isovar": __version__, "python": platform.python_version(),
                     "pysam": pysam.__version__, "samtools": pysam.__samtools_version__,
                     "pyensembl": version("pyensembl"), "varcode": version("varcode")},
        "read_collection_defaults": vars(ReadCollector()),
        "sequence_creation_defaults": vars(VariantSequenceCreator()),
        "protein_validation": {
            "scope": "one pinned transcript per variant; nominated edit only; not vaccine construct validation",
            "reference_manifest": json.loads((REFERENCE / "protein_reference_manifest.json").read_text()),
            "defaults": {k: v for k, v in vars(ProteinSequenceCreator()).items() if not k.startswith("_")},
        },
        "primary_exclude_flags": PRIMARY_EXCLUDE_FLAGS,
        "sources": {n: {"url": selection["sources"][n]["url"],
                        "sha256": d["source_sam_sha256"], "region_records": d["source_region_records"]}
                    for n, d in manifest["datasets"].items()},
        "rows": [],
    }
    with tempfile.TemporaryDirectory(prefix="isovar-sample-audit-") as directory:
        references = load_references()
        genome = reference_genome(Path(directory) / "reference")
        report["protein_validation"]["reference_dataset_identity"] = genome.reference_name
        for name in selection["sources"]:
            full = Path(directory) / (name + ".bam")
            primary = Path(directory) / (name + ".primary.bam")
            pysam.sort("--no-PG", "-o", str(full), str(source_dir / (name + ".sam")))
            pysam.index(str(full))
            pysam.view("--no-PG", "-b", "-F", str(PRIMARY_EXCLUDE_FLAGS), "-o",
                       str(primary), str(full), catch_stdout=False)
            pysam.index(str(primary))
            report["rows"].extend(audit_alignment(name, full, primary, selection["variants"], genome, references))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output already exists; choose a new report path")
    logging.disable(logging.CRITICAL)
    report = build_report(args.source_dir)
    with args.output.open("x") as handle:
        handle.write(json.dumps(report, indent=2) + "\n")
    print(f"Wrote {len(report['rows'])} variant/sample results to {args.output}")
