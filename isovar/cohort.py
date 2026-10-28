"""Reproducible, source-separated cohort reconstruction using the shared APIs."""

import gzip
import csv
import json
import logging
import os
from collections import Counter
from importlib.metadata import version
from pathlib import Path
import sqlite3
import tempfile
import traceback

import pysam
from pyensembl import EnsemblRelease
from varcode import Variant
from varcode.reference import infer_genome

from .contig_resolution import UnresolvedContigError, require_resolved_contig
from .evidence_report import (
    catalog_coding_index, changed_reference_interval, direct_read_context_proteins,
    overlapping_catalog_orfs, overlaps_exon,
)
from .isovar_result import IsovarResult
from .known_variants import KnownVariants
from .normal_evidence import literal_allele_key
from .protein_hypotheses import export_protein_hypotheses
from .protein_reporting import protein_comparison_status
from .protein_sequence_creator import ProteinSequenceCreator
from .read_collector import ReadCollector
from .rna_evidence import canonical_json, content_identifier
from .sv_rna_audit import file_sha256, read_json, run_sv_rna_audit, write_json
from .variant_rna_orfs import reconstruct_variant_orfs


logger = logging.getLogger(__name__)
SCHEMA = "isovar.cohort_manifest.v1"


def _write_manifest(path, manifest):
    """Publish a complete manifest even when separate source jobs share output."""
    descriptor, temporary = tempfile.mkstemp(prefix=".manifest-", suffix=".json", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w") as stream:
            stream.write(canonical_json(manifest) + "\n")
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def reconstruction_identity():
    """Pin scientific code separately from presentation-only changes."""
    from . import __version__

    root = Path(__file__).parent
    files = {p.relative_to(root).as_posix(): file_sha256(p) for p in sorted(root.rglob("*.py"))
             if p.name != "cohort_pdf.py"}
    return dict(isovar_version=__version__, implementation_files=files,
                dependencies={name: version(name) for name in
                              ("pysam", "edlib", "pyensembl", "varcode", "orfdata")})


def _inputs(manifest, root, keys=None):
    paths = {}
    for name, descriptor in manifest["files"].items():
        if keys is not None and name not in keys:
            continue
        path = (root / descriptor["path"]).resolve()
        if file_sha256(path) != descriptor["sha256"]:
            raise ValueError("Cohort input checksum mismatch: " + name)
        paths[name] = path
    return paths


def _selected_file_keys(manifest, source_ids):
    """A source worker validates its own inputs, without reading other BAMs."""
    sources = {s["source_id"]: s for s in manifest["sources"]}
    keys = set()
    for family in ("small_variants", "sv_audits", "fusions", "viral", "isoforms", "prior_assessments"):
        for task in manifest.get(family, []):
            if source_ids and task.get("source_id") not in source_ids:
                continue
            source = sources[task["source_id"]]
            if source.get("unavailable_reason"):
                continue
            keys.update(task.get("reference_files", []))
            for name in ("variants", "germline_variants", "normal_evidence", "nc_orf_catalog", "manifest",
                         "input", "original_fixture", "locus_annotations", "genome_fasta", "targets"):
                if task.get(name):
                    keys.add(task[name])
            if family in ("small_variants", "viral", "isoforms"):
                keys.add(source["bam"])
    if keys - manifest["files"].keys():
        raise ValueError("Task refers to undeclared input files: " + ", ".join(sorted(keys - manifest["files"].keys())))
    return keys


def _small_variant_file_keys(task, source):
    """One input contract for computation, resume identity and report verification."""
    return [source["bam"], task["variants"]] + [task[key] for key in
        ("germline_variants", "normal_evidence", "nc_orf_catalog") if task.get(key)] + task.get("reference_files", [])


def _report_identity(task, source, manifest, family):
    """Scientific declarations, independent of input locations and presentation."""
    task = {k: v for k, v in task.items() if k not in {"title", "summary", "locus_annotations"}}
    source = {k: v for k, v in source.items() if k != "title"}
    if family == "small_variants":
        keys = _small_variant_file_keys(task, source)
    else:
        keys = _selected_file_keys(dict(files=manifest["files"], sources=[source], **{family: [task]}), None)
    return dict(task=task, source=source,
                inputs={key: manifest["files"][key]["sha256"] for key in keys})


def _check_report_identity(receipt, expected):
    recorded = receipt["identity"]
    actual = dict(task={k: v for k, v in recorded["task"].items()
                        if k not in {"title", "summary", "locus_annotations"}},
                  source={k: v for k, v in recorded["source"].items() if k != "title"},
                  inputs=recorded["inputs"])
    if actual != expected:
        raise ValueError("Recorded analysis identity differs from the manifest; "
                         "restore the original declarations or compute into a new output directory")


def _artifact_path(directory, name):
    path = (directory / name).resolve()
    if not path.is_relative_to(directory.resolve()):
        raise ValueError("Recorded artifact is outside its analysis directory")
    return path


def _render_small_variant_receipt(task, source, manifest, directory):
    path = directory / "receipt.json"
    if not path.exists():
        if (directory / "variants.jsonl.gz").exists():
            raise ValueError("Missing analysis receipt; recompute into a new output directory")
        return dict(status="not_run")
    receipt = read_json(path)
    if receipt.get("schema") != "isovar.cohort_source.v1":
        raise ValueError("Unsupported cohort source receipt")
    _check_report_identity(receipt, _report_identity(task, source, manifest, "small_variants"))
    if file_sha256(_artifact_path(directory, receipt["results"])) != receipt["result_sha256"]:
        raise ValueError("Recorded small-variant result digest mismatch")
    return receipt


def _additional_artifacts(family, directory, target):
    """Pin every scientific artifact later consumed by the renderer."""
    paths = [target]
    if family == "sv_audits":
        paths.extend(_artifact_path(directory, row["artifact"]) for row in read_json(target)["runs"])
    elif family == "viral":
        paths.extend(directory / name for name in ("sequences.fasta", "proteins.fasta"))
    return {str(path.relative_to(directory)): file_sha256(path) for path in paths}


def _check_additional_receipt(path, expected, directory, target, *, resume=False):
    if not path.exists():
        if target.exists():
            raise ValueError("Missing additional analysis receipt; recompute into a new output directory")
        return None
    receipt = read_json(path)
    if receipt.get("schema") != "isovar.cohort_additional.v1":
        raise ValueError("Unsupported additional analysis receipt")
    _check_report_identity(receipt, expected)
    if resume and receipt["identity"]["software"] != reconstruction_identity():
        raise ValueError("Resume implementation changed; compute into a new output directory")
    if str(target.relative_to(directory)) not in receipt["artifacts"]:
        raise ValueError("Analysis receipt does not identify the result")
    for name, digest in receipt["artifacts"].items():
        if file_sha256(_artifact_path(directory, name)) != digest:
            raise ValueError("Recorded additional artifact digest mismatch: " + name)
    return receipt


def _literal(row, genome):
    return Variant(str(row["contig"]), int(row["start"]), row["ref"], row["alt"], genome,
                   convert_ucsc_contig_names=True)


def _normal_background(task, source, files, nominations, variants, genome):
    """Use one declared background source and retain its per-allele normal scope."""
    if not task.get("normal_evidence"):
        germline = ([_literal(row, genome) for row in read_json(files[task["germline_variants"]])]
                    if task.get("germline_variants") else [])
        return nominations, germline, None
    if task.get("germline_variants"):
        raise ValueError("Use normal_evidence or germline_variants, not both; reconcile the normal inventory first")
    normal = read_json(files[task["normal_evidence"]])
    if normal.get("schema") != "isovar.matched_normal.v1":
        raise ValueError("Unsupported matched-normal evidence receipt")
    if normal["sample_id"] != source["sample_id"] or normal["reference_name"] != genome.reference_name:
        raise ValueError("Matched-normal sample/reference differs from the RNA cohort")
    index = {tuple(row["normalized_allele"]): row for row in normal["alleles"]}
    if len(index) != len(normal["alleles"]):
        raise ValueError("Duplicate normal allele assessment")
    updated = []
    for nomination, variant in zip(nominations, variants):
        key = literal_allele_key(variant.contig, variant.start, variant.ref, variant.alt)
        evidence = index.get(key, dict(state="unavailable", source_assessments=[],
                                       reason="allele_not_assessed_in_normal_receipt", origin_conflict=False))
        conflict = nomination.get("origin_conflict") or (
            "Matched-normal assessment: " + evidence["state"] if evidence["origin_conflict"] else False)
        updated.append(dict(nomination, matched_normal_evidence=evidence,
            normal_evidence_status="Matched-normal evidence: " + evidence["state"],
            origin_conflict=conflict))
    scope = dict(receipt_input=task["normal_evidence"],
        receipt_sha256=file_sha256(files[task["normal_evidence"]]),
        reference_name=normal["reference_name"], sample_id=normal["sample_id"],
        source_inventory=normal["source_inventory"],
        allele_states=dict(Counter(row["matched_normal_evidence"]["state"] for row in updated)),
        background_flank=normal["background_flank"], background_alleles=len(normal["background_variants"]),
        interpretation=normal["interpretation"])
    return updated, [_literal(row, genome) for row in normal["background_variants"]], scope


def _reference_declaration(source, genome):
    """Check an explicit assembly declaration; contig spelling is not evidence."""
    declared = source.get("reference_name") or source.get("provenance", {}).get("assembly")
    if not declared:
        return dict(status="not_independently_declared", annotation_reference=genome.reference_name)
    declared_genome, _ = infer_genome(declared)
    if declared_genome.reference_name != genome.reference_name:
        raise ValueError(f"RNA reference {declared!r} differs from annotation {genome.reference_name}; "
                         "contig alias conversion does not perform liftover")
    return dict(status="assembly_declaration_matches", supplied=declared,
                annotation_reference=genome.reference_name,
                sequence_equivalence="not_established_by_assembly_name")


def _nc_orfs(variant, evidence, catalog, index, creator, sample, source):
    """Bind every overlapping catalogue model; retain failures independently."""
    if catalog is None:
        return dict(status="not_assessed", reason="No ncORF catalogue supplied")
    begin, end, _ = changed_reference_interval(variant.original_start, variant.original_ref,
                                              variant.original_alt)
    records = [catalog.get(key) for key in overlapping_catalog_orfs(index, variant.contig, begin - 1, end)]
    if not evidence.alt_reads:
        return dict(status="no_alt_reads", overlapping_models=len(records), bound_models=None,
                    unresolved_models=[], result=None)
    references, unresolved = [], []
    for record in records:
        try:
            transcript = variant.genome.transcript_by_id(record.transcript_id.split(".")[0])
            references.append(record.bind(transcript))
        except (ValueError, KeyError) as error:
            unresolved.append(dict(orf_id=record.orf_id, transcript_id=record.transcript_id,
                                   source_ids=list(record.source_ids), reason=str(error)))
    result = reconstruct_variant_orfs(variant, evidence, references, sample_id=sample, source=source,
                                      protein_sequence_creator=creator) if references else None
    return dict(status="assessed", overlapping_models=len(records), bound_models=len(references),
                unresolved_models=unresolved, result=result)


def _evaluate(variant, nomination, bam, source, collector, creator, catalog, index, known, dense_limit):
    evidence = collector.read_evidence_for_variant(variant, bam)
    effects = variant.effects() if evidence.alt_reads else []
    effect = effects.top_priority_effect() if effects else None
    dense = len(evidence.alt_reads) > dense_limit
    products = (direct_read_context_proteins(variant, evidence.alt_reads, creator) if dense else
                creator.sorted_protein_sequences_for_variant(variant, evidence)) if evidence.alt_reads else []
    products = [protein.with_known_variants(known) for protein in products]
    result = IsovarResult(variant=variant, read_evidence=evidence, predicted_effect=effect,
                          sorted_protein_sequences=products, protein_sequence_settings=creator.settings())
    export = export_protein_hypotheses([result], sample_id=source["sample_id"], source=source["source_id"],
                                      alignment_header=bam.header.to_dict(),
                                      cell_umi_alignment_file=bam if source.get("cell_umi_labels") else None,
                                      read_collector=collector if source.get("cell_umi_labels") else None)
    contig = require_resolved_contig(collector.resolve_variant_contig(variant, bam), bam)
    raw = None
    if not (evidence.ref_reads or evidence.alt_reads or evidence.other_reads):
        start = variant.original_start - 1
        raw = bam.count(contig, start, start + max(1, len(variant.original_ref)))
    gene_names = sorted({getattr(e, "gene_name", None) for e in effects} - {None, ""})
    if not gene_names:
        gene_names = sorted({g.name for g in variant.genome.genes_at_locus(
            variant.contig, variant.original_start,
            variant.original_start + max(1, len(variant.original_ref)) - 1) if g.name})
    length_delta = len(variant.original_alt) - len(variant.original_ref)
    category = ("SNV" if len(variant.original_ref) == len(variant.original_alt) == 1 else
                "Frameshift" if any(p.frameshift for p in products) or type(effect).__name__ == "FrameShift" else
                "In-frame indel" if products and length_delta and length_delta % 3 == 0 else "Other non-SNV")
    nc_creator = (ProteinSequenceCreator(**dict(creator.settings(), variant_sequence_assembly=False))
                  if dense else creator)
    nc = _nc_orfs(variant, evidence, catalog, index, nc_creator, source["sample_id"], source["source_id"])
    nc["assembly_route"] = "direct_read_contexts_high_depth" if dense else "overlap_assembly"
    return dict(status="assessed", nomination=nomination, category=category, gene_names=gene_names,
                igv=f"{contig}:{variant.original_start}-"
                    f"{variant.original_start + max(1, len(variant.original_ref)) - 1}",
                raw_alignment_overlap_when_no_allele=raw,
                exon_overlap=overlaps_exon(variant.genome, variant.contig, variant.original_start,
                                          variant.original_ref, variant.original_alt),
                assembly_route="direct_read_contexts_high_depth" if dense else "overlap_assembly",
                protein_hypotheses=export, nc_orfs=nc)


def run_small_variant_source(task, source, files, output, *, resume=False):
    """Reconstruct every nominated allele with persistent per-call checkpoints."""
    import orfdata

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    genome = EnsemblRelease(task["ensembl_release"])
    reference_declaration = _reference_declaration(source, genome)
    nominations = read_json(files[task["variants"]])
    variants = [_literal(row, genome) for row in nominations]
    if len(set(variants)) != len(variants):
        raise ValueError("Duplicate literal variants: consolidate their nomination provenance first")
    nominations, germline, normal_scope = _normal_background(task, source, files, nominations, variants, genome)
    conversions = Counter((v.original_contig, v.contig) for v in variants + germline
                          if v.original_contig != v.contig)
    annotation_diagnostics = dict(
        reference_name=genome.reference_name, ensembl_release=task["ensembl_release"],
        reference_declaration=reference_declaration,
        contig_conversions=[dict(supplied=source_name, annotation=target, alleles=count)
                            for (source_name, target), count in sorted(conversions.items())],
        scope="Name conversion only; assembly and sequence equivalence are not inferred from labels")
    if conversions:
        logger.warning("%s: converted supplied contig names through Varcode: %s. "
                       "Original names are preserved; this is not coordinate liftover.",
                       source["source_id"], "; ".join(f"{a} -> {b} ({n} alleles)"
                                                    for (a, b), n in sorted(conversions.items())))
    known = KnownVariants(somatic=[v for v, n in zip(variants, nominations)
                                    if not n.get("origin_conflict")], germline=germline)
    collector = ReadCollector(**task.get("read_collection", {}))
    settings = dict(task.get("protein_settings", {}))
    settings.setdefault("max_protein_sequences_per_variant", 0)
    creator = ProteinSequenceCreator(**settings)
    catalog = orfdata.load_catalog(files[task["nc_orf_catalog"]]) if task.get("nc_orf_catalog") else None
    index = catalog_coding_index(catalog.records) if catalog is not None else None
    if catalog and catalog.genome_build != genome.reference_name:
        raise ValueError("ncORF catalogue and variant reference differ")
    bam_path = files[source["bam"]]
    indexes = [path for path in (Path(str(bam_path) + ".bai"), bam_path.with_suffix(".bai"),
                                Path(str(bam_path) + ".csi")) if path.exists()]
    if not indexes:
        raise ValueError("Indexed genomic RNA BAM required")
    identity = dict(task=task, source=source, software=reconstruction_identity(),
                    inputs={key: file_sha256(files[key]) for key in _small_variant_file_keys(task, source)},
                    indexes={path.name: file_sha256(path) for path in indexes})
    key = content_identifier("cohort_source", identity)
    checkpoint = output / "checkpoint.sqlite"
    if checkpoint.exists() and not resume:
        raise ValueError("Checkpoint exists; use --resume or a new output directory")
    with sqlite3.connect(checkpoint) as db:
        db.execute("PRAGMA journal_mode=WAL")
        db.execute("CREATE TABLE IF NOT EXISTS metadata (identity TEXT NOT NULL)")
        db.execute("CREATE TABLE IF NOT EXISTS calls (ordinal INTEGER PRIMARY KEY, result TEXT NOT NULL)")
        prior = db.execute("SELECT identity FROM metadata").fetchall()
        if prior and prior != [(key,)]:
            raise ValueError("Resume inputs, parameters or implementation changed")
        if not prior:
            db.execute("INSERT INTO metadata VALUES (?)", (key,))
            db.commit()
        completed = {row[0] for row in db.execute("SELECT ordinal FROM calls")}
        with pysam.AlignmentFile(bam_path) as bam:
            if not bam.has_index():
                raise ValueError("RNA alignment requires an index")
            for i, (variant, nomination) in enumerate(zip(variants, nominations)):
                if i in completed:
                    continue
                try:
                    row = _evaluate(variant, nomination, bam, source, collector, creator, catalog, index, known,
                                    task.get("max_alt_observations_for_overlap_assembly", 750))
                except UnresolvedContigError as error:
                    row = dict(status="unresolved_chromosome", nomination=nomination,
                               message=str(error), contig_resolution=error.to_dict())
                except Exception as error:
                    row = dict(status="error", nomination=nomination, exception_type=type(error).__name__,
                               message=str(error), traceback=traceback.format_exc())
                    logger.exception("Cohort call failed: %s / %s", source["source_id"], variant)
                db.execute("INSERT INTO calls VALUES (?, ?)", (i, canonical_json(row)))
                if i % 25 == 0:
                    db.commit()
                    logger.info("%s: %d/%d calls", source["source_id"], i + 1, len(variants))
            db.commit()
        outcomes = Counter()
        temporary = output / "variants.jsonl.gz.tmp"
        with temporary.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as stream:
            for value, in db.execute("SELECT result FROM calls ORDER BY ordinal"):
                outcomes[json.loads(value)["status"]] += 1
                stream.write((value + "\n").encode())
        temporary.replace(output / "variants.jsonl.gz")
    unresolved = outcomes["unresolved_chromosome"]
    if unresolved:
        logger.warning("%s: %d nominations have unresolved chromosomes; RNA support was not measured. "
                       "See per-call contig_resolution diagnostics.", source["source_id"], unresolved)
    errors = outcomes["error"]
    receipt = dict(schema="isovar.cohort_source.v1", identity=identity, identity_key=key,
                   annotation_diagnostics=annotation_diagnostics,
                   matched_normal_scope=normal_scope,
                   status="completed_with_errors" if errors or unresolved else "completed", errors=errors,
                   unresolved_chromosomes=unresolved,
                   nominated_alleles=len(variants), results="variants.jsonl.gz",
                   result_sha256=file_sha256(output / "variants.jsonl.gz"))
    _write_companions(output / "variants.jsonl.gz", output)
    write_json(output / "receipt.json", receipt)
    return receipt



def _write_companions(path, output):
    """Export every retained protein and a small per-call table for inspection."""
    with gzip.open(path, "rt") as stream, (output / "proteins.fasta").open("w") as fasta, \
            (output / "variants.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["variant", "genes", "category", "status", "alt_fragments", "protein_hypotheses",
                         "changed_protein_hypotheses", "ncORF_candidates", "error"])
        for line in stream:
            row = json.loads(line)
            v = row["nomination"]
            name = f"{v['contig']}:{v['start']}:{v['ref']}>{v['alt']}"
            if row["status"] != "assessed":
                writer.writerow([name, "", "", row["status"], "", "", "", "", row["message"]])
                continue
            event = row["protein_hypotheses"]["events"][0]
            proteins = event["protein_hypotheses"]
            changed = [p for p in proteins if protein_comparison_status(p, focal=True) == "amino_acid_divergent"]
            nc = row["nc_orfs"].get("result")
            nc_proteins = nc["protein_candidates"] if nc else []
            writer.writerow([name, ";".join(row["gene_names"]), row["category"], row["status"],
                             event["allele_support"]["alt"]["fragments"], len(proteins), len(changed),
                             len(nc_proteins), ""])
            for protein in proteins:
                sequence = protein["amino_acids"]
                header = f"{protein['hypothesis_id']} {name} focal_change={protein in changed}"
                fasta.write(">" + header + "\n" + "\n".join(sequence[i:i + 80]
                                                               for i in range(0, len(sequence), 80)) + "\n")
            for protein in nc_proteins:
                sequence = protein["protein_sequence"]
                header = (f"{protein['candidate_id']} {name} ncORF={protein['orf_id']} "
                          f"focal_change={protein['focal_amino_acid_change']}")
                fasta.write(">" + header + "\n" + "\n".join(sequence[i:i + 80]
                                                               for i in range(0, len(sequence), 80)) + "\n")

def run_cohort(manifest_path, output_dir, *, resume=False, render_only=False, source_ids=None, render=True):
    """Run source tasks or render identity-verified recorded results.

    Rendering requires receipts matching the manifest's scientific declarations
    and verifies result digests. It does not require original BAMs, rerun analysis
    or relabel results when the manifest's sample, source or settings change.
    """
    manifest_path, output = Path(manifest_path).resolve(), Path(output_dir).resolve()
    manifest = read_json(manifest_path)
    if manifest.get("schema") != SCHEMA:
        raise ValueError("Unsupported cohort manifest schema")
    sources = {row["source_id"]: row for row in manifest["sources"]}
    sample_ids = [sample["sample_id"] for sample in manifest["samples"]]
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError("Duplicate cohort sample ID")
    if any(not key or Path(key).name != key or key in {".", ".."} for key in sample_ids):
        raise ValueError("Sample IDs must be single path components")
    if len(sources) != len(manifest["sources"]):
        raise ValueError("Duplicate cohort source ID")
    for key, source in sources.items():
        if not all(source.get(k) for k in ("sample_id", "library_id", "provenance", "scope")):
            raise ValueError("Cohort sources require sample, library, scope and provenance: " + key)
        if not key or Path(key).name != key or key in {".", ".."}:
            raise ValueError("Source IDs must be single path components")
        if source["sample_id"] not in sample_ids:
            raise ValueError("Source refers to an undeclared sample: " + key)
        if "input_sample_id" in source and (not isinstance(source["input_sample_id"], str) or not source["input_sample_id"].strip()):
            raise ValueError("Source input_sample_id must be a nonempty string")
    for sample in manifest["samples"]:
        for source_id in sample.get("source_ids", []):
            if source_id not in sources or sources[source_id]["sample_id"] != sample["sample_id"]:
                raise ValueError("Sample source_ids must name sources belonging to that sample")
    if source_ids and set(source_ids) - set(sources):
        raise ValueError("Requested source is absent from manifest")
    for family in ("small_variants", "sv_audits", "fusions", "viral", "isoforms", "prior_assessments"):
        tasks = manifest.get(family, [])
        identifiers = [t["id"] for t in tasks] if family != "small_variants" else []
        if len(identifiers) != len(set(identifiers)):
            raise ValueError("Duplicate task ID within " + family)
        for task in tasks:
            source = sources[task["source_id"]]
            if family != "small_variants" and not task.get("sample_id"):
                raise ValueError("Additional tasks require sample_id")
            if task.get("sample_id", source["sample_id"]) != source["sample_id"]:
                raise ValueError("Task and source refer to different samples")
    output.mkdir(parents=True, exist_ok=True)
    files = _inputs(manifest, manifest_path.parent, _selected_file_keys(manifest, source_ids)) if not render_only else {
        key: (manifest_path.parent / value["path"]).resolve() for key, value in manifest["files"].items()}
    outcomes = []
    tasks = manifest.get("small_variants", [])
    task_keys = [(t["source_id"], t.get("pass_id", "default")) for t in tasks]
    if len(set(task_keys)) != len(task_keys):
        raise ValueError("Duplicate source/pass task")
    for source_id, pass_id in task_keys:
        if not pass_id or Path(pass_id).name != pass_id or pass_id in {".", ".."}:
            raise ValueError("Pass IDs must be single path components")
    for task in tasks:
        source = sources[task["source_id"]]
        if source_ids and source["source_id"] not in source_ids:
            continue
        directory = output / "sources" / source["source_id"] / task.get("pass_id", "default")
        if source.get("unavailable_reason"):
            row = dict(status="unavailable", reason=source["unavailable_reason"])
        elif render_only:
            row = _render_small_variant_receipt(task, source, manifest, directory)
        else:
            row = run_small_variant_source(task, source, files, directory, resume=resume)
        outcomes.append(dict(source_id=source["source_id"], sample_id=source["sample_id"],
                             pass_id=task.get("pass_id", "default"), output=str(directory), **row))
        if not render_only:
            write_json(output / ("progress-" + source["source_id"] + ".json"), dict(small_variants=outcomes))
    additional = []
    for family in ("sv_audits", "fusions", "viral", "isoforms", "prior_assessments"):
        for task in manifest.get(family, []):
            if source_ids and task.get("source_id") not in source_ids:
                continue
            identifier = task["id"]
            if not identifier or Path(identifier).name != identifier or identifier in {".", ".."}:
                raise ValueError("Task IDs must be single path components")
            source = sources[task["source_id"]]
            if source.get("unavailable_reason"):
                additional.append(dict(family=family, task=task, status="unavailable",
                                       reason=source["unavailable_reason"], path=None, sha256=None))
                continue
            directory = output / family / identifier
            if not render_only:
                directory.mkdir(parents=True, exist_ok=True)
            target = directory / "result.json.gz"
            if family == "sv_audits":
                target = directory / "report.json"
            expected = _report_identity(task, source, manifest, family)
            receipt_path = directory / "receipt.json"
            receipt = None
            if render_only or resume:
                receipt = _check_additional_receipt(receipt_path, expected, directory, target, resume=resume and not render_only)
            elif target.exists() or receipt_path.exists():
                raise ValueError("Additional result exists; use --resume or a new output directory")
            if not render_only and receipt is None:
                identity = dict(expected, software=reconstruction_identity())
                if family == "sv_audits":
                    run_sv_rna_audit(files[task["manifest"]], directory, resume=resume)
                elif family == "fusions":
                    from .fusion import fusion_from_dict, reconstruct_fusion

                    fusion, references, reads = fusion_from_dict(read_json(files[task["input"]]))
                    if fusion.provenance["sample_id"] != source.get("input_sample_id", source["sample_id"]):
                        raise ValueError("Supplied fusion and cohort source refer to different samples")
                    result = reconstruct_fusion(fusion, references, reads, **task.get("parameters", {}))
                    write_json(target, result)
                elif family == "prior_assessments":
                    write_json(target, read_json(files[task["input"]]))
                elif family == "viral":
                    from .viral_sequences import reconstruct_viral_sequences, write_viral_fasta
                    from .viral_coding import write_viral_proteins

                    source = sources[task["source_id"]]
                    with pysam.AlignmentFile(files[source["bam"]]) as bam:
                        result = reconstruct_viral_sequences(
                            bam, genome_fasta=files[task["genome_fasta"]],
                            sample_id=source["sample_id"], source=source["source_id"],
                            **task.get("parameters", {}))
                    write_json(target, result)
                    write_viral_fasta(result, directory / "sequences.fasta")
                    write_viral_proteins(result, directory / "proteins.fasta")
                else:
                    from .antigen_isoforms import reconstruct_antigen_isoforms
                    from .antigen_reference import IsoformTarget
                    from .fusion import FusionReference

                    targets = [IsoformTarget(**dict(row, reference=FusionReference(**row["reference"])))
                               for row in read_json(files[task["targets"]])]
                    source = sources[task["source_id"]]
                    with pysam.AlignmentFile(files[source["bam"]]) as bam:
                        result = reconstruct_antigen_isoforms(
                            bam, targets=targets, reference_identity=task["reference_identity"],
                            sample_id=source["sample_id"], source=source["source_id"],
                            **task.get("parameters", {}))
                    write_json(target, result)
                receipt = dict(schema="isovar.cohort_additional.v1", identity=identity,
                               artifacts=_additional_artifacts(family, directory, target))
                write_json(receipt_path, receipt)
            additional.append(dict(family=family, task=task, path=str(target),
                                   status="recorded" if receipt else "not_run",
                                   sha256=receipt["artifacts"][target.name] if receipt else None,
                                   receipt=str(receipt_path) if receipt else None,
                                   identity=receipt["identity"] if receipt else None))
    summary = dict(schema="isovar.cohort_results.v1", manifest=str(manifest_path),
                   samples=manifest["samples"], sources=list(sources.values()), small_variants=outcomes, additional=additional)
    summary_name = ("cohort-" + content_identifier("selection", sorted(source_ids)) + ".json"
                    if source_ids else "cohort.json")
    _write_manifest(output / "manifest.json", manifest)
    write_json(output / summary_name, summary)
    if not render:
        return []
    from .cohort_pdf import write_cohort_reports

    return write_cohort_reports(manifest, summary, output)
