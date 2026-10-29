"""Manifest-driven, source-separated nominated-SV reconstruction and audit."""

import csv
from collections import defaultdict
from collections.abc import Iterable
from contextlib import ExitStack, closing
from dataclasses import dataclass
import gzip
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path
import sqlite3
from tempfile import TemporaryDirectory
from typing import Any

import pysam

from . import __version__
from ._audit_json import write_json_stream
from .audit_bundle import software_identity
from .fusion import fusion_from_dict, reconstruct_fusion
from .read_collector import ReadCollector
from .rna_library import RnaLibraryIndex, prepare_rna_library
from .sv_path_routing import PathEvidenceQuery, route_path_evidence
from .sv_rna import reconstruct_sv_rna, sv_rna_input_from_dict
from .sv_rna_assessment import PathEvidence, assess_sv_rna_paths
from .sv_rna_reconciliation import reconcile_sv_rna_hypotheses


def file_sha256(path):
    """Hash file bytes without materializing large original-read files."""
    digest = sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as stream:
        return json.load(stream)


def write_json(path, value):
    """Write deterministic bytes atomically, including compressed snapshots."""
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    try:
        with temporary.open("wb") as stream:
            write_json_stream(stream, value, compressed=str(path).endswith(".gz"))
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def raw_fastq_reads(path, name_format="sra_dot_mate", read_group=""):
    """Decode original read identities using an explicitly declared convention."""
    if name_format not in {"sra_dot_mate", "slash_mate", "unpaired"}:
        raise ValueError("Unknown original FASTQ name format")
    opener = gzip.open if str(path).endswith(".gz") else open
    with TemporaryDirectory(prefix='isovar-fastq-identities-') as scratch, \
            closing(sqlite3.connect(Path(scratch) / 'identities.sqlite')) as identities, opener(path, "rt") as stream:
        identities.execute('CREATE TABLE segments (name TEXT, mate INTEGER, PRIMARY KEY (name, mate)) WITHOUT ROWID')
        while name := stream.readline().strip():
            sequence, plus, quality = (stream.readline().strip() for _ in range(3))
            if not name.startswith("@") or not plus.startswith("+") or len(sequence) != len(quality):
                raise ValueError("Malformed original FASTQ")
            name = name[1:].split()[0]
            mate = 0
            if name_format != "unpaired":
                name, end = name.rsplit("." if name_format == "sra_dot_mate" else "/", 1)
                if end not in {"1", "2"}:
                    raise ValueError("Unexpected original FASTQ mate")
                mate = int(end)
            key = (read_group, name, mate)
            try:
                identities.execute('INSERT INTO segments VALUES (?, ?)', (name, mate))
            except sqlite3.IntegrityError:
                raise ValueError("Duplicate original FASTQ segment identity")
            yield *key, sequence.upper()


def _files(manifest, root):
    files = {}
    for identifier, descriptor in manifest["files"].items():
        if set(descriptor) - {"path", "sha256", "source_url", "description"}:
            raise ValueError("Unknown audit file descriptor field")
        path = (root / descriptor["path"]).resolve()
        expected = descriptor["sha256"]
        if len(expected) != 64 or file_sha256(path) != expected:
            raise ValueError("Audit input checksum mismatch: " + identifier)
        files[identifier] = path
    return files


def _summary(snapshot, artifact):
    assessment = snapshot.get("assessment", {})
    paths = assessment.get("paths", [])
    return dict(event_id=snapshot["event_id"], source_id=snapshot["source_id"],
                view_id=snapshot["view_id"], status=snapshot["status"],
                retained_paths=len(paths),
                paths_with_direct_junction=sum(any(j["direct_support"].get("fragments", 0)
                                                  and j["relation"] == "breakpoint_junction"
                                                  for j in p["junctions"]) for p in paths),
                paths_with_event_linked_protein=sum(any(c["event_linked"] for c in p["proteins"]) for p in paths),
                limitations=snapshot.get("limitations", assessment.get("limitations", [])),
                artifact=artifact)


@dataclass(frozen=True)
class _AuditContext:
    manifest: dict[str, Any]
    files: dict[str, Path]
    sources: dict[str, dict[str, Any]]
    implementation: dict[str, str]
    dependencies: dict[str, str]
    index_hashes: dict[str, dict[str, str]]


@dataclass(frozen=True)
class _AuditJob:
    event: dict[str, Any]
    source: dict[str, Any]
    view: dict[str, Any]
    identity: tuple[str, str, str]
    artifact: str
    target: Path
    receipt: Path
    fingerprint: str


def _load_audit_context(manifest_path: Path, event_ids: set[str] | None) -> _AuditContext:
    """Validate manifest identities and pin every source before planning work."""
    manifest = read_json(manifest_path)
    if manifest.get("schema") != "isovar.sv_rna_audit_manifest.v1":
        raise ValueError("Unsupported SV audit manifest schema")
    allowed = {"schema", "description", "files", "sources", "events", "parameters", "read_collection", "assessment"}
    if set(manifest) - allowed:
        raise ValueError("Unknown audit manifest field")
    thresholds = manifest.get('assessment', {})
    if set(thresholds) - {'min_path_match', 'min_exon_overlap', 'junction_flank', 'min_overlap', 'max_chain_reads'}:
        raise ValueError('Unknown path assessment parameter')
    if any(type(value) is not int or value < 1 for value in thresholds.values()):
        raise ValueError('Assessment thresholds must be positive integers')
    files = _files(manifest, manifest_path.parent)
    sources = {source["source_id"]: source for source in manifest["sources"]}
    if len(sources) != len(manifest["sources"]):
        raise ValueError("Duplicate source identity")
    ids = [event["event_id"] for event in manifest["events"]]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate event identity")
    if event_ids is not None and set(event_ids) - set(ids):
        raise ValueError("Requested event not in manifest")
    for source in sources.values():
        if not all(source.get(k) for k in ("sample_id", "library_id", "culture", "role", "provenance")):
            raise ValueError("Sources require sample, library, culture, role and provenance")
        if bool(source.get("unavailable_reason")) == bool(source.get("bam")):
            raise ValueError("A source needs either a BAM or an explicit unavailable reason")
    # All package modules and numerical/sequence dependencies participate in
    # cache identity. A helper-only change must not reuse stale scientific work.
    software = software_identity()
    implementation = software["engine_sha256"]
    dependencies = dict(software["dependencies"], **{name: version(name) for name in ("varcode", "pyensembl", "orfdata")})
    index_hashes = {}
    for source in sources.values():
        if source.get("bam"):
            bam_path = files[source["bam"]]
            indexes = [bam_path.with_suffix(".bai"), Path(str(bam_path) + ".bai"), Path(str(bam_path) + ".csi")]
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                if not bam.has_index():
                    raise ValueError("Audit BAM requires an index: " + source["source_id"])
            index_hashes[source["source_id"]] = {p.name: file_sha256(p) for p in indexes if p.exists()}
        if source.get("raw_fastq") not in (None, *files):
            raise ValueError("Unknown original FASTQ file")
    for event in manifest["events"]:
        if not event.get("sources") or len(event["sources"]) != len(set(event["sources"])):
            raise ValueError("Events require distinct source IDs")
        if set(event["sources"]) - set(sources) or not event.get("provenance"):
            raise ValueError("Event source or provenance missing")
        views = event.get("views", [])
        if not views and not event.get("unassessed_reason"):
            raise ValueError("Event requires views or an explicit unassessed reason")
        if len({v["view_id"] for v in views}) != len(views):
            raise ValueError("Duplicate view identity")
        for view in views:
            for key in ("references", "hypotheses", "supplied_fusion"):
                if view.get(key) not in (None, *files):
                    raise ValueError("Unknown view input file: " + key)
            if "references" not in view or "input" not in view:
                raise ValueError("View requires reference models and an oriented input")
    return _AuditContext(manifest, files, sources, implementation, dependencies, index_hashes)


def _plan_jobs(context: _AuditContext, output: Path, event_ids: set[str] | None) -> list[_AuditJob]:
    """Preserve inventory order while assigning stable per-view artifacts."""
    jobs = []
    for event in context.manifest['events']:
        if event_ids is not None and event['event_id'] not in event_ids:
            continue
        for source_id in event['sources']:
            source = context.sources[source_id]
            for view in event.get('views', []) or [dict(view_id='unresolved')]:
                identity = event['event_id'], source_id, view['view_id']
                artifact = sha256(json.dumps(identity).encode()).hexdigest()[:24] + '.json.gz'
                fingerprint = sha256(json.dumps(dict(manifest=context.manifest, event=event, source=source, view=view,
                    implementation=context.implementation, producer_version=__version__,
                    dependencies=context.dependencies, index_hashes=context.index_hashes), sort_keys=True).encode()).hexdigest()
                jobs.append(_AuditJob(event, source, view, identity, artifact, output / artifact,
                                      output / (artifact + '.receipt.json'), fingerprint))
    return jobs


def _read_verified_artifact(target: Path, receipt: Path, fingerprint: str, resume: bool) -> dict[str, Any] | None:
    """A partial or changed checkpoint is an error, never evidence of no support."""
    if resume and target.exists() and receipt.exists():
        prior = read_json(receipt)
        if prior['fingerprint'] != fingerprint or prior['sha256'] != file_sha256(target):
            raise ValueError('Cannot resume changed inputs, implementation or output: ' + target.name)
        return read_json(target)
    if target.exists() or receipt.exists():
        raise ValueError('Audit output already exists; use verified resume or a new directory')
    return None


def _snapshot(job: _AuditJob, context: _AuditContext) -> dict[str, Any]:
    return dict(schema='isovar.sv_rna_audit_event.v1', event_id=job.identity[0], source_id=job.identity[1],
                view_id=job.identity[2], source=job.source, event_provenance=job.event['provenance'],
                producer_version=__version__, implementation_sha256=context.implementation,
                dependency_versions=context.dependencies, index_sha256=context.index_hashes,
                input_manifest=context.manifest, fingerprint=job.fingerprint)


def _publish_snapshot(job: _AuditJob, snapshot: dict[str, Any]) -> dict[str, Any]:
    write_json(job.target, snapshot)
    write_json(job.receipt, dict(fingerprint=job.fingerprint, sha256=file_sha256(job.target)))
    print('%s / %s / %s: %s' % (*job.identity, snapshot['status']), flush=True)
    return _summary(snapshot, job.artifact)


def _view_inputs(job: _AuditJob, context: _AuditContext) -> tuple[dict[str, Any], list[str]]:
    models = read_json(context.files[job.view['references']])
    references = models['models'] if isinstance(models, dict) else models
    data = dict(job.view['input'], references=references, event_id=job.identity[0],
                sample_id=job.source['sample_id'], event_provenance=job.event['provenance'])
    return sv_rna_input_from_dict(data), models.get('limitations', []) if isinstance(models, dict) else []


def _stage_paths(job: _AuditJob) -> tuple[Path, Path]:
    directory = job.target.parent / '.reconstructions'
    return directory / job.artifact, directory / (job.artifact + '.receipt.json')


def _reconstruct_job(job: _AuditJob, context: _AuditContext, library: RnaLibraryIndex,
                     collector: ReadCollector, resume: bool) -> PathEvidenceQuery:
    """Checkpoint reconstruction before the library-wide assessment routing pass."""
    inputs, _ = _view_inputs(job, context)
    target, receipt = _stage_paths(job)
    prior = _read_verified_artifact(target, receipt, job.fingerprint, resume)
    if prior is not None:
        result = prior['reconstruction']
    else:
        parameters = dict(context.manifest.get('parameters', {}), **job.view.get('parameters', {}))
        if 'library_index' in parameters:
            raise ValueError('Use library_cache_dir for audit indexes, not a manifest library_index')
        if job.view.get('supplied_fusion'):
            supplied = fusion_from_dict(read_json(context.files[job.view['supplied_fusion']]))
            result = reconstruct_fusion(*supplied)
            if (result['sample_id'] != job.source['sample_id'] or result['event_id'] != job.identity[0]
                    or result['reference_name'] != inputs['reference_name']):
                raise ValueError('Supplied fusion identity mismatch')
            if tuple(inputs['references']) != supplied[1]:
                raise ValueError('Supplied fusion reference model mismatch')
            result['source'] = job.identity[1]
        else:
            if parameters.get('dense_support'):
                parameters['library_index'] = library
            with pysam.AlignmentFile(context.files[job.source['bam']], 'rb') as bam:
                result = reconstruct_sv_rna(bam, source=job.identity[1], read_collector=collector, **inputs, **parameters)
        target.parent.mkdir(parents=True, exist_ok=True)
        write_json(target, dict(schema='isovar.sv_rna_reconstruction_stage.v1', fingerprint=job.fingerprint,
                               reconstruction=result))
        write_json(receipt, dict(fingerprint=job.fingerprint, sha256=file_sha256(target)))
    print('%s / %s / %s: reconstructed %d paths' % (*job.identity, len(result['paths'])), flush=True)
    parameters = context.manifest.get('assessment', {})
    return PathEvidenceQuery(job.artifact, tuple(path['sequence'] for path in result['paths']),
                             tuple(inputs['references']), parameters.get('min_path_match', 20),
                             parameters.get('min_exon_overlap', 15))


def _assess_job(job: _AuditJob, context: _AuditContext, evidence: PathEvidence) -> dict[str, Any]:
    """Assess and reconcile one staged view without reopening its whole library."""
    inputs, model_limitations = _view_inputs(job, context)
    target, receipt = _stage_paths(job)
    staged = _read_verified_artifact(target, receipt, job.fingerprint, True)
    if staged is None:
        raise ValueError('Missing staged reconstruction: ' + job.artifact)
    result = staged['reconstruction']
    assessment = assess_sv_rna_paths(result, None, inputs['references'], evidence=evidence,
                                    **context.manifest.get('assessment', {}))
    hypotheses = read_json(context.files[job.view['hypotheses']]) if job.view.get('hypotheses') else []
    reconciliation = reconcile_sv_rna_hypotheses(assessment, hypotheses)
    assessment['limitations'] = list(dict.fromkeys(assessment['limitations']
        + job.source['provenance'].get('limitations', []) + model_limitations))
    reconciliation['limitations'] = assessment['limitations']
    snapshot = dict(_snapshot(job, context), status=result['status'], reconstruction=result,
                    assessment=assessment, reconciliation=reconciliation)
    row = _publish_snapshot(job, snapshot)
    # The final receipt is durable before removing this job's verified stage.
    target.unlink()
    receipt.unlink()
    return row


def run_sv_rna_audit(manifest_path: str | Path, output_dir: str | Path, *, resume: bool = False,
                      event_ids: Iterable[str] | None = None, library_cache_dir: str | Path | None = None
                      ) -> dict[str, Any]:
    """Reconstruct, route and assess an inventory through one source-specific flow.

    Inputs and completed/staged artifacts are verified before reuse. Each source
    and assembly context is handled serially: reconstruct all pending views,
    route original BAM/FASTQ sequences once, then assess each view with its full
    template alternatives and reference-exon evidence. Unavailable sources and
    unresolved events remain explicit. Libraries are never pooled.

    ``library_cache_dir`` optionally persists verified original-template indexes.
    Its location does not affect scientific fingerprints or discovery budgets.
    Interrupted reconstructions remain in ``.reconstructions`` for verified
    resume; completed final snapshots retain the existing audit schema.
    """
    manifest_path, output_dir = Path(manifest_path).resolve(), Path(output_dir)
    selected = None if event_ids is None else set(event_ids)
    context = _load_audit_context(manifest_path, selected)
    output_dir.mkdir(parents=True, exist_ok=True)
    jobs = _plan_jobs(context, output_dir, selected)
    rows: dict[str, dict[str, Any]] = {}
    pending: dict[tuple[str, str], list[_AuditJob]] = defaultdict(list)
    # Validate every final checkpoint before starting any RNA work.
    completed = {}
    for job in jobs:
        prior = _read_verified_artifact(job.target, job.receipt, job.fingerprint, resume)
        if prior is not None:
            completed[job.artifact] = _summary(prior, job.artifact)
    for job in jobs:
        if job.artifact in completed:
            rows[job.artifact] = completed[job.artifact]
            continue
        missing = job.source.get('unavailable_reason') or job.event.get('unassessed_reason')
        if missing:
            rows[job.artifact] = _publish_snapshot(job, dict(_snapshot(job, context), status='unassessed', limitations=[missing]))
        else:
            pending[job.identity[1], job.view['input']['reference_name']].append(job)
    del completed
    collector = ReadCollector(**context.manifest.get('read_collection', {}))
    for (source_id, reference_name), group in pending.items():
        with ExitStack() as resources:
            cache = (Path(library_cache_dir).expanduser().resolve() if library_cache_dir is not None
                     else Path(resources.enter_context(TemporaryDirectory(prefix='isovar-audit-libraries-'))))
            source = context.sources[source_id]
            prepared = prepare_rna_library(context.files[source['bam']], cache,
                                           read_collector=collector, reference_name=reference_name)
            if prepared.identity['source_sha256'] != context.manifest['files'][source['bam']]['sha256']:
                raise ValueError('Original RNA changed after manifest verification: ' + source_id)
            with prepared.open() as library, TemporaryDirectory(prefix='isovar-path-routing-') as routing:
                queries = [_reconstruct_job(job, context, library, collector, resume) for job in group]
                raw = (raw_fastq_reads(context.files[source['raw_fastq']], source.get('fastq_name_format', 'sra_dot_mate'),
                                      source.get('raw_read_group', '')) if source.get('raw_fastq') else ())
                with route_path_evidence(library, queries, Path(routing) / 'routes.sqlite', raw_reads=raw) as routed:
                    for job in group:
                        rows[job.artifact] = _assess_job(job, context, routed.evidence(job.artifact))
    ordered = [rows[job.artifact] for job in jobs]
    report = dict(schema='isovar.sv_rna_audit_report.v1', manifest_sha256=file_sha256(manifest_path),
                  producer_version=__version__, implementation_sha256=context.implementation,
                  dependency_versions=context.dependencies, index_sha256=context.index_hashes,
                  selected_event_ids=sorted(selected) if selected is not None else None, runs=ordered)
    write_json(output_dir / 'report.json', report)
    fields = list(ordered[0]) if ordered else ['event_id', 'source_id', 'view_id', 'status']
    with (output_dir / 'report.tsv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(dict(row, limitations=';'.join(row['limitations'])) for row in ordered)
    return report
