"""Manifest-driven, source-separated nominated-SV reconstruction and audit."""

import csv
from contextlib import ExitStack
import gzip
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import pysam

from . import __version__
from ._audit_json import write_json_stream
from .audit_bundle import software_identity
from .fusion import fusion_from_dict, reconstruct_fusion
from .read_collector import ReadCollector
from .rna_library import prepare_rna_library
from .sv_rna import reconstruct_sv_rna, sv_rna_input_from_dict
from .sv_rna_assessment import assess_sv_rna_paths
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
    seen = set()
    with opener(path, "rt") as stream:
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
            if key in seen:
                raise ValueError("Duplicate original FASTQ segment identity")
            seen.add(key)
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


def run_sv_rna_audit(manifest_path, output_dir, *, resume=False, event_ids=None, library_cache_dir=None):
    """Run all declared events/views/sources, preserving unresolved outcomes.

    Relative input paths resolve against the manifest. Checksums are verified
    before any computation. Resume requires identical inputs, parameters,
    producer implementation and a verified output digest. Libraries are never
    pooled, and unavailable sources produce explicit unassessed rows.

    Dense views automatically share one original-template index per source and
    reference context. ``library_cache_dir`` optionally retains verified indexes
    across invocations; by default they have this audit's temporary lifetime.
    Cache location does not change the scientific fingerprint or budgets.
    """
    manifest_path, output_dir = Path(manifest_path).resolve(), Path(output_dir)
    manifest = read_json(manifest_path)
    if manifest.get("schema") != "isovar.sv_rna_audit_manifest.v1":
        raise ValueError("Unsupported SV audit manifest schema")
    allowed = {"schema", "description", "files", "sources", "events", "parameters", "read_collection", "assessment"}
    if set(manifest) - allowed:
        raise ValueError("Unknown audit manifest field")
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
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    with ExitStack() as resources:
        library_directory = (Path(library_cache_dir).expanduser().resolve() if library_cache_dir is not None
                             else Path(resources.enter_context(TemporaryDirectory(prefix='isovar-audit-libraries-'))))
        libraries = {}
        for event in manifest["events"]:
            if event_ids is not None and event["event_id"] not in event_ids:
                continue
            for source_id in event["sources"]:
                source = sources[source_id]
                for view in event.get("views", []) or [dict(view_id="unresolved")]:
                    identity = [event["event_id"], source_id, view["view_id"]]
                    artifact = sha256(json.dumps(identity).encode()).hexdigest()[:24] + ".json.gz"
                    target, receipt = output_dir / artifact, output_dir / (artifact + ".receipt.json")
                    fingerprint = sha256(json.dumps(dict(manifest=manifest, event=event, source=source, view=view,
                                                         implementation=implementation, producer_version=__version__,
                                                         dependencies=dependencies, index_hashes=index_hashes), sort_keys=True).encode()).hexdigest()
                    if resume and target.exists() and receipt.exists():
                        prior = read_json(receipt)
                        if prior["fingerprint"] != fingerprint or prior["sha256"] != file_sha256(target):
                            raise ValueError("Cannot resume changed inputs, implementation or output: " + artifact)
                        rows.append(_summary(read_json(target), artifact))
                        continue
                    if target.exists() or receipt.exists():
                        raise ValueError("Audit output already exists; use verified resume or a new directory")
                    snapshot = dict(schema="isovar.sv_rna_audit_event.v1", event_id=event["event_id"],
                                    source_id=source_id, view_id=view["view_id"], source=source,
                                    event_provenance=event["provenance"], producer_version=__version__,
                                    implementation_sha256=implementation, dependency_versions=dependencies,
                                    index_sha256=index_hashes, input_manifest=manifest,
                                    fingerprint=fingerprint)
                    missing = source.get("unavailable_reason") or event.get("unassessed_reason")
                    if missing:
                        snapshot.update(status="unassessed", limitations=[missing])
                    else:
                        parameters = dict(manifest.get("parameters", {}), **view.get("parameters", {}))
                        models = read_json(files[view["references"]])
                        references = models["models"] if isinstance(models, dict) else models
                        data = dict(view["input"], references=references, event_id=event["event_id"],
                                    sample_id=source["sample_id"], event_provenance=event["provenance"])
                        inputs = sv_rna_input_from_dict(data)
                        raw = tuple(raw_fastq_reads(files[source["raw_fastq"]],
                                                   source.get("fastq_name_format", "sra_dot_mate"),
                                                   source.get("raw_read_group", ""))) if source.get("raw_fastq") else ()
                        with pysam.AlignmentFile(str(files[source["bam"]]), "rb") as bam:
                            if view.get("supplied_fusion"):
                                supplied = fusion_from_dict(read_json(files[view["supplied_fusion"]]))
                                result = reconstruct_fusion(*supplied)
                                if (result["sample_id"] != source["sample_id"] or result["event_id"] != event["event_id"]
                                        or result["reference_name"] != inputs["reference_name"]):
                                    raise ValueError("Supplied fusion identity mismatch")
                                if tuple(inputs["references"]) != supplied[1]:
                                    raise ValueError("Supplied fusion reference model mismatch")
                                result["source"] = source_id
                            else:
                                collector = ReadCollector(**manifest.get("read_collection", {}))
                                if 'library_index' in parameters:
                                    raise ValueError('Use library_cache_dir for audit indexes, not a manifest library_index')
                                if parameters.get('dense_support'):
                                    library_key = source_id, inputs['reference_name']
                                    if library_key not in libraries:
                                        prepared = prepare_rna_library(files[source['bam']], library_directory,
                                            read_collector=collector, reference_name=inputs['reference_name'],
                                            read_filter_identity=parameters.get('read_filter_identity'))
                                        if prepared.identity['source_sha256'] != manifest['files'][source['bam']]['sha256']:
                                            raise ValueError('Original RNA changed after manifest verification: ' + source_id)
                                        libraries[library_key] = resources.enter_context(prepared.open())
                                    parameters['library_index'] = libraries[library_key]
                                result = reconstruct_sv_rna(bam, source=source_id, read_collector=collector,
                                                            **inputs, **parameters)
                            print("%s / %s / %s: assessing %d paths" % (*identity, len(result["paths"])), flush=True)
                            assessment = assess_sv_rna_paths(result, bam, inputs["references"], raw_reads=raw,
                                                            **manifest.get("assessment", {}))
                        hypotheses = read_json(files[view["hypotheses"]]) if view.get("hypotheses") else []
                        reconciliation = reconcile_sv_rna_hypotheses(assessment, hypotheses)
                        assessment["limitations"] = list(dict.fromkeys(
                            assessment["limitations"] + source["provenance"].get("limitations", [])
                            + (models.get("limitations", []) if isinstance(models, dict) else [])))
                        reconciliation["limitations"] = assessment["limitations"]
                        snapshot.update(status=result["status"], reconstruction=result, assessment=assessment,
                                        reconciliation=reconciliation)
                    write_json(target, snapshot)
                    write_json(receipt, dict(fingerprint=fingerprint, sha256=file_sha256(target)))
                    rows.append(_summary(snapshot, artifact))
                    print("%s / %s / %s: %s" % (*identity, snapshot["status"]), flush=True)
    report = dict(schema="isovar.sv_rna_audit_report.v1", manifest_sha256=file_sha256(manifest_path),
                  producer_version=__version__, implementation_sha256=implementation, dependency_versions=dependencies,
                  index_sha256=index_hashes, selected_event_ids=sorted(event_ids) if event_ids is not None else None, runs=rows)
    write_json(output_dir / "report.json", report)
    fields = list(rows[0]) if rows else ["event_id", "source_id", "view_id", "status"]
    with (output_dir / "report.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(dict(row, limitations=";".join(row["limitations"])) for row in rows)
    return report
