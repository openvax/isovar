"""Offline replay of structured SV RNA plans from an :class:`AuditBundle`."""

import gzip
import json
from pathlib import Path

import pysam

from .audit_bundle import AuditBundleError, canonical, identity, publish_json, software_identity
from .read_collector import ReadCollector
from .sv_rna import reconstruct_sv_rna, sv_rna_input_from_dict
__all__ = ["replay_sv_rna"]


def read_json(path):
    path = Path(path)
    with gzip.open(path, "rt") if path.suffix == ".gz" else path.open() as handle:
        return json.load(handle)


def replay_sv_rna(bundle, run_id, output_directory, *, allow_software_change=False, scratch_dir=None):
    """Reconstruct one frozen SV view, compare its result and safely resume.

    Parameters
    ----------
    bundle : AuditBundle
        Verified local files and a ``metadata.sv_rna_runs`` plan dictionary.
    run_id : str
        Exact source/event/orientation run identifier in that dictionary.
    output_directory : path-like
        Checkpoint directory, outside the immutable bundle. Existing matching
        checkpoints are reused only after their content is verified.
    allow_software_change : bool
        Permit a different engine/dependency identity. The frozen result must
        still match; a changed implementation is never silently accepted.
    scratch_dir : path-like, optional
        Disposable scratch storage for the existing dense reconstruction engine.

    Returns
    -------
    dict
        Checkpoint containing request, full reconstruction, result checksum and
        comparison. A mismatch raises AuditBundleError and retains its checkpoint.
    """
    bundle.verify()
    plans = bundle.manifest["metadata"].get("sv_rna_runs", {})
    if not isinstance(plans, dict):
        raise AuditBundleError("SV RNA replay runs must be a mapping of run IDs to plans")
    if not isinstance(run_id, str) or run_id not in plans:
        raise AuditBundleError("Unknown SV RNA replay run: " + str(run_id))
    plan = plans[run_id]
    required = {"bam", "index", "input", "source", "parameters"}
    allowed = required | {"expected", "read_collection", "input_limitations", "upstream_acquisition", "original_request"}
    if not isinstance(plan, dict) or required - set(plan) or set(plan) - allowed:
        raise AuditBundleError("Missing or unexpected SV RNA replay plan fields")
    if not isinstance(plan["source"], str) or not plan["source"]:
        raise AuditBundleError("Replay source must be an explicit nonempty identity")
    if not isinstance(plan["parameters"], dict) or not isinstance(plan.get("read_collection", {}), dict):
        raise AuditBundleError("Replay parameters and read collection settings must be objects")
    if not isinstance(plan.get("input_limitations", []), list) or any(
            not isinstance(flag, str) for flag in plan.get("input_limitations", [])):
        raise AuditBundleError("Replay input limitations must be a list of strings")
    for key in ("bam", "index", "input", "expected"):
        if key in plan:
            bundle.path(plan[key])
    software = software_identity()
    if not allow_software_change and software != bundle.manifest["software"]:
        raise AuditBundleError("Replay software differs from the bundle; use allow_software_change explicitly")
    output = Path(output_directory).resolve()
    if output.is_relative_to(bundle.directory):
        raise AuditBundleError("Replay output must be outside the immutable bundle")
    request = dict(bundle_id=bundle.bundle_id, run_id=run_id, plan_sha256=identity(plan), software=software)
    path = output / (identity(run_id)[:32] + ".json")
    expected = read_json(bundle.path(plan["expected"])) if plan.get("expected") else None
    if path.exists():
        checkpoint = read_json(path)
        if checkpoint.get("request") != request:
            raise AuditBundleError("Replay checkpoint request drift: " + str(path))
        if "result" not in checkpoint or checkpoint.get("result_sha256") != identity(checkpoint["result"]):
            raise AuditBundleError("Replay checkpoint result checksum mismatch: " + str(path))
    else:
        inputs = sv_rna_input_from_dict(read_json(bundle.path(plan["input"])))
        collector = ReadCollector(**plan.get("read_collection", {}))
        with pysam.AlignmentFile(bundle.path(plan["bam"]), index_filename=bundle.path(plan["index"])) as bam:
            result = reconstruct_sv_rna(bam, **inputs, source=plan["source"],
                                        read_collector=collector, scratch_dir=scratch_dir, **plan["parameters"])
        result["limitations"] = sorted(set(result["limitations"]) | set(plan.get("input_limitations", [])))
        if "upstream_acquisition" in plan:
            result["upstream_acquisition"] = plan["upstream_acquisition"]
        # The public checkpoint is JSON data on both first run and resume.
        # Tuples in native engine output serialize as lists in frozen evidence.
        result = json.loads(canonical(result))
        comparison = "not_requested" if expected is None else "matched" if result == expected else "mismatch"
        checkpoint = dict(request=request, result=result, result_sha256=identity(result), comparison=comparison)
        publish_json(path, checkpoint)
    comparison = "not_requested" if expected is None else "matched" if checkpoint["result"] == expected else "mismatch"
    if checkpoint.get("comparison") != comparison:
        raise AuditBundleError("Replay checkpoint comparison was changed: " + str(path))
    if comparison == "mismatch":
        raise AuditBundleError("Replay differs from the frozen result; checkpoint retained at " + str(path))
    return checkpoint
