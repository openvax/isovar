"""Import existing Sid audit evidence into the public portable bundle format."""

from copy import deepcopy
from pathlib import Path
import tempfile

from .audit_bundle import AuditBundle, AuditBundleError, digest, identity, publish_json
from .audit_replay import read_json
__all__ = ["bundle_sid_audit"]


def _pins(value):
    if isinstance(value, dict):
        if "path" in value and "sha256" in value:
            yield value
        for item in value.values():
            yield from _pins(item)
    elif isinstance(value, list):
        for item in value:
            yield from _pins(item)


def _replay_parameters(result):
    parameters = deepcopy(result["parameters"])
    if parameters.pop("genetic_code") != 1:
        raise AuditBundleError("Unsupported genetic code in archived reconstruction")
    collector = parameters.pop("read_collection")
    if collector.pop("custom_read_filter") or collector.pop("read_end_profile_sha256") is not None:
        raise AuditBundleError("Custom read filters/profiles require an explicit portable replay plan")
    parameters.update({"inclusion_" + key: value for key, value in parameters.pop("inclusion").items()})
    parameters.setdefault("dense_support", "discovery" in result)
    return parameters, collector


def bundle_sid_audit(directory, destination, *, extra_files=None):
    """Freeze an existing Sid audit and generate installed offline SV replay plans.

    Parameters
    ----------
    directory : path-like
        Existing audit with inventory, frozen transcript models, acquisition
        receipts, original regional BAMs/indexes and completed checkpoints.
    destination : path-like
        New bundle directory. Original receipts and results are copied verbatim.
    extra_files : mapping of str to path-like, optional
        Additional explicitly selected annotations, provenance or analysis files.

    Returns
    -------
    AuditBundle
        Portable evidence, SV replay plans, and all preserved outcomes. Missing
        input pins fail export rather than fetching or discarding evidence.
    """
    directory = Path(directory).resolve()
    if Path(destination).resolve().is_relative_to(directory):
        raise AuditBundleError("Create the bundle outside the original audit directory")
    files, aliases, documents, checksums, verified = {}, {}, {}, {}, {}

    def add(path, checksum=None, base=None):
        path = Path(path)
        absolute = path if path.is_absolute() else (directory if base is None else base) / path
        absolute = absolute.resolve()
        if absolute in verified:
            name, observed = verified[absolute]
        else:
            if not absolute.is_file():
                raise AuditBundleError("Missing archived input: " + str(absolute))
            observed = digest(absolute)
            name = ("audit/" + absolute.relative_to(directory).as_posix() if absolute.is_relative_to(directory)
                    else "assets/" + observed + "/" + absolute.name)
            verified[absolute] = name, observed
        if checksum is not None and observed != checksum:
            raise AuditBundleError("Archived input checksum mismatch: " + str(absolute))
        files[name] = absolute
        checksums[name] = observed
        aliases[str(path)] = name
        aliases[str(absolute)] = name
        return name

    root_names = ("inventory.json.gz", "inventory-pin.json", "references.json.gz", "references-pin.json", "batches.json")
    for name in root_names:
        add(directory / name)
    inventory, reference = read_json(directory / root_names[0]), read_json(directory / root_names[2])
    for name, pin_name in ((root_names[0], root_names[1]), (root_names[2], root_names[3])):
        if digest(directory / name) != read_json(directory / pin_name)["sha256"]:
            raise AuditBundleError("Audit root pin mismatch: " + name)
    if reference["inventory_sha256"] != digest(directory / root_names[0]):
        raise AuditBundleError("References belong to another inventory")
    paths = {p for p in directory.glob("*.json*") if p.is_file()}
    for tree in ("sources", "small-variants", "results", "reconstructions", "report", "search"):
        paths.update(p for p in (directory / tree).rglob("*.json*") if p.is_file())
    # Reports/searches declare their generated assets relative to their own
    # output root. Preserve non-JSON tables/FASTA too, without copying cache spools.
    for tree in ("report", "search"):
        for path in sorted((directory / tree).rglob("*")):
            if path.is_file() and "__pycache__" not in path.parts:
                add(path)
    for path in sorted(paths):
        add(path)
        data = read_json(path)
        documents[path] = data
        base = next((directory / tree for tree in ("report", "search")
                     if path.is_relative_to(directory / tree)), directory)
        for pin in _pins(data):
            add(pin["path"], pin["sha256"], base=base)
    acquired = {identity(row): row for path, row in documents.items()
                if path.is_relative_to(directory / "sources") and "batch_id" in row}
    plans, outcomes = {}, []
    with tempfile.TemporaryDirectory() as temporary:
        for path, row in documents.items():
            if not path.is_relative_to(directory / "results"):
                continue
            request = row["request"]
            sid, gid, orientation = (request[k] for k in ("source_id", "geometry_id", "orientation"))
            run_id = "/".join((sid, gid, orientation))
            outcomes.append(dict(run_id=run_id, status=row["status"],
                                 acquisition_status=row["acquisition_status"],
                                 input_limitations=row.get("input_limitations", []),
                                 checkpoint=add(path)))
            if "reconstruction" not in row:
                continue  # Preserve an unassessable/timeout outcome, never invent a replay.
            receipt = acquired.get(request["acquisition_sha256"])
            if receipt is None or receipt["source_id"] != sid or gid not in receipt["request"]["batch"]["geometries"]:
                raise AuditBundleError("Reconstruction has no matching acquisition: " + run_id)
            if request["inventory_sha256"] != digest(directory / root_names[0]) or request["references_sha256"] != digest(directory / root_names[2]):
                raise AuditBundleError("Reconstruction root pins changed: " + run_id)
            if receipt["request"]["source_identity"] != identity(inventory["sources"][sid]):
                raise AuditBundleError("Acquisition source identity changed: " + run_id)
            expected_path = directory / row["reconstruction"]["path"]
            expected = read_json(expected_path)
            if expected["event_id"] != gid or expected["event_provenance"]["orientation"] != orientation:
                raise AuditBundleError("Reconstruction event/orientation changed: " + run_id)
            models = []
            for description in expected["reference_models"]:
                model = deepcopy(reference["models"][description["transcript_id"]])
                model["contig"] = description["contig"]
                models.append(model)
            inputs = {key: deepcopy(expected[key]) for key in
                      ("event_id", "reference_name", "sample_id", "donor", "acceptor", "event_provenance")}
            inputs["references"] = models
            input_path = Path(temporary) / (identity(run_id) + ".json")
            publish_json(input_path, inputs)
            input_name = "replay-inputs/" + input_path.name
            files[input_name] = input_path
            parameters, collector = _replay_parameters(expected)
            plans[run_id] = dict(input=input_name, bam=aliases[receipt["bam"]["path"]],
                                 index=aliases[receipt["index"]["path"]], expected=add(expected_path),
                                 source=expected["source"], parameters=parameters, read_collection=collector,
                                 input_limitations=row.get("input_limitations", []),
                                 upstream_acquisition=expected["upstream_acquisition"], original_request=request)
        for name, path in (extra_files or {}).items():
            if name in files:
                raise AuditBundleError("Additional file name already exists: " + name)
            files[name] = Path(path)
        metadata = dict(kind="sid-sv-audit", original_directory=str(directory), original_paths=aliases,
                        outcomes=outcomes, sv_rna_runs=plans,
                        scope="preserved_existing_audit_inputs_and_outcomes; no new acquisition")
        return AuditBundle.create(destination, files, metadata, expected_sha256=checksums)
