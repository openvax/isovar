"""Offline reconstruction with a checkpoint for every geometry/product/view."""

from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter
from contextlib import contextmanager
from copy import deepcopy
from hashlib import sha256
import inspect
from pathlib import Path
import signal

import pysam

import isovar
from isovar import export_sv_rna_orfs
from isovar.sv_rna import reconstruct_sv_rna, sv_rna_input_from_dict
from .acquisition import USABLE_STATUSES, acquisition_tree, input_limitations, selected_sources
from .inventory import canonical, digest, identity, load_inventory, oriented_breakpoints, read_json, write_json
from . import references


PARAMETERS = dict(assemble=True, breakpoint_window=2000, min_orf_amino_acids=8,
                  max_records=50000, max_queries=10000, max_paths=200,
                  max_extension_segments=500, max_orf_candidates=1000,
                  include_regional_candidates=False, dense_support=True)

EVENT_RELATIONS = {"breakpoint_junction", "event_compatible_junction",
                   "splice_ambiguous_event_junction", "breakpoint_clip_partner_unplaced"}


def implementation_identity():
    """Pin actual engine files, including uncommitted development changes."""
    root = Path(inspect.getfile(isovar)).parent
    files = {"isovar/" + str(path.relative_to(root)): digest(path) for path in sorted(root.rglob("*.py"))}
    audit = Path(__file__).parent
    files.update({"audit/" + path.name: digest(path) for path in sorted(audit.glob("*.py"))})
    return identity(files)


def reconstruct_one(bam, group_id, group, source, reference, orientation, parameters):
    models = [deepcopy(reference["models"][tid])
              for tid in reference["assignments"][group_id]["transcripts"]]
    donor, acceptor = oriented_breakpoints(group["breakends"], reverse=orientation == "reverse")
    names = set(bam.references)

    def resolve(contig):
        matches = {contig, contig.removeprefix("chr"), "chr" + contig.removeprefix("chr")} & names
        if len(matches) != 1:
            raise ValueError("Missing or ambiguous BAM contig: " + contig)
        return next(iter(matches))

    for row in [donor, acceptor, *models]:
        row["contig"] = resolve(row["contig"])
    inputs = sv_rna_input_from_dict(dict(
        event_id=group_id, reference_name="GRCh38", sample_id="product:" + source["id"],
        donor=donor, acceptor=acceptor, references=models,
        event_provenance=dict(nominations=group["nominations"], orientation=orientation,
                              source_claims=source["claims"], support_scope="within_product_only")))
    return reconstruct_sv_rna(bam, source=source["url"], **inputs, **parameters)


def summarize(result):
    """Full witness counts remain separate from local junction and assembly support."""
    # Filter occurrences before unioning witnesses. An identical sequence at a
    # regional junction must not supply witnesses for the nominated adjacency.
    event_paths = []
    for path in result["paths"]:
        orfs = path["exploratory_orfs"]
        candidates = [c for c in orfs["candidates"] if any(
            path["junctions"][j]["relation"] in EVENT_RELATIONS for j in c["crossed_junctions"])]
        event_paths.append(dict(path, exploratory_orfs=dict(orfs, candidates=candidates)))
    export = export_sv_rna_orfs(dict(result, paths=event_paths))
    for candidate in export["candidates"]:
        candidate["event_relations"] = sorted({j["relation"] for o in candidate["occurrences"]
                                                for j in o["junctions"] if j["relation"] in EVENT_RELATIONS})
    return dict(status=result["status"], limitations=result["limitations"],
                records=result["acquisition"]["records"], paths=len(result["paths"]),
                event_paths=sum(any(j["relation"] in EVENT_RELATIONS for j in p["junctions"])
                                for p in result["paths"]),
                candidates=export["candidates"], reference_models=result["reference_models"],
                excluded_records=result["excluded_records"],
                discovery=result.get("discovery"),
                support_acquisition=result.get("support_acquisition"),
                orf_search_truncated=any(p["exploratory_orfs"]["candidate_limit_reached"] for p in result["paths"]))


def result_path(directory, source_id, group_id, orientation):
    return Path(directory) / "results" / source_id / (group_id + "." + orientation + ".json.gz")


def _expired(signum, frame):
    raise TimeoutError("SV reconstruction exceeded the recorded per-view wall-time limit")


@contextmanager
def reconstruction_deadline(seconds):
    """Bound computation without interrupting publication of a checkpoint."""
    previous = signal.signal(signal.SIGALRM, _expired)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)


def run_batch(directory, source_id, batch_id, parameters, seconds, engine_id):
    """One process per acquisition batch; all reconstruction reads are local."""
    leaves, receipts = acquisition_tree(directory, source_id, batch_id)
    source_identity = identity(load_inventory(directory)["sources"][source_id])
    if any(row["request"]["source_identity"] != source_identity for row in receipts.values()):
        raise ValueError("Acquisition belongs to another source")
    counts = Counter()
    for bid in leaves:
        result = _run_leaf(directory, source_id, bid, parameters, seconds, engine_id)
        counts.update(result["outcomes"])
    return dict(source=source_id, batch=batch_id, outcomes=dict(counts), acquisition_leaves=len(leaves))


def _run_leaf(directory, source_id, batch_id, parameters, seconds, engine_id):
    directory = Path(directory)
    manifest, reference = load_inventory(directory), references.load(directory)
    source = manifest["sources"][source_id]
    acquired = read_json(directory / "sources" / source_id / (batch_id + ".json"))
    batch = acquired["request"]["batch"]
    for key in ("bam", "index"):
        if key in acquired and digest(acquired[key]["path"]) != acquired[key]["sha256"]:
            raise ValueError("Acquired %s checksum mismatch" % key)
    request_base = dict(inventory_sha256=digest(directory / "inventory.json.gz"),
                        references_sha256=digest(directory / "references.json.gz"),
                        acquisition_sha256=identity(acquired), parameters=parameters,
                        isovar_version=isovar.__version__, engine_id=engine_id, max_seconds=seconds)
    counts = {}
    for group_id in batch["geometries"]:
        group = manifest["geometries"][group_id]
        needed_contigs = {end["contig"] for end in group["breakends"]}
        needed_contigs.update(reference["models"][tid]["contig"]
                              for tid in reference["assignments"][group_id]["transcripts"])
        missing_contigs = needed_contigs & set(acquired.get("unavailable_contigs", []))
        for orientation in ("forward", "reverse"):
            request = dict(request_base, source_id=source_id, geometry_id=group_id, orientation=orientation)
            path = result_path(directory, source_id, group_id, orientation)
            if path.exists():
                saved = read_json(path)
                if saved["request"] != request:
                    raise ValueError("Reconstruction request drift: " + str(path))
                counts["reused"] = counts.get("reused", 0) + 1
                continue
            result = dict(request=request, acquisition_status=acquired["status"],
                          input_limitations=input_limitations(acquired))
            if acquired["status"] not in USABLE_STATUSES:
                result.update(status="not_assessable", reason=acquired["status"])
            elif missing_contigs:
                result.update(status="not_assessable", reason="missing_or_ambiguous_contig",
                              contigs=sorted(missing_contigs))
            else:
                try:
                    with reconstruction_deadline(seconds):
                        with pysam.AlignmentFile(acquired["bam"]["path"]) as bam:
                            reconstruction = reconstruct_one(bam, group_id, manifest["geometries"][group_id],
                                                             source, reference, orientation, parameters)
                        reconstruction["limitations"] = sorted(set(reconstruction["limitations"]) |
                                                               set(result["input_limitations"]))
                        reconstruction["upstream_acquisition"] = dict(
                            status=acquired["status"], receipt_sha256=identity(acquired),
                            limits=acquired.get("receipt", {}).get("limits", []),
                            failed_queries=acquired.get("receipt", {}).get("failed_queries", []))
                        summary = summarize(reconstruction)
                except TimeoutError as error:
                    result.update(status="reconstruction_timeout", error=str(error))
                else:
                    result.update(summary)
                    details = directory / "reconstructions" / source_id / path.name
                    write_json(details, reconstruction)
                    result["reconstruction"] = dict(path=str(details.relative_to(directory)), sha256=digest(details))
                    result["reference_exclusions"] = reference["assignments"][group_id]["excluded_transcripts"]
                    result["unsupported_reference_cds"] = sorted(
                        set(reference.get("unsupported_cds", {})) &
                        set(reference["assignments"][group_id]["transcripts"]))
            write_json(path, result)
            counts[result["status"]] = counts.get(result["status"], 0) + 1
    return dict(source=source_id, batch=batch_id, outcomes=counts)


def run(directory, cohort="tumor_candidate", workers=2, source_id=None, seconds=180):
    if not 1 <= seconds <= 3600:
        raise ValueError("Per-view time limit must be 1 through 3600 seconds")
    directory = Path(directory)
    manifest = load_inventory(directory)
    batches = read_json(directory / "batches.json")
    sources = selected_sources(manifest, cohort)
    if source_id:
        if source_id not in sources:
            raise ValueError("Source not in cohort")
        sources = {source_id: sources[source_id]}
    engine = implementation_identity()
    run_definition = dict(parameters=PARAMETERS, max_seconds=seconds, engine_id=engine,
                          isovar_version=isovar.__version__, cohort=cohort, sources=sorted(sources))
    write_json(directory / ("run-" + sha256(canonical(run_definition)).hexdigest()[:24] + ".json"), run_definition)
    tasks = [(sid, batch["id"]) for sid in sources for batch in batches
             if (directory / "sources" / sid / (batch["id"] + ".json")).exists()]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(run_batch, directory, sid, bid, PARAMETERS, seconds, engine) for sid, bid in tasks]
        for future in as_completed(futures):
            print(future.result(), flush=True)
