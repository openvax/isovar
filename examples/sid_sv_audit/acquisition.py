"""Indexed regional acquisition through osteosarc, with original-record recovery."""

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
import shutil
import subprocess

from osteosarc import Cache, File, RecoveryPolicy, Region, extract_reads, inspect_alignment, read_receipt_files
from osteosarc.models import SampleClaim
from osteosarc.errors import OsteosarcError, RecordLimitError

from tests.data.osteosarc.expansion.acquire import assembly_from_header
from .inventory import digest, identity, load_inventory, read_json, write_json
from . import references

USABLE_STATUSES = ("acquired", "incomplete", "truncated")


def receipt_summary(receipt):
    """Keep reporting fields inline; pin the complete upstream provenance once."""
    return {key: receipt[key] for key in ("status", "scope", "records", "complete_template",
                                        "limits", "failed_queries") if key in receipt}


def upstream_provenance(subset):
    path = subset.receipt_path
    return [dict(path=str(p), sha256=checksum) for p, checksum in
            read_receipt_files(path, path.parents[2]).items()]


def verify_upstream(acquired):
    for pin in acquired.get("upstream_receipt_files", []):
        if digest(pin["path"]) != pin["sha256"]:
            raise ValueError("Upstream acquisition provenance checksum mismatch")


def input_limitations(acquired):
    """Keep upstream limits visible even when a usable BAM was returned."""
    flags = {"input_acquisition_" + acquired["status"]} if acquired["status"] != "acquired" else set()
    flags.update("input_acquisition_limit:" + limit for limit in acquired.get("receipt", {}).get("limits", []))
    if acquired.get("unavailable_contigs"):
        flags.add("input_acquisition_missing_contigs")
    return sorted(flags)


def source_file(source):
    fields = {key: source[key] for key in File.__dataclass_fields__}
    fields["claims"] = tuple(SampleClaim(**row) for row in fields["claims"])
    fields["index_urls"] = tuple(fields["index_urls"])
    return File(**fields)


def selected_sources(manifest, cohort):
    if cohort not in ("tumor_candidate", "blood_control", "unresolved_specimen"):
        raise ValueError("Choose an explicit RNA cohort")
    return {key: row for key, row in manifest["sources"].items() if row["selection"]["cohort"] == cohort}


def header(directory, source, cache, timeout=90):
    directory = Path(directory) / "sources" / source["id"]
    path = directory / "header.json"
    request = identity(source)
    if path.exists():
        result = read_json(path)
        if result["source_identity"] != request:
            raise ValueError("Source drift in saved header")
        return result
    result = dict(source_identity=request, source_id=source["id"])
    if not source["index_urls"]:
        result.update(status="no_listed_index")
    elif source["format"] != "bam":
        result.update(status="unsupported_alignment_format")
    else:
        try:
            info = inspect_alignment(source_file(source), cache=cache, timeout=timeout)
            assembly = assembly_from_header(info.header)
            result.update(status="ready" if assembly == "GRCh38" else "reference_unavailable",
                          assembly=assembly, header=info.header, receipt=info.receipt)
        except (OSError, ValueError, OsteosarcError, subprocess.SubprocessError) as error:
            result.update(status="header_error", error=dict(type=type(error).__name__, message=str(error)))
    write_json(path, result)
    return result


def survey(directory, cache, cohort="tumor_candidate", workers=4):
    manifest = load_inventory(directory)
    cache = Cache(cache)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(header, directory, source, cache): source
                   for source in selected_sources(manifest, cohort).values()}
        for future in as_completed(futures):
            source = futures[future]
            result = future.result()
            print(source["id"][:12], result["status"], source["key"], flush=True)


def merge_regions(regions):
    by_contig = defaultdict(list)
    for contig, start, end in regions:
        by_contig[contig].append((start, end))
    result = []
    for contig, intervals in sorted(by_contig.items()):
        merged = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
            else:
                merged.append((start, end))
        result.extend((contig, a, b) for a, b in merged)
    return result


def make_batch(manifest, reference, group_ids, window):
    regions = []
    for name in group_ids:
        for end in manifest["geometries"][name]["breakends"]:
            regions.append((end["contig"], max(0, end["position"] - window), end["position"] + window))
        for tid in reference["assignments"][name]["transcripts"]:
            model = reference["models"][tid]
            regions.extend((model["contig"], a, b) for a, b in model["exons"])
    batch = dict(geometries=group_ids, regions=[list(r) for r in merge_regions(regions)], window=window)
    return dict(batch, id=identity(batch)[:24])


def planned_batches(manifest, reference, batch_size=50, window=2000):
    """Batch overlapping requests, without capping or selecting supporting reads."""
    names = sorted(manifest["geometries"])
    return [make_batch(manifest, reference, names[start:start + batch_size], window)
            for start in range(0, len(names), batch_size)]


def acquisition_tree(directory, source_id, batch_id):
    """Verify an exact partition; missing or changed children cannot disappear."""
    leaves, receipts = {}, {}

    def visit(bid):
        if bid in receipts or Path(bid).name != bid:
            raise ValueError("Repeated or invalid acquisition batch")
        result = read_json(Path(directory) / "sources" / source_id / (bid + ".json"))
        if result["source_id"] != source_id or result["batch_id"] != bid:
            raise ValueError("Acquisition is filed under the wrong identity")
        verify_upstream(result)
        receipts[bid] = result
        if result["status"] != "partitioned":
            leaves[bid] = result
            return result
        request, observed, regions = result["request"], [], []
        for child in result["children"]:
            saved = visit(child["batch"]["id"])
            if identity(saved) != child["receipt_sha256"]:
                raise ValueError("Partition child checksum mismatch")
            if saved["request"]["batch"] != child["batch"]:
                raise ValueError("Partition child batch changed")
            if dict(saved["request"], batch=request["batch"]) != request:
                raise ValueError("Partition child source or acquisition policy changed")
            observed.extend(child["batch"]["geometries"])
            regions.extend(child["batch"]["regions"])
            if child["batch"]["window"] != request["batch"]["window"]:
                raise ValueError("Partition changed event window")
        if len(observed) != len(set(observed)) or sorted(observed) != sorted(request["batch"]["geometries"]):
            raise ValueError("Partition does not cover every geometry exactly once")
        if merge_regions(regions) != merge_regions(request["batch"]["regions"]):
            raise ValueError("Partition changed requested context")
        return result

    visit(batch_id)
    return leaves, receipts


def acquire_batch(directory, source, batch, cache, policy, timeout=300,
                  manifest=None, reference=None, min_free_gib=0):
    """Keep context and verified mate/SA partners; limits stay in the receipt."""
    directory = Path(directory)
    request = dict(source_identity=identity(source), batch=batch, policy=asdict(policy), timeout=timeout,
                   partition_on_record_limit=manifest is not None)
    path = directory / "sources" / source["id"] / (batch["id"] + ".json")
    if path.exists():
        result = read_json(path)
        if result["request"] != request:
            raise ValueError("Acquisition request drift")
        verify_upstream(result)
        if result["status"] in USABLE_STATUSES:
            for key in ("bam", "index"):
                if digest(result[key]["path"]) != result[key]["sha256"]:
                    raise ValueError("Acquired input checksum mismatch")
        elif result["status"] == "partitioned":
            leaves, _ = acquisition_tree(directory, source["id"], batch["id"])
            for leaf in leaves.values():
                for key in ("bam", "index"):
                    if key in leaf and digest(leaf[key]["path"]) != leaf[key]["sha256"]:
                        raise ValueError("Acquired input checksum mismatch")
        return result
    if cache is not None and min_free_gib and shutil.disk_usage(cache.root).free < min_free_gib * 1024 ** 3:
        # No receipt: this is pending work, not an unavailable specimen.
        raise OSError("Acquisition paused: cache has less than %s GiB free" % min_free_gib)
    info = header(directory, source, cache)
    result = dict(request=request, source_id=source["id"], batch_id=batch["id"])
    if info["status"] != "ready":
        result.update(status=info["status"])
    else:
        lengths = {row["SN"]: row["LN"] for row in info["header"]["SQ"]}
        aliases, unavailable = {}, []
        for contig, _, _ in batch["regions"]:
            matches = {contig, contig.removeprefix("chr"), "chr" + contig.removeprefix("chr")} & lengths.keys()
            if len(matches) != 1:
                unavailable.append(contig)
            else:
                aliases[contig] = next(iter(matches))
        unavailable.extend(c for c, a, _ in batch["regions"] if c in aliases and a >= lengths[aliases[c]])
        # A missing contig must not prevent acquisition for the other targets.
        result["unavailable_contigs"] = sorted(set(unavailable))
        regions = [Region(aliases[c], a, min(b, lengths[aliases[c]]), "GRCh38")
                   for c, a, b in batch["regions"] if c in aliases and a < lengths[aliases[c]]]
        if not regions:
            result.update(status="missing_or_ambiguous_contig")
        else:
            try:
                subset = extract_reads(source_file(source), regions, cache=cache,
                                       recovery=policy, timeout=timeout)
                status = subset.receipt.get("status")
                result.update(status=status if status in ("incomplete", "truncated") else "acquired",
                              receipt=receipt_summary(subset.receipt), upstream_receipt_files=upstream_provenance(subset),
                              contig_aliases=aliases,
                              bam=dict(path=str(subset.path), sha256=digest(subset.path)),
                              index=dict(path=str(subset.index_path), sha256=digest(subset.index_path)))
            except (OSError, ValueError, OsteosarcError, subprocess.SubprocessError) as error:
                if isinstance(error, RecordLimitError) and manifest is not None and len(batch["geometries"]) > 1:
                    names = batch["geometries"]
                    midpoint = len(names) // 2
                    children = []
                    for group_ids in (names[:midpoint], names[midpoint:]):
                        child = make_batch(manifest, reference, group_ids, batch["window"])
                        saved = acquire_batch(directory, source, child, cache, policy, timeout,
                                              manifest, reference, min_free_gib)
                        children.append(dict(batch=child, receipt_sha256=identity(saved)))
                    result.update(status="partitioned", reason="record_limit", children=children)
                else:
                    result.update(status="acquisition_error", error=dict(type=type(error).__name__, message=str(error)))
                    stderr = getattr(error, "stderr", None)
                    if stderr:
                        result["error"]["stderr"] = (stderr.decode(errors="replace") if isinstance(stderr, bytes)
                                                     else str(stderr))
    write_json(path, result)
    return result


def acquire(directory, cache, cohort="tumor_candidate", workers=2, source_id=None,
            partner_batch_size=64, max_partner_queries=128, partner_timeout=30,
            batch_id=None, min_free_gib=8):
    if min_free_gib < 0:
        raise ValueError("Minimum free disk space cannot be negative")
    manifest, reference = load_inventory(directory), references.load(directory)
    batches = planned_batches(manifest, reference)
    write_json(Path(directory) / "batches.json", batches)
    if batch_id:
        batches = [batch for batch in batches if batch["id"] == batch_id]
        if not batches:
            raise ValueError("Unknown acquisition batch")
    sources = selected_sources(manifest, cohort)
    if source_id:
        if source_id not in sources:
            raise ValueError("Source not in selected cohort")
        sources = {source_id: sources[source_id]}
    policy = RecoveryPolicy(max_rounds=4, max_intervals=20000, max_bases=20_000_000,
                            max_records=500_000, on_timeout="incomplete",
                            partner_batch_size=partner_batch_size, max_partner_queries=max_partner_queries,
                            partner_timeout=partner_timeout)
    cache = Cache(cache)
    cache.root.mkdir(parents=True, exist_ok=True)
    # A source's batches run serially: bounded memory and no shared-index races.
    def run_source(source):
        for batch in batches:
            result = acquire_batch(directory, source, batch, cache, policy,
                                   manifest=manifest, reference=reference, min_free_gib=min_free_gib)
            print(source["id"][:12], batch["id"], result["status"], flush=True)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(run_source, source) for source in sources.values()]
        for future in as_completed(futures):
            future.result()
