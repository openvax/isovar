"""Indexed regional acquisition through osteosarc, with original-record recovery."""

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
import subprocess

from osteosarc import Cache, File, RecoveryPolicy, Region, extract_reads, inspect_alignment
from osteosarc.models import SampleClaim
from osteosarc.errors import OsteosarcError

from tests.data.osteosarc.expansion.acquire import assembly_from_header
from .inventory import digest, identity, load_inventory, read_json, write_json
from . import references


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


def planned_batches(manifest, reference, batch_size=50, window=2000):
    """Batch overlapping requests, without capping or selecting supporting reads."""
    names = sorted(manifest["geometries"])
    batches = []
    for start in range(0, len(names), batch_size):
        group_ids = names[start:start + batch_size]
        regions = []
        for name in group_ids:
            for end in manifest["geometries"][name]["breakends"]:
                regions.append((end["contig"], max(0, end["position"] - window), end["position"] + window))
            for tid in reference["assignments"][name]["transcripts"]:
                model = reference["models"][tid]
                regions.extend((model["contig"], a, b) for a, b in model["exons"])
        batch = dict(geometries=group_ids, regions=[list(r) for r in merge_regions(regions)], window=window)
        batches.append(dict(batch, id=identity(batch)[:24]))
    return batches


def acquire_batch(directory, source, batch, cache, policy, timeout=300):
    """Keep context and verified mate/SA partners; limits stay in the receipt."""
    directory = Path(directory)
    request = dict(source_identity=identity(source), batch=batch, policy=asdict(policy))
    path = directory / "sources" / source["id"] / (batch["id"] + ".json")
    if path.exists():
        result = read_json(path)
        if result["request"] != request:
            raise ValueError("Acquisition request drift")
        if result["status"] in ("acquired", "incomplete"):
            for key in ("bam", "index"):
                if digest(result[key]["path"]) != result[key]["sha256"]:
                    raise ValueError("Acquired input checksum mismatch")
        return result
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
        # Do not send one bad contig to samtools and silently lose a batch.
        if unavailable:
            result.update(status="missing_or_ambiguous_contig", contigs=sorted(set(unavailable)))
        else:
            regions = [Region(aliases[c], a, min(b, lengths[aliases[c]]), "GRCh38")
                       for c, a, b in batch["regions"]]
            try:
                subset = extract_reads(source_file(source), regions, cache=cache,
                                       recovery=policy, timeout=timeout)
                result.update(status="incomplete" if subset.receipt.get("status") == "incomplete" else "acquired",
                              receipt=subset.receipt, contig_aliases=aliases,
                              bam=dict(path=str(subset.path), sha256=digest(subset.path)),
                              index=dict(path=str(subset.index_path), sha256=digest(subset.index_path)))
            except (OSError, ValueError, OsteosarcError, subprocess.SubprocessError) as error:
                result.update(status="acquisition_error", error=dict(type=type(error).__name__, message=str(error)))
    write_json(path, result)
    return result


def acquire(directory, cache, cohort="tumor_candidate", workers=2, source_id=None):
    manifest, reference = load_inventory(directory), references.load(directory)
    batches = planned_batches(manifest, reference)
    write_json(Path(directory) / "batches.json", batches)
    sources = selected_sources(manifest, cohort)
    if source_id:
        if source_id not in sources:
            raise ValueError("Source not in selected cohort")
        sources = {source_id: sources[source_id]}
    policy = RecoveryPolicy(max_rounds=4, max_intervals=20000, max_bases=20_000_000,
                            max_records=500_000, on_timeout="incomplete")
    cache = Cache(cache)
    # A source's batches run serially: bounded memory and no shared-index races.
    def run_source(source):
        for batch in batches:
            result = acquire_batch(directory, source, batch, cache, policy)
            print(source["id"][:12], batch["id"], result["status"], flush=True)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(run_source, source) for source in sources.values()]
        for future in as_completed(futures):
            future.result()
