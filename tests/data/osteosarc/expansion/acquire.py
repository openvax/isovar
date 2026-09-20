"""Bounded header/index/region acquisition; never download whole BAMs.

Run from the repository root with a snapshot directory made by inventory.py.
Original headers and immutable request receipts are retained per source.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from tests.data.osteosarc.expansion.inventory import digest, write_json  # noqa: E402
from isovar.sid_data import open_dataset, extract_regions  # noqa: E402


ASSEMBLY_LENGTHS = {
    "GRCh38": {"1": 248956422, "2": 242193529, "3": 198295559, "X": 156040895},
    "GRCh37": {"1": 249250621, "2": 243199373, "3": 198022430, "X": 155270560},
}


def assembly_from_header(header):
    """Require at least two distinguishing contigs and no contradictory lengths."""
    sequences = {s["SN"].removeprefix("chr"): s["LN"] for s in header.get("SQ", [])}
    if not sequences:
        return None
    matches = []
    for assembly, lengths in ASSEMBLY_LENGTHS.items():
        present = set(lengths) & set(sequences)
        if len(present) >= 2 and all(sequences[c] == lengths[c] for c in present):
            matches.append(assembly)
    return matches[0] if len(matches) == 1 else None


def error_record(error, stage):
    return dict(status="error", stage=stage,
                error=dict(type=type(error).__name__, message=str(error)))


def survey_header(source, destination):
    """Inspect reference/processing metadata before constructing region queries."""
    directory = Path(destination) / "alignments" / source["source_id"]
    directory.mkdir(parents=True, exist_ok=True)
    result_path = directory / "header.json"
    header_path = directory / "header.sam"
    if result_path.exists():
        result = json.loads(result_path.read_text())
        if result["source_url"] != source["url"]:
            raise ValueError("Header receipt source mismatch")
        if result["status"] == "ok" and digest(header_path) != result["header_sha256"]:
            raise ValueError("Header snapshot checksum mismatch")
        return result
    attempts = []
    result = dict(source_id=source["source_id"], source_url=source["url"])
    for attempt in range(3):
        try:
            info = open_dataset().inspect_alignment(source["url"], timeout=90)
            raw_header = info.path.read_bytes()
            header = info.header
            assembly = assembly_from_header(header)
            if header_path.exists():
                if header_path.read_bytes() != raw_header:
                    raise ValueError("Existing header differs from source")
            else:
                with header_path.open("xb") as handle:
                    handle.write(raw_header)
            result.update(status="ok", header_sha256=sha256(raw_header).hexdigest(),
                          header=header, assembly=assembly, osteosarc=info.receipt,
                          coordinate_status=("genomic" if assembly else
                                             "unmapped_reference" if not header.get("SQ") else
                                             "unresolved_reference"))
            break
        except (subprocess.SubprocessError, ValueError) as error:
            detail = error_record(error, "header_acquisition")
            if isinstance(error, subprocess.CalledProcessError):
                detail["stderr"] = error.stderr.decode(errors="replace")
            attempts.append(detail)
    else:
        result.update(attempts[-1])
    result["failed_attempts"] = attempts
    write_json(result_path, result)
    return result


def regions_for_header(variants, header):
    """Resolve exact contig aliases, retaining explicit missing-contig outcomes."""
    names = {r["SN"] for r in header.get("SQ", [])}
    regions, missing = [], []
    for variant in variants:
        canonical = variant["chrom"].removeprefix("chr")
        aliases = ({"M", "MT", "chrM", "chrMT"} if canonical in ("M", "MT") else
                   {canonical, "chr" + canonical})
        candidates = names & aliases
        if not candidates:
            missing.append(variant["variant_id"])
            continue
        if len(candidates) != 1:
            raise ValueError(f"Ambiguous BAM contig identity: {sorted(candidates)}")
        contig = next(iter(candidates))
        regions.append(f"{contig}:{max(1, variant['pos'] - 1)}-{variant['pos'] + len(variant['ref'])}")
    return sorted(set(regions)), missing


def acquire_regions(source, variants, destination, assembly, label=None, timeout=600):
    """Fetch the indexed union once, preserving original record multiplicity."""
    directory = Path(destination) / "alignments" / source["source_id"]
    header = survey_header(source, destination)
    label = assembly if label is None else label
    if not label or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_" for c in label):
        raise ValueError("Invalid acquisition label")
    result_path = directory / f"regions-{label}.json"
    output = directory / f"regions-{label}.bam"
    request = dict(source_url=source["url"], assembly=assembly,
                   variants_sha256=sha256(json.dumps(variants, sort_keys=True).encode()).hexdigest())
    if result_path.exists():
        result = json.loads(result_path.read_text())
        if result["request"] != request:
            raise ValueError("Regional request drift")
        if result["status"] == "ok" and digest(output) != result["bam_sha256"]:
            raise ValueError("Regional BAM checksum mismatch")
        return result
    result = dict(source_id=source["source_id"], request=request, variant_ids=[v["variant_id"] for v in variants])
    if header["status"] != "ok":
        result.update(status="error", stage="header_acquisition", error=header["error"])
    elif header["assembly"] != assembly:
        return dict(result, status="not_requested", observed_assembly=header["assembly"])
    else:
        stage = "index_acquisition"
        try:
            mito = [v for v in variants if v["chrom"].removeprefix("chr") in ("M", "MT")]
            if mito:
                expected = 16571 if mito[0].get("genomic_validation", {}).get("genome") == "hg19" else 16569
                lengths = [r["LN"] for r in header["header"]["SQ"] if r["SN"] in ("M", "MT", "chrM", "chrMT")]
                if lengths and lengths != [expected]:
                    raise ValueError("Mitochondrial reference length does not match validated variant coordinate system")
            if not source["listed_indexes"]:
                raise ValueError("No genomic BAM index in pinned inventory; whole-BAM fetching is prohibited")
            regions, missing = regions_for_header(variants, header["header"])
            if not regions:
                raise ValueError("None of the requested contigs are in this alignment")
            result.update(regions=regions, missing_contigs_for=missing)
            stage = "regional_acquisition"
            lengths = {r["SN"]: r["LN"] for r in header["header"]["SQ"]}
            subset = extract_regions(source["url"], regions, assembly, output,
                                     reference_lengths=lengths, timeout=timeout)
            result.update(status="ok", bam_sha256=digest(output),
                          region_records=subset.receipt["records"], bam_bytes=output.stat().st_size,
                          failed_attempts=[], index_sha256=digest(str(output) + ".bai"),
                          osteosarc=subset.receipt)
        except (subprocess.SubprocessError, ValueError, OSError) as error:
            result.update(error_record(error, stage))
    write_json(result_path, result)
    return result


def survey(destination, workers=4, mode="headers", assembly="GRCh38", inventory_name="inventory.json",
           subset="all", retry_failed=False, source_id=None, timeout=600):
    inventory = json.loads((Path(destination) / inventory_name).read_text())
    sources = [r for r in inventory["alignments"] if r["scope"] in ("rna_candidate", "unresolved_assay")]
    if source_id:
        sources = [s for s in sources if s["source_id"] == source_id]
        if not sources:
            raise ValueError("Unknown/ineligible source identity")
    if retry_failed:
        selected = []
        for source in sources:
            receipt = Path(destination) / "alignments" / source["source_id"] / f"regions-{assembly}.json"
            if receipt.exists():
                old = json.loads(receipt.read_text())
                if old["status"] == "error" and old["stage"] == "regional_acquisition":
                    selected.append(source)
        sources = selected
    variants = inventory["variants"]
    if subset != "all":
        variants = [v for v in variants if (v["chrom"] == "chrM") == (subset == "mitochondrial")]
    if not variants or any(v["assembly"] != assembly for v in variants):
        raise ValueError("Variant list must be nonempty and explicitly match requested assembly")
    label = assembly if subset == "all" else f"{assembly}-{subset}"
    with ThreadPoolExecutor(max_workers=workers) as executor:
        if mode == "headers":
            futures = {executor.submit(survey_header, source, destination): source for source in sources}
        else:
            if assembly == "GRCh37" and any("coordinate_mapping" not in v or "genomic_validation" not in v for v in variants):
                raise ValueError("GRCh37 acquisition requires independently validated lifted variant definitions")
            futures = {executor.submit(acquire_regions, source, variants, destination, assembly, label, timeout): source
                       for source in sources if source["scope"] == "rna_candidate"}
        for future in as_completed(futures):
            source = futures[future]
            result = future.result()
            print(source["source_id"], result["status"], result.get("assembly", result.get("request", {}).get("assembly")),
                  source["key"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--mode", choices=("headers", "regions"), default="headers")
    parser.add_argument("--assembly", choices=tuple(ASSEMBLY_LENGTHS), default="GRCh38")
    parser.add_argument("--inventory", default="inventory.json")
    parser.add_argument("--subset", choices=("all", "nuclear", "mitochondrial"), default="all")
    parser.add_argument("--retry-failed", action="store_true")
    parser.add_argument("--source-id")
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    if not 1 <= args.workers <= 8:
        parser.error("Use 1 through 8 bounded acquisition workers")
    if not 1 <= args.timeout <= 1800:
        parser.error("Timeout must be 1 through 1800 seconds")
    survey(args.destination, args.workers, args.mode, args.assembly, args.inventory, args.subset,
           args.retry_failed, args.source_id, args.timeout)
