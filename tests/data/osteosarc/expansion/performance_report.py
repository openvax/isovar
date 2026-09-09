"""Package #227 full-source reruns and benchmark provenance, losslessly."""

import argparse
import json
from pathlib import Path

from tests.data.osteosarc.expansion.inventory import digest, write_json
from tests.data.osteosarc.expansion.report import compact_row, write_compressed


SOURCE_IDS = ["53f498a544883d51", "1120a096937e29e1", "2bc1fc291308debb",
              "58c4d68e5f7e55b5", "62a9d3e67e1b564e"]
MT_VARIANT = "MT_ND5-chrM-12994"


def read_json(path):
    return json.loads(Path(path).read_text())


def collect_reruns(audit, baseline):
    rows, runs, controls, checkpoints = [], {}, 0, {}
    for source_id in SOURCE_IDS:
        directory = Path(audit) / source_id
        runs[source_id] = read_json(directory / "run.json")
        paths = sorted(p for p in directory.glob("*.json") if p.name != "run.json")
        if len(paths) != 44:
            raise ValueError(f"Expected all 44 source loci: {source_id}")
        for path in paths:
            current = read_json(path)
            previous_path = Path(baseline) / source_id / path.name
            previous = read_json(previous_path)
            checkpoints[f"{source_id}/{path.name}"] = dict(current=digest(path), baseline=digest(previous_path))
            for mode in ("defaults", "primary_only"):
                if current[mode]["counts"] != previous[mode]["counts"]:
                    raise ValueError(f"Read counts changed: {path} {mode}")
                if current[mode]["status"] != "ok" or current[mode]["uncapped_validation_status"] != "ok":
                    raise ValueError(f"Unresolved rerun: {path} {mode}")
            if current["variant_id"] == MT_VARIANT:
                if any(previous[m]["outcome"] != "audit_timeout" for m in ("defaults", "primary_only")):
                    raise ValueError("Expected an originally interrupted mitochondrial row")
                if any(not current[m]["proteins"] for m in ("defaults", "primary_only")):
                    raise ValueError("Recovered mitochondrial modes must yield independently checked proteins")
                metadata = lambda row: {k: v for k, v in row.items() if k not in ("defaults", "primary_only", "software")}
                if metadata(current) != metadata(previous):
                    raise ValueError("Mitochondrial input/independent evidence changed")
                rows.append(compact_row(current))
            else:
                # Compare all completed nuclear results, including every
                # pre-cap protein/support set. Only software version changes.
                left, right = compact_row(current), compact_row(previous)
                left.pop("software")
                right.pop("software")
                if left != right:
                    raise ValueError(f"Completed nuclear audit changed: {path}")
                controls += 1
    return dict(rows=rows, runs=runs, checkpoints=checkpoints,
                unchanged_nuclear_rows=controls, unchanged_nuclear_modes=2 * controls,
                recovered_mitochondrial_modes=2 * len(rows))


def collect_benchmarks(entries, comparisons):
    runs = {}
    for entry in entries:
        label, directory = entry.split("=", 1)
        if label in runs:
            raise ValueError("Duplicate benchmark label")
        directory = Path(directory)
        run = dict(identity=read_json(directory / "identity.json"), result=read_json(directory / "result.json"))
        run["measurements"] = {name: read_json(directory / (name + ".json"))
                               for name in ("production", "validation") if (directory / (name + ".json")).exists()}
        if not run["measurements"]:  # Original profiling pilot wrote one combined file.
            run["measurements"] = read_json(directory / "measurements.json")
        run["profiles"] = {name: (directory / (name + ".txt")).read_text()
                           for name in ("production", "validation") if (directory / (name + ".txt")).exists()}
        run["artifacts"] = {p.name: digest(p) for p in sorted(directory.iterdir()) if p.is_file()}
        runs[label] = run
    checked = []
    for pair in comparisons:
        before, after = pair.split("=", 1)
        left, right = runs[before], runs[after]
        for key in ("input_bam_sha256", "input_index_sha256", "variant_id", "mode", "reference_manifest_sha256"):
            if left["identity"][key] != right["identity"][key]:
                raise ValueError(f"Different comparison inputs: {pair} {key}")
        for run in (left, right):
            if run["result"]["status"] != "ok" or run["result"]["validation_status"] != "ok":
                raise ValueError(f"Incomplete comparison: {pair}")
        for key in ("rna_sha256", "proteins_sha256", "counts", "public_protein_count", "checked_protein_count"):
            if left["result"][key] != right["result"][key]:
                raise ValueError(f"Results changed: {pair} {key}")
        checked.append(dict(before=before, after=after, exact_rna_and_ranked_proteins=True))
    return dict(runs=runs, comparisons=checked)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", required=True)
    parser.add_argument("--baseline-audit", required=True)
    parser.add_argument("--benchmark", action="append", default=[], help="label=directory")
    parser.add_argument("--compare", action="append", default=[], help="before-label=after-label")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    reruns = collect_reruns(args.audit, args.baseline_audit)
    write_compressed(output / "reruns.json.gz", reruns)
    write_json(output / "benchmarks.json", collect_benchmarks(args.benchmark, args.compare))
    write_json(output / "manifest.json", dict(
        files={p.name: digest(p) for p in sorted(output.iterdir())},
        generator_sha256=digest(__file__),
        unchanged_nuclear_rows=reruns["unchanged_nuclear_rows"],
        recovered_mitochondrial_modes=reruns["recovered_mitochondrial_modes"]))
