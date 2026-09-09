"""Assemble a complete, explicitly missing-aware variant/RNA-product matrix.

This only reads completed checkpoints and original acquisition receipts.
It never imputes counts from a different product, assembly, or fixture.
"""

import argparse
from collections import Counter, defaultdict
from copy import deepcopy
import gzip
from hashlib import sha256
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from tests.data.osteosarc.expansion.inventory import digest, write_json  # noqa: E402


def json_bytes(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def write_compressed(path, value):
    with Path(path).open("xb") as handle:
        with gzip.GzipFile(filename="", fileobj=handle, mode="wb", mtime=0) as archive:
            archive.write(json_bytes(value))


def compact_row(row):
    """Losslessly intern support names shared by many ranked alternatives."""
    row = deepcopy(row)
    proteins = [p for mode in ("defaults", "primary_only") for field in ("proteins", "uncapped_ranked_proteins")
                for p in row.get(mode, {}).get(field, [])]
    names = sorted({name for p in proteins for name in p["supporting_read_names"]})
    indices = {name: index for index, name in enumerate(names)}
    row["supporting_read_name_table"] = names
    for protein in proteins:
        support = sorted(protein.pop("supporting_read_names"))
        protein["supporting_read_name_indices"] = [indices[name] for name in support]
        protein["supporting_read_names_sha256"] = sha256(json_bytes(support)).hexdigest()
        protein["checks"].sort(key=json_bytes)
    return row


def protein_supporting_read_names(row, protein):
    """Decode the lossless name dictionary and validate the recorded set hash."""
    indices = protein["supporting_read_name_indices"]
    table = row["supporting_read_name_table"]
    if any(type(i) is not int or not 0 <= i < len(table) for i in indices) or indices != sorted(set(indices)):
        raise ValueError("Invalid protein support-name indices")
    names = [table[i] for i in indices]
    if sha256(json_bytes(names)).hexdigest() != protein["supporting_read_names_sha256"]:
        raise ValueError("Protein support-name checksum mismatch")
    return names


def source_disposition(source, header, catalogue):
    """Refine conservative path classification only with explicit evidence."""
    category = catalogue.get(source["key"])
    if category in ("WGS", "WES"):
        return "excluded_dna_alignment", f"BAM catalogue category {category}"
    if header and header.get("status") == "ok":
        programs = header["header"].get("PG", [])
        if source["key"].startswith("kamil/basespace/results/") and any(
                "DRAGEN SW" in p.get("ID", "") and re.search(r"Mapper\.rna-mode:\s*0(?:\s|$)", p.get("CL", ""))
                for p in programs):
            return "excluded_dna_alignment", "Original DRAGEN alignment header: Mapper.rna-mode: 0; WGS processing family"
        sequences = header["header"].get("SQ", [])
        if sequences and all(s["SN"].startswith("HLA") for s in sequences):
            return "excluded_receptor_or_targeted_alignment", "Original header contains only HLA reference sequences"
    return source["scope"], "Conservative inventory path/category rule; unresolved assays remain unresolved"


def missing_row(source, variant, header, acquisitions):
    """Absent computation is never biological zero, even after a good download."""
    vid = variant["variant_id"]
    applicable = [r for r in acquisitions if vid in r.get("variant_ids", [vid])]
    if header is None or header.get("status") != "ok":
        state = "header_unavailable"
    elif header.get("coordinate_status") == "unmapped_reference":
        state = "no_genomic_coordinates"
    elif header.get("assembly") not in ("GRCh37", "GRCh38"):
        state = "unresolved_reference"
    elif any(r["status"] == "ok" for r in applicable):
        # Important: the report cannot declare completion with this state.
        state = "audit_not_completed"
    elif not source.get("listed_indexes"):
        state = "missing_index"
    elif applicable:
        state = "acquisition_failed"
    else:
        state = "acquisition_not_attempted"
    return dict(source_id=source["source_id"], variant_id=vid, status=state,
                defaults=dict(status="unavailable", outcome=state, counts=None),
                primary_only=dict(status="unavailable", outcome=state, counts=None),
                independent_primary=None)


def summarize(rows, variants):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["variant_id"]].append(row)
    results = []
    for variant in variants:
        records = grouped[variant["variant_id"]]
        modes = {}
        for mode in ("defaults", "primary_only"):
            outcomes = Counter(r.get(mode, {}).get("outcome", r["status"]) for r in records)
            with_counts = [r for r in records if r.get(mode, {}).get("counts") is not None]
            alt = [r for r in with_counts if r[mode]["counts"]["reads"]["alt"] > 0]
            proteins = [r for r in records if r.get(mode, {}).get("proteins")]
            modes[mode] = dict(outcomes=dict(sorted(outcomes.items())), products_with_counts=len(with_counts),
                               products_with_alt=len(alt), products_with_ranked_protein=len(proteins),
                               products_with_any_expected_protein=sum(any(p["matches_expected"] for p in r[mode]["proteins"])
                                                                     for r in proteins),
                               products_with_mutation_containing_peptide_window=sum(
                                   any(p.get("mutation_containing_peptide_windows", 0) > 0 for p in r[mode]["proteins"])
                                   for r in proteins),
                               products_with_uncapped_validation_error=sum(
                                   r.get(mode, {}).get("uncapped_validation_status") == "error" for r in records),
                               max_alt_alignment_count=max((r[mode]["counts"]["reads"]["alt"] for r in alt), default=0),
                               distinct_top_sequences=sorted({r[mode]["proteins"][0]["amino_acids"] for r in proteins}))
        results.append(dict(variant_id=variant["variant_id"], gene=variant["gene"], **modes))
    return results


def collect(destination, audit_names, allow_incomplete=False):
    destination = Path(destination)
    inventory = json.loads((destination / "inventory-live-v2.json").read_text())
    metadata = json.loads((destination / "source-metadata.json").read_text())
    catalogue = json.loads((destination / "bams.json").read_text())
    categories = {b["url"]: c["name"] for c in catalogue["categories"] for b in c["bams"]}
    rows_by_key, runs = {}, {}
    for name in audit_names:
        directory = destination / name
        if not directory.is_dir():
            raise ValueError(f"Audit directory missing: {name}")
        for run_path in sorted(directory.glob("*/run.json")):
            identity = json.loads(run_path.read_text())
            run_id = sha256(json_bytes(identity)).hexdigest()
            runs[run_id] = identity
            for path in sorted(run_path.parent.glob("*.json")):
                if path.name == "run.json":
                    continue
                row = json.loads(path.read_text())
                key = row["source_id"], row["variant_id"]
                if key in rows_by_key:
                    raise ValueError(f"Duplicate source/variant audit: {key}")
                if row.get("source_bam_sha256") != identity["source_bam_sha256"]:
                    raise ValueError("Row/run source identity mismatch")
                rows_by_key[key] = dict(compact_row(row), audit_run=run_id)
    variants = inventory["variants"]
    sources, rows, provenance = [], [], {}
    for source in inventory["alignments"]:
        sid = source["source_id"]
        directory = destination / "alignments" / sid
        header_path = directory / "header.json"
        header = json.loads(header_path.read_text()) if header_path.exists() else None
        disposition, reason = source_disposition(source, header, categories)
        item = dict(source, scope=disposition, disposition_evidence=reason,
                    published_md5_claim=metadata["alignment_md5_claims"].get(sid))
        item["published_md5_verified_against_full_remote_file"] = False
        if header:
            if header["status"] == "ok" and digest(directory / "header.sam") != header["header_sha256"]:
                raise ValueError("Original header integrity failure")
            item["header"] = {k: v for k, v in header.items() if k != "header"}
            if header["status"] == "ok":
                # Vendor CL payloads can include credential-like configuration.
                # Publish identifiers/versions and a checksum, not arbitrary CL.
                item["programs"] = [dict({k: p[k] for k in ("ID", "PN", "VN", "PP") if k in p},
                                         command_sha256=sha256(p.get("CL", "").encode()).hexdigest())
                                    for p in header["header"].get("PG", [])]
                item["read_groups"] = header["header"].get("RG", [])
                item["reference_sequences"] = header["header"].get("SQ", [])
        acquisitions = []
        for path in sorted(directory.glob("regions-*.json")):
            if "stress" in path.name:
                continue
            receipt = json.loads(path.read_text())
            if receipt["status"] == "ok":
                if digest(path.with_suffix(".bam")) != receipt["bam_sha256"]:
                    raise ValueError("Full regional source integrity failure")
            acquisitions.append(receipt)
        item["acquisition_receipts"] = acquisitions
        sources.append(item)
        if disposition == "rna_candidate":
            provenance[sid] = dict(header_receipt_sha256=digest(header_path) if header else None)
            for variant in variants:
                key = sid, variant["variant_id"]
                rows.append(rows_by_key.pop(key) if key in rows_by_key else missing_row(source, variant, header, acquisitions))
    if rows_by_key:
        raise ValueError("Audit contains variants or sources outside the resolved RNA inventory")
    incomplete = [r for r in rows if r["status"] in ("audit_not_completed", "acquisition_not_attempted", "header_unavailable")]
    if incomplete and not allow_incomplete:
        raise ValueError(f"{len(incomplete)} rows have unfinished audit/acquisition; use --allow-incomplete only for diagnostics")
    payload = dict(schema_version=1, variants=variants, sources=sources, rows=rows, audit_runs=runs,
                   support_names_encoding="Each row has a sorted supporting_read_name_table; each protein indexes its exact subset",
                   summary=summarize(rows, variants), source_metadata=metadata,
                   inventory_metadata={k: v for k, v in inventory.items() if k not in ("variants", "alignments")},
                   inputs={n: digest(destination / n) for n in ("inventory-live-v2.json", "source-metadata.json", "bams.json",
                                                               "mitochondrial-diagnostics.json")},
                   mitochondrial_diagnostics=json.loads((destination / "mitochondrial-diagnostics.json").read_text()),
                   provenance=provenance, incomplete_row_count=len(incomplete),
                   scope="alignment products, not independent patients; no cross-product counts are summed")
    # Preserve original messages/receipts but replace this machine's cache root.
    return json.loads(json.dumps(payload).replace(str(destination), "<SOURCE_CACHE>"))


def markdown(payload):
    rows, sources = payload["rows"], payload["sources"]
    rna = [s for s in sources if s["scope"] == "rna_candidate"]
    outcomes = Counter(r.get("defaults", {}).get("outcome", r["status"]) for r in rows)
    lines = ["# Full-region vaccine RNA/protein audit", "",
             f"{len(payload['variants'])} vaccine variants × {len(rna)} RNA alignment products = {len(rows)} rows.", "",
             "These are products, **not independent samples**. Technical reprocessings and unresolved pooled/blood libraries are not combined.", "",
             "Counts below count source/variant rows, never pooled biological read totals. Missing acquisition and interrupted computation are not zero support.", "",
             "## Default-path outcomes", "", "| Outcome | Source/variant rows |", "|---|---:|"]
    lines += [f"| {name} | {count} |" for name, count in sorted(outcomes.items())]
    lines += ["", "## Per-variant coverage", "",
              "Each count is the number of RNA products with that result; one product can represent pooled donors or a duplicate processing.", "",
              "| Variant | Products with counts | With alt | With ranked protein | With any expected protein | With mutant 25-mer |",
              "|---|---:|---:|---:|---:|---:|"]
    for item in payload["summary"]:
        d = item["defaults"]
        lines.append(f"| {item['variant_id']} | {d['products_with_counts']} | {d['products_with_alt']} | "
                     f"{d['products_with_ranked_protein']} | {d['products_with_any_expected_protein']} | "
                     f"{d['products_with_mutation_containing_peptide_window']} |")
    lines += ["", "## Reading the machine report", "",
              "`matrix.json.gz` includes every source disposition, assembly/header/program/read-group evidence, acquisition receipts, "
              "all ranked amino-acid sequences, mutation intervals, support units, independent frame/stop/sequence checks, "
              "default and primary-only results, exact-CIGAR quality sensitivity, and mitochondrial caveats.", "",
              "A validated protein is a local RNA-derived window, not a reconstructed full-length protein or vaccine eligibility claim. "
              "`different_top_protein` means correct translation of RNA that differs from the reference-plus-nominated-edit expectation; "
              "inspect all checks and adjacent alleles. `protein_validation_error` is a distinct failure.", "",
              "The original RNA BAMs are checksum-verified, full expanded-locus extractions. The offline corpus is deliberately selected, "
              "and its counts must not replace this matrix. Full supporting-name sets are retained losslessly via each row's "
              "`supporting_read_name_table` and each protein's `supporting_read_name_indices`, with subset SHA256. "
              "`report.protein_supporting_read_names(row, protein)` decodes and verifies them.", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--audit", action="append", required=True)
    parser.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args()
    payload = collect(args.destination, args.audit, args.allow_incomplete)
    args.output.mkdir(parents=True, exist_ok=False)
    write_compressed(args.output / "matrix.json.gz", payload)
    with (args.output / "SUMMARY.md").open("x") as handle:
        handle.write(markdown(payload))
    write_json(args.output / "manifest.json", dict(files={p.name: digest(p) for p in sorted(args.output.iterdir())},
               inputs=payload["inputs"], row_count=len(payload["rows"]), incomplete_row_count=payload["incomplete_row_count"],
               generator_sources={p.name: digest(p) for p in sorted(Path(__file__).parent.glob("*.py"))}))
    print("REPORT", len(payload["rows"]), "rows; incomplete", payload["incomplete_row_count"])
