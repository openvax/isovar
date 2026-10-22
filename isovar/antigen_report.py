"""Compose distinct Isovar antigen-source evidence without mixing units.

Inputs are versioned outputs from the source-specific Isovar workflows. A
missing input remains unassessed, and reference annotations are never counted
as observed protein expression.
"""

import csv
import gzip
import json
from pathlib import Path

from .evidence_report import variant_type


SCHEMAS = {
    "sv": "isovar.sv_rna_candidates.v5",
    "fusion": "isovar.fusion_rna.v3",
    "viral": "isovar.viral_sequences.v1",
    "isoforms": "isovar.antigen_isoforms.v1",
}
CATEGORIES = ("small_variants", "sv", "fusion", "viral", "isoforms")
EVENT_LINKED = {"breakpoint_junction", "event_compatible_junction"}


def load_json(path, schema, sample_id):
    """Read a versioned Isovar JSON export, rejecting a different sample."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        result = json.load(handle)
    if not isinstance(result, dict) or result.get("schema") != schema:
        raise ValueError(f"{path}: expected schema {schema}")
    if result.get("sample_id") != sample_id:
        raise ValueError(f"{path}: sample_id does not match {sample_id}")
    return result


def load_small_variants(path, sample_id):
    """Read the per-call TSV from `isovar evidence-report`."""
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        needed = {"sample_id", "call_id", "ref", "alt", "exon_overlap",
                  "alt_reads", "protein_count", "orfdata_coding_block_overlap"}
        if not needed.issubset(reader.fieldnames or ()):
            raise ValueError(f"{path}: not an Isovar evidence-report per-call TSV")
        rows = list(reader)
    if any(row["sample_id"] != sample_id for row in rows):
        raise ValueError(f"{path}: sample_id does not match {sample_id}")
    if len({row["call_id"] for row in rows}) != len(rows):
        raise ValueError(f"{path}: duplicate small-variant call_id")
    return rows


def small_variant_summary(rows):
    """Summarize literal DNA calls and original RNA ALT support."""
    exonic = [row for row in rows if row["exon_overlap"] == "True"]
    by_class = {}
    for row in exonic:
        by_class.setdefault(variant_type(row), []).append(int(row["alt_reads"]))
    return dict(input_records=len(rows), exon_involved=len(exonic),
                rna_alt_observed=sum(int(row["alt_reads"]) > 0 for row in rows),
                local_products=sum(int(row["protein_count"]) > 0 for row in rows),
                ncorf_reference_candidates=sum(
                    row["orfdata_coding_block_overlap"] == "True" for row in rows),
                exonic_reads_by_class=by_class)


def _direct_fragments(path, *, relations=None):
    return max((int(junction.get("direct_support", {}).get("fragments", 0))
                for junction in path.get("junctions", ())
                if relations is None or junction.get("relation") in relations), default=0)


def sv_summary(result):
    """Summarize one nominated DNA SV without asserting tumor protein expression."""
    paths = result["paths"]
    linked = [path for path in paths if path.get("event_linkage", {}).get("status")
              in EVENT_LINKED]
    supported = [path for path in linked
                 if _direct_fragments(path, relations=EVENT_LINKED) > 0]
    translated = [path for path in supported
                  if any(EVENT_LINKED.intersection(row.get("departure_relations", ()))
                         for row in path.get("translations", ()))]
    exploratory = sum(len(path.get("exploratory_orfs", {}).get("candidates", ()))
                      for path in paths)
    return dict(event_id=result["event_id"], source=result.get("source", ""),
                input_events=1, candidate_paths=len(paths),
                linked_paths=len(linked),
                max_direct_junction_fragments=max(
                    (_direct_fragments(path, relations=EVENT_LINKED) for path in linked),
                    default=0),
                rna_linked_events=bool(supported),
                event_linked_translated_paths=len(translated),
                exploratory_orf_candidates=exploratory,
                status=result["status"])


def fusion_summary(result):
    """Summarize one supplied fusion transcript and its validated junction."""
    path = result["paths"][0]
    fragments = _direct_fragments(path)
    translated = bool(path.get("translations") and result["status"] in
                      {"translated", "ambiguous"} and fragments > 0)
    return dict(event_id=result["event_id"],
                source=result.get("provenance", {}).get("source", ""),
                input_transcripts=1, direct_junction_fragments=fragments,
                rna_supported=result["status"] != "insufficient_support" and fragments > 0,
                translated_candidates=len(path.get("translations", ())) if translated else 0,
                status=result["status"])


def viral_summary(result):
    """Summarize the selected viral panel with fragment and coding units."""
    contigs = result["contigs"]
    sequences = result["sequences"]
    coding = [coding for row in sequences for coding in row.get("coding", ())
              if coding.get("amino_acids")]
    return dict(selected_contigs=len(contigs),
                aligned_contigs=sum(int(row["aligned_fragments"]) > 0 for row in contigs),
                assigned_contigs=sum(int(row["sequence_assigned_fragments"]) > 0
                                     for row in contigs),
                observed_paths=len(sequences), coding_predictions=len(coding),
                complete_cds_predictions=sum(row.get("status") == "complete_CDS"
                                             for row in coding),
                contigs=[{key: row[key] for key in ("contig", "aligned_fragments",
                        "sequence_assigned_fragments", "status")} for row in contigs],
                source=result.get("source", ""),
                discovery_truncated=bool(result.get("discovery", {}).get("truncated")))


def isoform_summary(result):
    """Summarize nominated genes and RNA-covered coding isoform groups."""
    targets = {target["gene_id"] for target in result["targets"]}
    paths = result["isoforms"] + result["partial_sequences"]
    observed = [row for row in paths if int(row.get("compatible_fragments", 0)) > 0]
    genes = {gene for row in observed for gene in row.get("gene_ids", ())}
    groups = [row for row in result["protein_groups"]
              if row.get("fully_covered_isoform_ids")
              and int(row.get("compatible_fragments", 0)) > 0]
    return dict(target_genes=len(targets), rna_supported_genes=len(genes & targets),
                observed_paths=len(observed), fully_covered_protein_groups=len(groups),
                source=result.get("source", ""),
                discovery_truncated=bool(result.get("discovery", {}).get("truncated")))


def compose(sample_id, *, small_variants=None, sv=(), fusion=(), viral=(), isoforms=()):
    """Build a source-separated inventory, retaining every supplied export."""
    output = dict(schema="isovar.antigen_evidence_report.v1", sample_id=sample_id,
                  categories={name: dict(status="not_assessed", records=[])
                              for name in CATEGORIES})
    if small_variants is not None:
        rows = load_small_variants(small_variants, sample_id)
        output["categories"]["small_variants"] = dict(
            status="assessed", source=str(small_variants),
            records=[small_variant_summary(rows)])
    for category, paths, summarizer in (
            ("sv", sv, sv_summary), ("fusion", fusion, fusion_summary),
            ("viral", viral, viral_summary), ("isoforms", isoforms, isoform_summary)):
        seen = set()
        records = []
        for path in paths:
            result = load_json(path, SCHEMAS[category], sample_id)
            item = summarizer(result)
            identity = (item.get("event_id", ""), item.get("source", ""))
            if category in {"sv", "fusion"} and identity in seen:
                raise ValueError(f"Duplicate {category} event and source: {identity}")
            seen.add(identity)
            records.append(dict(source_file=str(path), **item))
        if records:
            output["categories"][category] = dict(status="assessed", records=records)
    if not any(category["status"] == "assessed"
               for category in output["categories"].values()):
        raise ValueError("Supply at least one Isovar evidence export")
    return output
