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
from .rna_evidence import input_sample_id, rna_support
from .sv_rna_relations import ADJACENCY_JOIN_RELATIONS


SCHEMAS = {
    "sv": "isovar.sv_rna_candidates.v5",
    "fusion": "isovar.fusion_rna.v3",
    "viral": "isovar.viral_sequences.v1",
    "isoforms": "isovar.antigen_isoforms.v1",
}
CATEGORIES = ("small_variants", "sv", "fusion", "viral", "isoforms")


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
        try:
            if not needed.issubset(reader.fieldnames or ()):
                raise ValueError(f"{path}: not an Isovar evidence-report per-call TSV")
            rows = list(reader)
        except csv.Error as error:
            raise ValueError(f"{path}: not an Isovar evidence-report per-call TSV ({error})") from None
    if any(row["sample_id"] != sample_id for row in rows):
        raise ValueError(f"{path}: sample_id does not match {sample_id}")
    if len({row["call_id"] for row in rows}) != len(rows):
        raise ValueError(f"{path}: duplicate small-variant call_id")
    return rows



def _measurement(name, subject, support=None, *, status=None, scope=None, **metadata):
    support = support or {}
    return dict(name=name, subject=subject, reads=support.get("reads"), fragments=support.get("fragments"),
                status=status, scope=scope, **metadata)


def _path_measurements(paths, names, scope):
    """Project native support and original witness IDs through the public counter."""
    rows = []
    for path in paths:
        identifier = path["id"]
        name = names(path) + " / " + identifier
        status = path.get("sequence_status", path.get("coding_status"))
        support = path.get("rna_support", dict(reads=path.get("supporting_segments"),
                                               fragments=path.get("supporting_fragments")))
        rows.append(_measurement(name, "RNA path contributors", support, status=status, scope=scope, path_id=identifier))
        for key, subject in (("compatible_fragments", "Compatible with RNA path"),
                             ("unique_fragments", "Unique among retained paths")):
            rows.append(_measurement(name, subject, dict(fragments=path.get(key)), status=status, scope=scope, path_id=identifier))
        witnesses = path.get("whole_path_segments")
        rows.append(_measurement(name, "Whole retained RNA path",
            None if witnesses is None else rna_support(witnesses), status=status, scope=scope, path_id=identifier))
        for index, coding in enumerate(path.get("coding", ())):
            witnesses = coding.get("whole_cds_segments") if coding.get("status") == "complete_CDS" else None
            rows.append(_measurement(name, "Whole complete CDS %d" % (index + 1),
                None if witnesses is None else rna_support(witnesses), status=coding.get("status"), scope=scope,
                path_id=identifier, coding_index=index))
    return rows


def _junction_measurements(result, scope):
    rows = []
    for index, path in enumerate(result["paths"]):
        identifier = path.get("path_id", str(index + 1))
        name = result["event_id"] + " / " + identifier
        for junction_index, junction in enumerate(path.get("junctions", ())):
            rows.append(_measurement(name, "Direct RNA junction %d" % (junction_index + 1),
                junction.get("direct_support"), status=junction.get("relation"), scope=scope,
                path_id=identifier, junction_index=junction_index))
        if "full_interval_support" in path:
            rows.append(_measurement(name, "Whole retained RNA path", path["full_interval_support"],
                status=path.get("event_linkage", {}).get("status"), scope=scope, path_id=identifier))
    return rows

def small_variant_summary(rows):
    """Summarize literal DNA calls and original RNA ALT support."""
    exonic = [row for row in rows if row["exon_overlap"] == "True"]
    ncorf_values = [row["orfdata_coding_block_overlap"] for row in rows]
    ncorf_assessed = bool(rows) and all(value in {"True", "False"}
                                        for value in ncorf_values)
    # Blank counts are unmeasured (e.g. a chromosome absent from the reads), not zero.
    measured = [row for row in rows if row["alt_reads"] != ""]
    by_class = {}
    for row in exonic:
        if row["alt_reads"] != "":
            by_class.setdefault(variant_type(row), []).append(int(row["alt_reads"]))
    return dict(input_records=len(rows), exon_involved=len(exonic),
                rna_unmeasured=len(rows) - len(measured),
                rna_alt_observed=sum(int(row["alt_reads"]) > 0 for row in measured),
                local_products=sum(int(row["protein_count"] or 0) > 0 for row in measured),
                ncorf_reference_candidates=(sum(value == "True" for value in ncorf_values)
                                            if ncorf_assessed else None),
                exonic_reads_by_class=by_class, calls=rows)


def _direct_fragments(path, *, relations=None):
    return max((int(junction.get("direct_support", {}).get("fragments", 0))
                for junction in path.get("junctions", ())
                if relations is None or junction.get("relation") in relations), default=0)


def sv_summary(result):
    """Summarize one nominated DNA SV without asserting tumor protein expression."""
    paths = result["paths"]
    linked = [path for path in paths if path.get("event_linkage", {}).get("status")
              in ADJACENCY_JOIN_RELATIONS]
    supported = [path for path in linked
                 if _direct_fragments(path, relations=ADJACENCY_JOIN_RELATIONS) > 0]
    translated = [path for path in supported
                  if any(any(relation in ADJACENCY_JOIN_RELATIONS for relation in row.get("departure_relations", ()))
                         for row in path.get("translations", ()))]
    exploratory = sum(len(path.get("exploratory_orfs", {}).get("candidates", ()))
                      for path in paths)
    return dict(event_id=result["event_id"], source=result.get("source", ""),
                input_events=1, candidate_paths=len(paths),
                linked_paths=len(linked),
                max_direct_junction_fragments=max(
                    (_direct_fragments(path, relations=ADJACENCY_JOIN_RELATIONS) for path in linked),
                    default=0),
                rna_linked_events=bool(supported),
                event_linked_translated_paths=len(translated),
                exploratory_orf_candidates=exploratory,
                status=result["status"],
                measurements=_junction_measurements(result,
                    result.get("support_acquisition", {}).get("scope", "collected_regional_records")))


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
                status=result["status"], reasons=result.get("reasons", []),
                measurements=_junction_measurements(result, "supplied_fusion_observations"))


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
                        "sequence_assigned_fragments", "status")} | {"aligned_segments": row.get("aligned_segments")}
                         for row in contigs],
                source=result.get("source", ""),
                discovery_truncated=bool(result.get("discovery", {}).get("truncated")),
                measurements=_path_measurements(sequences, lambda row: ", ".join(row.get("viral_contigs", ())) or "Viral path",
                    result.get("support", {}).get("scope", "scope_not_recorded")))


def isoform_summary(result):
    """Summarize nominated genes and RNA-covered coding isoform groups."""
    gene_names = {target["gene_id"]: target.get("gene_name", target["gene_id"]) for target in result["targets"]}
    targets = set(gene_names)
    paths = result["isoforms"] + result["partial_sequences"]
    observed = [row for row in paths if int(row.get("compatible_fragments", 0)) > 0]
    genes = {gene for row in observed for gene in row.get("gene_ids", ())}
    groups = [row for row in result["protein_groups"]
              if row.get("fully_covered_isoform_ids")
              and int(row.get("compatible_fragments", 0)) > 0]
    return dict(target_genes=len(targets), rna_supported_genes=len(genes & targets),
                observed_paths=len(observed), fully_covered_protein_groups=len(groups),
                source=result.get("source", ""),
                discovery_truncated=bool(result.get("discovery", {}).get("truncated")),
                targets=result["targets"],
                measurements=_path_measurements(paths,
                    lambda row: ", ".join(gene_names.get(gene, gene) for gene in row.get("gene_ids", ())) or "Unassigned gene",
                    result.get("support", {}).get("scope", "scope_not_recorded")))



def compose(sample_id, *, small_variants=None, sv=(), fusion=(), viral=(), isoforms=(), input_sample_ids=None):
    """Build a source-separated inventory with explicit input-to-report mapping.

    input_sample_ids maps input file paths to their expected native sample IDs.
    Unmapped inputs must match sample_id. Every export keeps its original sample,
    source and measurements; the report label never pools their RNA counts.
    small_variants accepts either one existing TSV path or an iterable of paths.
    """
    small = [] if small_variants is None else ([small_variants] if isinstance(small_variants, (str, Path)) else list(small_variants))
    inputs = dict(small_variants=small, sv=list(sv), fusion=list(fusion), viral=list(viral), isoforms=list(isoforms))
    paths = [Path(path).resolve() for paths in inputs.values() for path in paths]
    if len(set(paths)) != len(paths):
        raise ValueError("Duplicate input export path")
    mapping = {Path(path).resolve(): value for path, value in (input_sample_ids or {}).items()}
    if mapping.keys() - set(paths):
        raise ValueError("Input sample mapping names an unsupplied export: " + ", ".join(map(str, sorted(mapping.keys() - set(paths)))))
    input_sample_id(dict(sample_id=sample_id))
    def expected(path):
        return input_sample_id(dict(sample_id=sample_id, **(
            {"input_sample_id": mapping[Path(path).resolve()]} if Path(path).resolve() in mapping else {})))

    output = dict(schema="isovar.antigen_evidence_report.v1", sample_id=sample_id,
                  categories={name: dict(status="not_assessed", records=[]) for name in CATEGORIES})
    if small:
        records = []
        for path in small:
            original_sample = expected(path)
            rows = load_small_variants(path, original_sample)
            records.append(dict(source_file=str(path), input_sample_id=original_sample, **small_variant_summary(rows)))
        output["categories"]["small_variants"] = dict(status="assessed", records=records)
        if len(small) == 1:
            output["categories"]["small_variants"]["source"] = str(small[0])
    for category, summarizer in (("sv", sv_summary), ("fusion", fusion_summary),
                                 ("viral", viral_summary), ("isoforms", isoform_summary)):
        seen, records = set(), []
        for path in inputs[category]:
            original_sample = expected(path)
            result = load_json(path, SCHEMAS[category], original_sample)
            item = summarizer(result)
            identity = (original_sample, item.get("event_id", ""), item.get("source", ""))
            if category in {"sv", "fusion"} and identity in seen:
                raise ValueError(f"Duplicate {category} event and source: {identity}")
            seen.add(identity)
            records.append(dict(source_file=str(path), input_sample_id=original_sample, **item))
        if records:
            output["categories"][category] = dict(status="assessed", records=records)
    if not any(category["status"] == "assessed" for category in output["categories"].values()):
        raise ValueError("Supply at least one Isovar evidence export")
    return output
