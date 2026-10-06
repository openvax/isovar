"""Compare acquired Sid SV/indel ORF hypotheses with annotated human proteins.

This is a candidate screen, not evidence of translation or tumor specificity.
The exact and I/L-collapsed comparisons cover a pinned Ensembl protein FASTA.
"""

import argparse
from collections import defaultdict
import gzip
from pathlib import Path

from Bio import SeqIO

from examples.sid_sv_audit.inventory import digest, identity, load_inventory, read_json, write_json


def proteome_screen(sequences, path, length=8):
    """Scan every annotated protein; peptide novelty and full-sequence hits differ."""
    sequences = set(sequences)
    queries = {s[i:i + length] for s in sequences for i in range(len(s) - length + 1)}
    collapse = lambda s: s.replace("I", "L")
    il_queries = {collapse(s) for s in queries}
    exact, il, full, il_full = set(), set(), defaultdict(list), defaultdict(list)
    prefixes, il_prefixes = defaultdict(set), defaultdict(set)
    for s in sequences:
        if len(s) >= length:
            prefixes[s[:length]].add(s)
            il_prefixes[collapse(s[:length])].add(s)
    proteins = bases = 0
    opener = gzip.open if Path(path).suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            proteins += 1
            protein = str(record.seq)
            bases += len(protein)
            collapsed = collapse(protein)
            hits, il_hits = set(), set()
            for i in range(len(protein) - length + 1):
                peptide, normalized = protein[i:i + length], collapsed[i:i + length]
                if peptide in queries:
                    exact.add(peptide)
                    hits.add(peptide)
                if normalized in il_queries:
                    il.add(normalized)
                    il_hits.add(normalized)
            for s in {s for h in hits for s in prefixes[h]}:
                if s in protein and len(full[s]) < 5:
                    full[s].append(record.id)
            for s in {s for h in il_hits for s in il_prefixes[h]}:
                if collapse(s) in collapsed and len(il_full[s]) < 5:
                    il_full[s].append(record.id)
    if proteins == 0:
        raise ValueError("Empty annotated proteome")
    return dict(proteome_sha256=digest(path), proteins_scanned=proteins, amino_acids_scanned=bases,
        peptide_length=length, sequences={s: dict(
            novel_exact_peptides=sorted({s[i:i + length] for i in range(len(s) - length + 1)} - exact),
            novel_il_collapsed_peptides=sorted({s[i:i + length] for i in range(len(s) - length + 1)
                                               if collapse(s[i:i + length]) not in il}),
            exact_full_sequence_hits=full[s], il_collapsed_full_sequence_hits=il_full[s]) for s in sorted(sequences)})


def collect(audits, variants):
    """Retain each source, geometry and sequence; no support is pooled."""
    rows, pins, outcomes = [], {}, []
    for directory in audits:
        manifest = load_inventory(directory)
        pins[str(directory / "inventory.json.gz")] = digest(directory / "inventory.json.gz")
        for path in sorted((directory / "results").glob("*/*.json.gz")):
            result = read_json(path)
            pins[str(path)] = digest(path)
            req = result["request"]
            gid, sid = req["geometry_id"], req["source_id"]
            names = sorted({manifest["nominations"][n]["original"].get("name", n)
                            for n in manifest["geometries"][gid]["nominations"]})
            for c in result.get("candidates", []):
                rows.append(dict(kind="SV_RNA_ORF", target=names, geometry_id=gid, source_id=sid,
                    orientation=req["orientation"], candidate_id=c["candidate_id"],
                    amino_acids=c["amino_acids"], nucleotide_sequence=c["nucleotide_sequence"],
                    sequence_group_id=identity([sid, c["nucleotide_sequence"]]),
                    ends_with_stop_codon=c["ends_with_stop_codon"], event_relations=c["event_relations"],
                    fragments=c["rna_support"]["fragments"], cells=c["rna_support"]["cells"],
                    missing_quality_reads=c["rna_support"]["missing_quality_reads"],
                    start_evidence=c["start_evidence_summary"],
                    frame_statuses=sorted({o["frame_status"] for o in c["occurrences"]}),
                    uncertainty_flags=c["uncertainty_flags"],
                    input_limitations=result.get("input_limitations", []),
                    discovery=result.get("discovery"), support_acquisition=result.get("support_acquisition"),
                    result_path=str(path)))
    for path in sorted(variants.glob("*/*.json.gz")):
        result = read_json(path)
        if not result["request"].get("convert_ucsc_contig_names"):
            raise ValueError("Unvalidated annotation chromosome labels")
        pins[str(path)] = digest(path)
        for pin in result["request"]["reference_source_pins"]:
            verified = pins.get(pin["path"])
            if (verified if verified is not None else digest(pin["path"])) != pin["sha256"]:
                raise ValueError("Annotation source changed")
            pins[pin["path"]] = pin["sha256"]
        outcomes.append(dict(target_id=result["request"]["target_id"], source_id=result["request"]["source_id"],
                             status=result["status"], proteins=len(result.get("proteins", [])),
                             allele_support=result.get("allele_support"), error=result.get("error"),
                             request=result["request"], acquisition_status=result["acquisition_status"], path=str(path)))
        for p in result.get("proteins", []):
            rows.append(dict(kind="indel_frameshift" if p["frameshift"] else "inframe_coding_window",
                target=[result["request"]["target_id"]], source_id=result["request"]["source_id"],
                amino_acids=p["amino_acids"], mutant_amino_acids=p["mutant_amino_acids"],
                mutation_interval=p["mutation_interval"], frameshift=p["frameshift"],
                ends_with_stop_codon=p["ends_with_stop_codon"],
                fragments=p["full_window_fragments"], contributing_fragments=p["contributing_fragments"],
                q20_full_window_fragments=p["q20_full_window_fragments"],
                fragments_with_unavailable_qualities=p["fragments_with_unavailable_qualities"],
                transcript_ids=p["transcript_ids"], independent_translation_checks=p["independent_translation_checks"],
                acquisition_status=result["acquisition_status"], result_path=str(path)))
    return rows, pins, outcomes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--audit", type=Path, action="append", required=True)
    parser.add_argument("--variants", type=Path, required=True)
    parser.add_argument("--proteome", type=Path, required=True)
    args = parser.parse_args()
    rows, pins, outcomes = collect(args.audit, args.variants)
    comparison = proteome_screen((r["amino_acids"] for r in rows), args.proteome)
    for row in rows:
        row["annotated_protein_comparison"] = comparison["sequences"][row["amino_acids"]]
    pins[str(args.proteome)] = comparison["proteome_sha256"]
    write_json(args.output, dict(scope="RNA_candidates_in_acquired_selected_regions", candidates=rows,
        comparison={k: v for k, v in comparison.items() if k != "sequences"}, input_pins=pins,
        small_variant_outcomes=outcomes, code_sha256=digest(__file__),
        translation_observed=False, tumor_specificity_assessed=False,
        all_SVs_or_neoorfs_searched=False))
    print("Screened", len(rows), "source-specific hypotheses against", comparison["proteins_scanned"], "proteins")


if __name__ == "__main__":
    main()
