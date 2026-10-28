"""Recount pinned atlas exemplars from original RNA using the library export.

Run ``python -m reports.hcc1395_colo829.recount_protein_evidence``. The compact
receipt keeps one displayed nucleotide path, its witness sets and comparisons
across all retained models. Full native exports are written beside the receipt.
The historical atlas defines the target protein; it does not define its counts.
"""

import argparse
import csv
import gzip
import json
from pathlib import Path

import pysam
from pyensembl import EnsemblRelease
from varcode import Variant

from isovar import ReadCollector, ProteinSequenceCreator, IsovarResult, export_protein_hypotheses
from isovar.evidence_report import direct_read_context_proteins
from isovar.protein_reporting import best_supported_translation, protein_report_record
from isovar.sv_rna_audit import file_sha256


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "reports/hcc1395_colo829"
CASES = (
    ("COLO829", "CDKN2A", "chr9", 21971154, "CCG", "C"),
    ("COLO829", "FHL2", "chr2", 105361307, "G", "GCCCC"),
    ("COLO829", "GOLGA5", "chr14", 92811641, "GAGA", "G"),
    ("COLO829", "TUBB2A", "chr6", 3153911, "G", "A"),
    ("HCC1395", "C8orf76", "chr8", 123220221, "A", "AGGGCTAC"),
)


def table(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def input_paths(sample: str) -> tuple[Path, Path, Path]:
    prefix = "colo829" if sample == "COLO829" else "hcc1395_direct_rna"
    bam = ("colo829_whole_run/SRR8615617.all_sites.grch38.bam" if sample == "COLO829" else
           "hcc1395_whole_run/SRR23456709.fixture.grch38.bam")
    return BASE / (prefix + "_atlas_calls.tsv.gz"), BASE / (prefix + "_atlas_proteins.tsv.gz"), ROOT / "tests/data" / bam


def compact_receipt(export: dict, sample: str, gene: str, call_id: str,
                    inputs: dict, settings: dict) -> dict:
    """Keep display selection separate from all-model classification and support."""
    protein, = export["events"][0]["protein_hypotheses"]
    selected = best_supported_translation(protein)
    coverage = selected["coding_support"]
    supports = [protein["rna_support"], selected["rna_support"],
                *[coverage[key] for key in ("allele_observation_span_support",
                                           "single_read_span_support", "paired_template_span_support")]]
    ids = {support["evidence_set_id"] for support in supports}
    return dict(sample=sample, gene=gene, inputs=inputs, settings=settings,
        metrics=protein_report_record(protein, call_id), displayed_translation=selected,
        contributor_support=protein["rna_support"],
        model_comparisons=[dict(translation_id=t["translation_id"],
            reference_comparison=t["reference_comparison"]) for t in protein["translations"]],
        evidence_sets={key: export["evidence_sets"][key] for key in sorted(ids)})


def recount(case: tuple, genome: EnsemblRelease, creator: ProteinSequenceCreator) -> tuple[dict, dict]:
    sample, gene, chrom, position, ref, alt = case
    calls_path, proteins_path, bam_path = input_paths(sample)
    call, = [row for row in table(calls_path) if (row["chrom"], int(row["pos"]), row["ref"], row["alt"]) ==
             (chrom, position, ref, alt)]
    target, = [row for row in table(proteins_path) if row["call_id"] == call["call_id"] and row["rank"] == "1"]
    variant = Variant(chrom, position, ref, alt, genome, convert_ucsc_contig_names=True)
    with pysam.AlignmentFile(bam_path) as alignment:
        evidence = ReadCollector().read_evidence_for_variant(variant, alignment)
        direct = len(evidence.alt_reads) > 750
        proteins = (direct_read_context_proteins(variant, evidence.alt_reads, creator) if direct else
                    creator.sorted_protein_sequences_for_variant(variant, evidence))
        selected = [protein for protein in proteins if protein.amino_acids == target["protein_sequence"] and
                    (protein.mutation_start_idx, protein.mutation_end_idx) ==
                    (int(target["mutation_start_idx"]), int(target["mutation_end_idx"]))]
        if len(selected) != 1:
            raise ValueError(f"Expected one exact historical target for {gene}; found {len(selected)}")
        result = IsovarResult(variant=variant, read_evidence=evidence, predicted_effect=None,
                             sorted_protein_sequences=selected, protein_sequence_settings=creator.settings())
        export = export_protein_hypotheses([result], sample_id=sample, source=call["rna_run"],
                                           alignment_header=alignment.header)
    inputs = {str(path.relative_to(ROOT)): file_sha256(path) for path in
              (calls_path, proteins_path, bam_path, Path(str(bam_path) + ".bai"))}
    settings = dict(reference="GRCh38", annotation="Ensembl 87", library=call["rna_run"],
                    read_collection="ReadCollector defaults", protein_sequence=creator.settings(),
                    context_method="direct_read_context_proteins" if direct else "assembled_variant_sequences")
    return compact_receipt(export, sample, gene, call["call_id"], inputs, settings), export


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=BASE / "focused_proteins/protein_evidence_recount.json.gz")
    parser.add_argument("--native-directory", type=Path, default=ROOT / "output/protein-evidence-recount")
    args = parser.parse_args()
    genome = EnsemblRelease(87)
    creator = ProteinSequenceCreator(max_protein_sequences_per_variant=0)
    receipts = []
    args.native_directory.mkdir(parents=True, exist_ok=True)
    for case in CASES:
        receipt, export = recount(case, genome, creator)
        native = args.native_directory / (receipt["gene"] + ".json.gz")
        with gzip.open(native, "wt") as stream:
            json.dump(export, stream, sort_keys=True)
        receipts.append(receipt)
        print(json.dumps(dict(gene=receipt["gene"], **receipt["metrics"])), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.output, "wt") as stream:
        json.dump(dict(schema="isovar.exemplar_protein_recount.v1", isovar_version=export["isovar_version"],
                       cases=receipts), stream, sort_keys=True)


if __name__ == "__main__":
    main()
