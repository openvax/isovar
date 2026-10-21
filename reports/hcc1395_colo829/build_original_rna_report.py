"""Rebuild source-separated Isovar calls from bounded original RNA BAMs.

The BAMs contain original SRA reads realigned to GRCh38. TAGET polished
representatives are included only as a separate sequence comparison.
"""

import argparse
import csv
import gzip
from pathlib import Path

import pysam
from pyensembl import cached_release
from varcode import Variant

from isovar import ProteinSequenceCreator, ReadCollector


HERE = Path(__file__).parent
DEFAULT_BAMS = Path(__file__).parents[2] / "tests/data/hcc1395_colo829_original"
ASSAYS = (
    ("HCC1395", "SRR23456709", "ONT_direct_RNA", "hcc1395_srr23456709.bam"),
    ("HCC1395", "SRR23456710", "ONT_direct_cDNA", "hcc1395_srr23456710.bam"),
    ("COLO829", "SRR8615617", "Illumina_paired_RNA", "colo829_srr8615617.bam"),
)
POLISHED_FIELDS = (
    "tumor_alt_isoforms", "tumor_ref_isoforms", "normal_alt_isoforms",
    "normal_ref_isoforms", "tumor_alt_fl_count", "tumor_ref_fl_count",
    "normal_alt_fl_count", "normal_ref_fl_count",
)
FIELDS = (
    "sample", "gene", "variant_id", "dna_source", "rna_run", "rna_assay",
    "annotation", "ref_reads", "alt_reads", "other_reads", "alt_source_reads",
    "protein_status", "protein_sequence", "mutation_start_idx",
    "mutation_end_idx", "protein_supporting_fragments", *POLISHED_FIELDS,
)


def polished_rows(path):
    """Index the existing TAGET representative-sequence summary."""
    with gzip.open(path, "rt", newline="") as source:
        return {row["variant_id"]: row for row in csv.DictReader(source, delimiter="\t")}


def report_rows(candidates, bam_dir, polished):
    """Yield per-source calls and their original alternate-read witnesses."""
    genome = cached_release(87)
    collector = ReadCollector()
    creator = ProteinSequenceCreator()
    for sample, accession, assay, bam_name in ASSAYS:
        with pysam.AlignmentFile(str(bam_dir / bam_name), "rb") as alignment:
            for candidate in candidates:
                if candidate["sample"] != sample:
                    continue
                variant = Variant(candidate["chrom"].removeprefix("chr"),
                                  int(candidate["pos"]), candidate["ref"],
                                  candidate["alt"], genome)
                variant_id = (f"{candidate['chrom']}:{candidate['pos']}:"
                              f"{candidate['ref']}>{candidate['alt']}")
                evidence = collector.read_evidence_for_variant(variant, alignment)
                proteins = creator.sorted_protein_sequences_for_variant(variant, evidence)
                protein = proteins[0] if proteins else None
                row = {
                    "sample": sample, "gene": candidate["gene"],
                    "variant_id": variant_id, "dna_source": candidate["dna_source"],
                    "rna_run": accession, "rna_assay": assay,
                    "annotation": variant.effects().top_priority_effect().short_description,
                    "ref_reads": len(evidence.ref_reads),
                    "alt_reads": len(evidence.alt_reads),
                    "other_reads": len(evidence.other_reads),
                    "alt_source_reads": sum(read.source_read_count for read in evidence.alt_reads),
                    "protein_status": "translated" if protein else "no_translation",
                    "protein_sequence": protein.amino_acids if protein else "NA",
                    "mutation_start_idx": protein.mutation_start_idx if protein else "NA",
                    "mutation_end_idx": protein.mutation_end_idx if protein else "NA",
                    "protein_supporting_fragments": protein.num_supporting_fragments if protein else "NA",
                }
                row.update({field: polished.get(variant_id, {}).get(field, "NA")
                            for field in POLISHED_FIELDS})
                witnesses = [(sample, accession, variant_id, read.name, read.source_read_count)
                             for read in evidence.alt_reads]
                yield row, witnesses


def write_report(candidates_path, bam_dir, output_dir):
    """Write calls and witness IDs; return the number of evaluated calls."""
    output_dir.mkdir(parents=True, exist_ok=True)
    with candidates_path.open(newline="") as source:
        candidates = list(csv.DictReader(source, delimiter="\t"))
    polished = polished_rows(HERE / "colo829_polished_snv_summary.tsv.gz")
    witness_fields = ("sample", "rna_run", "variant_id", "read_id", "source_read_count")
    count = 0
    with (output_dir / "original_rna_protein_calls.tsv").open("w", newline="") as calls, \
            (output_dir / "original_rna_alt_witnesses.tsv").open("w", newline="") as reads:
        call_writer = csv.DictWriter(calls, fieldnames=FIELDS, delimiter="\t",
                                     lineterminator="\n")
        read_writer = csv.writer(reads, delimiter="\t", lineterminator="\n")
        call_writer.writeheader()
        read_writer.writerow(witness_fields)
        for row, witnesses in report_rows(candidates, bam_dir, polished):
            call_writer.writerow(row)
            read_writer.writerows(witnesses)
            count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, default=HERE / "original_rna_candidates.tsv")
    parser.add_argument("--bam-dir", type=Path, default=DEFAULT_BAMS)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    print(f"Evaluated {write_report(args.candidates, args.bam_dir, args.output_dir)} source calls")


if __name__ == "__main__":
    main()
