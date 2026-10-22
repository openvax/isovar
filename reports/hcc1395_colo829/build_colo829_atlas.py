"""Rebuild the COLO829 all-small-variant original-RNA atlas from an indexed BAM.

The BAM must come from the complete SRR8615617 run aligned to GRCh38 by
``align_full_colo829_rna.sh``. Protein sequences are local Isovar translations,
not reconstructed full-length mutant isoforms.
"""

import argparse
import csv
import gzip
import logging
from pathlib import Path

import pysam
from pyensembl import cached_release
from varcode import Variant

from isovar import ProteinSequenceCreator, ReadCollector
from isovar.evidence_report import direct_read_context_proteins
from isovar.read_identity import fragment_ids


HERE = Path(__file__).parent
DEFAULT_BAM = HERE.parents[1] / "tests/data/colo829_whole_run/SRR8615617.all_sites.grch38.bam"
MAX_ALT_READS_FOR_OVERLAP_ASSEMBLY = 750
CALL_FIELDS = (
    "sample", "dna_source", "class", "call_id", "chrom", "pos", "ref", "alt",
    "dna_filter", "rna_run", "raw_alignment_overlap", "ref_reads",
    "alt_reads", "other_reads", "ref_fragments", "alt_fragments",
    "other_fragments", "alt_source_reads", "chimeric_path_warning_reads",
    "rna_status",
    "genomic_overlapping_transcript_ids", "assembly_status", "protein_status",
    "protein_count",
)
PROTEIN_FIELDS = (
    "dna_source", "class", "call_id", "rna_run", "rank", "protein_sequence",
    "mutation_start_idx", "mutation_end_idx", "supporting_fragments",
    "supporting_reads", "compatible_transcript_ids", "sequence_scope",
)
WITNESS_FIELDS = ("dna_source", "class", "call_id", "rna_run", "read_id",
                  "source_read_count")


class AssemblyFallbackTracker(logging.Handler):
    """Note when Isovar bypasses greedy merge at a complex allele locus."""

    def __init__(self):
        super().__init__(level=logging.WARNING)
        self.skipped = False

    def emit(self, record):
        if record.getMessage().startswith("Too many variant sequences"):
            self.skipped = True


class ChimericPathWarningTracker(logging.Handler):
    """Count distinct reads whose supplementary path failed validation."""

    def __init__(self):
        super().__init__(level=logging.WARNING)
        self.read_ids = set()

    def emit(self, record):
        if record.getMessage().startswith("Ignoring chimeric-path evidence"):
            self.read_ids.add(record.args[0])


def source_calls(inventory, selected_ids=None, sample="COLO829",
                 dna_source="SMaHT_COLO829BLT50"):
    """Iterate exact small-variant source records for one call set."""
    with gzip.open(inventory, "rt", newline="") as source:
        for row in csv.DictReader(source, delimiter="\t"):
            if row["sample"] != sample or row["source"] != dna_source:
                continue
            if row["class"] not in ("SNV", "indel"):
                continue
            if selected_ids is None or row["call_id"] in selected_ids:
                yield row


def evaluate_call(source, alignment, genome, collector, creator, rna_run="SRR8615617"):
    """Return one call, zero or more local proteins, and ALT read witnesses."""
    chrom = source["chrom1"]
    pos = int(source["pos1"])
    variant = Variant(chrom.removeprefix("chr"), pos, source["ref"],
                      source["alt"], genome)
    path_tracker = ChimericPathWarningTracker()
    path_logger = logging.getLogger("isovar.chimeric_alignment")
    path_logger.addHandler(path_tracker)
    try:
        evidence = collector.read_evidence_for_variant(variant, alignment)
    finally:
        path_logger.removeHandler(path_tracker)
    ref, alt, other = evidence.ref_reads, evidence.alt_reads, evidence.other_reads
    if ref or alt or other:
        raw_overlap = "NA"
    else:
        # A spliced or filtered alignment can overlap the genomic locus without
        # yielding an allele call. Keep this distinct from genuine no coverage.
        raw_overlap = alignment.count(chrom, pos - 1,
                                      pos - 1 + max(1, len(source["ref"])))
    if alt:
        rna_status = "alternate_observed"
    elif ref and other:
        rna_status = "reference_and_other_only"
    elif ref:
        rna_status = "reference_only"
    elif other:
        rna_status = "other_allele_only"
    elif raw_overlap:
        rna_status = "alignment_without_allele_call"
    else:
        rna_status = "no_alignment_overlap"

    tracker = AssemblyFallbackTracker()
    direct_context_only = len(alt) > MAX_ALT_READS_FOR_OVERLAP_ASSEMBLY
    if direct_context_only:
        proteins = direct_read_context_proteins(variant, alt, creator)
    elif alt:
        assembly_logger = logging.getLogger("isovar.assembly")
        assembly_logger.addHandler(tracker)
        try:
            proteins = creator.sorted_protein_sequences_for_variant(variant, evidence)
        finally:
            assembly_logger.removeHandler(tracker)
    else:
        proteins = []
    protein_status = ("multiple_local_contexts" if len(proteins) > 1 else
                      "translated_local_context" if proteins else
                      "no_translation" if alt else "not_attempted_no_alt")
    overlapping_ids = genome.transcript_ids_at_locus(variant.contig, pos)
    call = {
        "sample": source["sample"], "dna_source": source["source"],
        "class": source["class"], "call_id": source["call_id"],
        "chrom": chrom, "pos": pos, "ref": source["ref"], "alt": source["alt"],
        "dna_filter": source["filter"], "rna_run": rna_run,
        "raw_alignment_overlap": raw_overlap,
        "ref_reads": len(ref), "alt_reads": len(alt), "other_reads": len(other),
        "ref_fragments": len(fragment_ids(ref)),
        "alt_fragments": len(fragment_ids(alt)),
        "other_fragments": len(fragment_ids(other)),
        "alt_source_reads": sum(read.source_read_count for read in alt),
        "chimeric_path_warning_reads": len(path_tracker.read_ids),
        "rna_status": rna_status,
        "genomic_overlapping_transcript_ids": ";".join(sorted(overlapping_ids)),
        "assembly_status": ("direct_read_contexts_high_depth" if direct_context_only else
                            "greedy_merge_skipped_complexity" if tracker.skipped else
                            "default" if alt else "not_attempted_no_alt"),
        "protein_status": protein_status, "protein_count": len(proteins),
    }
    protein_rows = [{
        "dna_source": source["source"], "class": source["class"],
        "call_id": source["call_id"], "rna_run": rna_run, "rank": rank,
        "protein_sequence": protein.amino_acids,
        "mutation_start_idx": protein.mutation_start_idx,
        "mutation_end_idx": protein.mutation_end_idx,
        "supporting_fragments": protein.num_supporting_fragments,
        "supporting_reads": protein.num_supporting_reads,
        "compatible_transcript_ids": ";".join(sorted(protein.transcript_ids)),
        "sequence_scope": "local_context",
    } for rank, protein in enumerate(proteins, 1)]
    witnesses = [{
        "dna_source": source["source"], "class": source["class"],
        "call_id": source["call_id"], "rna_run": rna_run,
        "read_id": read.name, "source_read_count": read.source_read_count,
    } for read in alt]
    return call, protein_rows, witnesses


def write_atlas(inventory, bam, output_dir, selected_ids=None,
                sample="COLO829", dna_source="SMaHT_COLO829BLT50",
                rna_run="SRR8615617", prefix="colo829"):
    """Write three gzip TSVs, returning their row counts."""
    output_dir.mkdir(parents=True, exist_ok=True)
    genome = cached_release(87)
    collector = ReadCollector()
    creator = ProteinSequenceCreator(max_protein_sequences_per_variant=0)
    counts = [0, 0, 0]
    with pysam.AlignmentFile(str(bam), "rb") as alignment, \
            gzip.open(output_dir / f"{prefix}_atlas_calls.tsv.gz", "wt", newline="") as calls, \
            gzip.open(output_dir / f"{prefix}_atlas_proteins.tsv.gz", "wt", newline="") as proteins, \
            gzip.open(output_dir / f"{prefix}_atlas_alt_witnesses.tsv.gz", "wt", newline="") as reads:
        call_writer = csv.DictWriter(calls, CALL_FIELDS, delimiter="\t", lineterminator="\n")
        protein_writer = csv.DictWriter(proteins, PROTEIN_FIELDS, delimiter="\t",
                                        lineterminator="\n")
        witness_writer = csv.DictWriter(reads, WITNESS_FIELDS, delimiter="\t",
                                        lineterminator="\n")
        for writer in (call_writer, protein_writer, witness_writer):
            writer.writeheader()
        for source in source_calls(inventory, selected_ids, sample, dna_source):
            call, protein_rows, witnesses = evaluate_call(
                source, alignment, genome, collector, creator, rna_run)
            call_writer.writerow(call)
            protein_writer.writerows(protein_rows)
            witness_writer.writerows(witnesses)
            counts[0] += 1
            counts[1] += len(protein_rows)
            counts[2] += len(witnesses)
    return tuple(counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, default=HERE / "somatic_call_inventory.tsv.gz")
    parser.add_argument("--bam", type=Path, default=DEFAULT_BAM)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    print("calls, proteins, ALT witnesses:",
          write_atlas(args.inventory, args.bam, args.output_dir))


if __name__ == "__main__":
    main()
