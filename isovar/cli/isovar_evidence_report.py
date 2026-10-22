"""Build a visual evidence report from a variant VCF and indexed tumor RNA BAM.

The PDF and per-call TSV distinguish original-read allele support, reference
annotation and local mutant protein contexts. Optional ncORF and reference-self
inputs add candidate coding-block and exact 9-mer comparisons.
"""

import csv
import hashlib
from pathlib import Path

from ..evidence_report import (
    catalog_coding_index, changed_reference_interval, direct_read_context_proteins,
    overlapping_catalog_orfs, overlaps_exon, reference_ninemer_matches, self_novelty,
    summarize,
)
from ..read_identity import fragment_ids
from .commands import parser_for_program
from .protein_sequence_args import make_protein_sequences_arg_parser, protein_sequence_creator_from_args
from .rna_args import alignment_file_from_args, read_collector_from_args
from .validation import CommandInputError, variant_collection_from_args


parser = make_protein_sequences_arg_parser(description=__doc__)
parser.set_defaults(max_protein_sequences_per_variant=0)
parser.add_argument("--sample-id", required=True, help="Identity of the RNA sample")
parser.add_argument("--output", type=Path, required=True, help="Output PDF")
parser.add_argument("--per-call", type=Path, help="Per-call TSV (default: PDF stem + .tsv)")
parser.add_argument("--dna-source", help="DNA call-set identity (default: VCF path)")
parser.add_argument("--rna-source", help="RNA read-set identity (default: BAM path)")
parser.add_argument("--nc-orf-catalog", type=Path,
                    help="Optional orfdata catalogue for reference coding-block candidates")
parser.add_argument("--self-fasta", type=Path,
                    help="Optional reference protein FASTA for exact mutant-overlapping 9-mer screen")

FIELDS = (
    "sample_id", "dna_source", "rna_source", "call_id", "class", "chrom", "pos",
    "ref", "alt", "exon_overlap", "reference_coding_effect",
    "orfdata_coding_block_overlap", "top_effect_type", "rna_status", "ref_reads",
    "alt_reads", "other_reads", "alt_fragments", "protein_count",
    "top_changed_aa_positions", "top_local_protein", "novel_9mers", "total_9mers",
    "self_fasta_sha256",
)


def _rna_status(ref, alt, other, raw_overlap):
    if alt:
        return "alternate_observed"
    if ref and other:
        return "reference_and_other_only"
    if ref:
        return "reference_only"
    if other:
        return "other_allele_only"
    return "alignment_without_allele_call" if raw_overlap else "no_alignment_overlap"


def _raw_overlap(alignment, variant):
    start = variant.original_start - 1
    end = start + max(1, len(variant.original_ref))
    for contig in (variant.original_contig,
                   variant.original_contig.removeprefix("chr"),
                   "chr" + variant.original_contig.removeprefix("chr")):
        if contig in alignment.references:
            return alignment.count(contig, start, end)
    raise ValueError(f"Variant contig {variant.original_contig} is absent from RNA alignment")


def collect(variants, alignment, collector, creator, *, catalog=None):
    """Evaluate every loaded small variant, retaining zero-support records."""
    calls, annotations, proteins = [], {}, []
    index = catalog_coding_index(catalog.records) if catalog is not None else None
    genome = variants[0].genome
    if catalog is not None and (catalog.metadata.get("genome_build") != genome.reference_name
                                or catalog.metadata.get("species") != genome.species.latin_name):
        raise ValueError("ncORF catalogue and VCF genome have incompatible references")
    for variant in variants:
        if variant.genome is not genome:
            raise ValueError("All variants must use the same reference genome")
        ident = variant.short_description
        if ident in annotations:
            raise ValueError(f"Duplicate variant identity: {ident}")
        evidence = collector.read_evidence_for_variant(variant, alignment)
        ref, alt, other = (list(getattr(evidence, field))
                           for field in ("ref_reads", "alt_reads", "other_reads"))
        raw_overlap = _raw_overlap(alignment, variant) if not (ref or alt or other) else 0
        if len(alt) > 750:
            products = direct_read_context_proteins(variant, alt, creator)
        elif alt:
            products = creator.sorted_protein_sequences_for_variant(variant, evidence)
        else:
            products = []
        pos = variant.original_start
        ref_base, alt_base = variant.original_ref, variant.original_alt
        cls = ("SNV" if len(ref_base) == len(alt_base) == 1 else
               "indel" if len(ref_base) != len(alt_base) else "replacement")
        call = dict(call_id=ident, **{"class": cls}, chrom=variant.original_contig,
                    pos=str(pos), ref=ref_base, alt=alt_base,
                    rna_status=_rna_status(ref, alt, other, raw_overlap),
                    ref_reads=str(len(ref)), alt_reads=str(len(alt)),
                    other_reads=str(len(other)),
                    alt_fragments=str(len(fragment_ids(alt))),
                    protein_count=str(len(products)))
        calls.append(call)
        effects = variant.effects()
        top_effect = effects.top_priority_effect()
        begin, end, _ = changed_reference_interval(pos, ref_base, alt_base)
        ids = (overlapping_catalog_orfs(index, variant.contig, begin - 1, end)
               if index is not None else ())
        annotations[ident] = dict(
            exon_overlap=overlaps_exon(genome, variant.contig, pos, ref_base, alt_base),
            reference_coding_effect=any(e.modifies_coding_sequence for e in effects),
            orfdata_coding_block_overlap=bool(ids) if catalog is not None else None,
            top_effect_type=type(top_effect).__name__)
        for rank, protein in enumerate(products, 1):
            proteins.append(dict(call_id=ident, rank=str(rank),
                                 protein_sequence=protein.amino_acids,
                                 mutation_start_idx=str(protein.mutation_start_idx),
                                 mutation_end_idx=str(protein.mutation_end_idx),
                                 supporting_reads=str(protein.num_supporting_reads)))
    return calls, annotations, proteins


def write_per_call(path, calls, annotations, summary, sources, self_scores):
    """Write the exact records underlying the report figures."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for call in calls:
            ident = call["call_id"]
            top = summary["top"].get(ident)
            score = self_scores.get(ident)
            writer.writerow(dict(
                sample_id=sources["sample_id"], dna_source=sources["dna_source"],
                rna_source=sources["rna_source"], **call,
                exon_overlap=annotations[ident]["exon_overlap"],
                reference_coding_effect=annotations[ident]["reference_coding_effect"],
                orfdata_coding_block_overlap=annotations[ident]["orfdata_coding_block_overlap"],
                top_effect_type=annotations[ident]["top_effect_type"],
                top_changed_aa_positions=(int(top["mutation_end_idx"]) -
                                          int(top["mutation_start_idx"])) if top else "",
                top_local_protein=top["protein_sequence"] if top else "",
                novel_9mers=score[0] if score else "",
                total_9mers=score[1] if score else "",
                self_fasta_sha256=sources.get("self_fasta_sha256", "")))


def run(args=None, *, prog=None):
    parsed = parser_for_program(parser, prog).parse_args(args)
    if parsed.output.suffix.lower() != ".pdf":
        parser.error("--output must end in .pdf")
    try:
        variants = tuple(variant_collection_from_args(parsed))
        if not variants:
            raise CommandInputError("VCF contains no loaded variants")
        collector = read_collector_from_args(parsed)
        creator = protein_sequence_creator_from_args(parsed)
        alignment = alignment_file_from_args(parsed)
        catalog = None
        if parsed.nc_orf_catalog is not None:
            import orfdata
            catalog = orfdata.load_catalog(path=parsed.nc_orf_catalog)
        with alignment:
            calls, annotations, proteins = collect(
                variants, alignment, collector, creator, catalog=catalog)
        summary = summarize(calls, annotations, proteins)
        self_scores = {}
        self_fasta_sha256 = ""
        if parsed.self_fasta is not None:
            with parsed.self_fasta.open("rb") as source:
                self_fasta_sha256 = hashlib.file_digest(source, "sha256").hexdigest()
            matches = reference_ninemer_matches(summary["top"].values(), parsed.self_fasta)
            self_scores = {key: self_novelty(value, matches)
                           for key, value in summary["top"].items()}
        sources = dict(sample_id=parsed.sample_id,
                       dna_source=parsed.dna_source or parsed.vcf,
                       rna_source=parsed.rna_source or parsed.bam,
                       self_fasta=str(parsed.self_fasta) if parsed.self_fasta else "",
                       self_fasta_sha256=self_fasta_sha256)
        parsed.output.parent.mkdir(parents=True, exist_ok=True)
        per_call = parsed.per_call or parsed.output.with_suffix(".tsv")
        per_call.parent.mkdir(parents=True, exist_ok=True)
        write_per_call(per_call, calls, annotations, summary, sources, self_scores)
        from ..evidence_pdf import write_generic_report
        write_generic_report(parsed.output, calls, annotations, summary, sources,
                             self_scores, catalog_assessed=catalog is not None)
    except (OSError, ValueError, CommandInputError) as error:
        parser_for_program(parser, prog).error(str(error))
    print(parsed.output)
    print(per_call)
