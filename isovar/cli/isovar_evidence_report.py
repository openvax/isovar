"""Build a visual evidence report from a variant VCF and indexed tumor RNA BAM.

The PDF and per-call TSV distinguish original-read allele support, reference
annotation and local mutant protein contexts. Optional ncORF and reference-self
inputs add candidate coding-block and exact 9-mer comparisons.
"""

import csv
import hashlib
import json
from pathlib import Path

from ..contig_resolution import UnresolvedContigError, require_resolved_contig
from ..effect_prediction import effect_outcomes
from ..isovar_result import IsovarResult
from ..protein_hypotheses import export_protein_hypotheses
from ..protein_reporting import protein_report_record
from ..evidence_report import (
    catalog_coding_index, changed_reference_interval,
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
parser.add_argument("--evidence-json", type=Path,
                    help="Complete per-call protein paths and witness sets (default: PDF stem + .json)")
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
    "self_fasta_sha256", "contig_resolution", "top_marked_aa_positions",
    "reference_identical_contexts", "amino_acid_divergent_contexts", "unresolved_contexts",
    "focal_changed_contexts", "top_contributor_fragments", "top_window_fragments",
    "top_single_read_fragments", "top_paired_template_fragments", "top_translation_id",
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


def _raw_overlap(alignment, variant, collector):
    start = variant.original_start - 1
    end = start + max(1, len(variant.original_ref))
    contig = require_resolved_contig(collector.resolve_variant_contig(variant, alignment), alignment)
    return alignment.count(contig, start, end)


def collect(variants, alignment, collector, creator, *, sample_id, source, catalog=None):
    """Retain assessed zero-support calls and unmeasured chromosome failures.

    Unresolved rows have null support and annotation fields. No annotation
    absence or reference allele is inferred from an unavailable RNA locus.
    """
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
        pos = variant.original_start
        ref_base, alt_base = variant.original_ref, variant.original_alt
        cls = ("SNV" if len(ref_base) == len(alt_base) == 1 else
               "indel" if len(ref_base) != len(alt_base) else "replacement")
        call = dict(call_id=ident, **{"class": cls}, chrom=variant.original_contig,
                    pos=str(pos), ref=ref_base, alt=alt_base)
        try:
            evidence = collector.read_evidence_for_variant(variant, alignment)
        except UnresolvedContigError as error:
            call.update(rna_status="unresolved_chromosome", contig_resolution=error.to_dict(),
                        ref_reads=None, alt_reads=None, other_reads=None, alt_fragments=None, protein_count=None)
            calls.append(call)
            # Annotation availability is separate from RNA availability. Avoid
            # querying an unknown annotation contig while preserving the input.
            annotations[ident] = dict(exon_overlap=None, reference_coding_effect=None,
                                      orfdata_coding_block_overlap=None, top_effect_type=None)
            continue
        ref, alt, other = (list(getattr(evidence, field))
                           for field in ("ref_reads", "alt_reads", "other_reads"))
        raw_overlap = _raw_overlap(alignment, variant, collector) if not (ref or alt or other) else 0
        assessment = creator.assess_variant_reads(
            variant, alt, assembly_mode="read_contexts" if len(alt) > 750 else "configured")
        products = creator.protein_sequences_from_assessment(assessment)
        call.update(rna_status=_rna_status(ref, alt, other, raw_overlap),
                    ref_reads=str(len(ref)), alt_reads=str(len(alt)), other_reads=str(len(other)),
                    alt_fragments=str(len(fragment_ids(alt))), protein_count=str(len(products)))
        calls.append(call)
        effects = variant.effects()
        top_effect = effects.top_priority_effect()
        result = IsovarResult(variant=variant, read_evidence=evidence, predicted_effect=top_effect,
                              sorted_protein_sequences=products, protein_sequence_settings=creator.settings(),
                              rna_assessment=assessment)
        export = export_protein_hypotheses([result], sample_id=sample_id, source=source,
                                            alignment_header=alignment.header)
        call["protein_hypotheses"] = export
        begin, end, _ = changed_reference_interval(pos, ref_base, alt_base)
        ids = (overlapping_catalog_orfs(index, variant.contig, begin - 1, end)
               if index is not None else ())
        annotations[ident] = dict(
            exon_overlap=overlaps_exon(genome, variant.contig, pos, ref_base, alt_base),
            reference_coding_effect=any(
                outcome is not None and outcome.modifies_coding_sequence is True
                for effect in effects for outcome in effect_outcomes(effect)),
            orfdata_coding_block_overlap=bool(ids) if catalog is not None else None,
            top_effect_type=type(top_effect).__name__)
        proteins.extend(protein_report_record(protein, ident)
                        for protein in export["events"][0]["protein_hypotheses"])
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
            serialized = dict({key: value for key, value in call.items() if key in FIELDS},
                              contig_resolution=json.dumps(call["contig_resolution"], sort_keys=True)
                              if call.get("contig_resolution") else "")
            context_counts = summary["context_counts_by_call"].get(ident, {})
            writer.writerow(dict(
                sample_id=sources["sample_id"], dna_source=sources["dna_source"],
                rna_source=sources["rna_source"], **serialized,
                exon_overlap=annotations[ident]["exon_overlap"],
                reference_coding_effect=annotations[ident]["reference_coding_effect"],
                orfdata_coding_block_overlap=annotations[ident]["orfdata_coding_block_overlap"],
                top_effect_type=annotations[ident]["top_effect_type"],
                top_changed_aa_positions=top.get("changed_residue_count") if top else "",
                top_marked_aa_positions=(int(top["mutation_end_idx"]) - int(top["mutation_start_idx"])) if top else "",
                **{name + "_contexts": context_counts.get(name) for name in
                   ("reference_identical", "amino_acid_divergent", "unresolved", "focal_changed")},
                **{"top_" + key: top.get(key) if top else None for key in
                   ("contributor_fragments", "window_fragments", "single_read_fragments", "paired_template_fragments")},
                top_translation_id=top.get("displayed_translation_id", "") if top else "",
                top_local_protein=top["protein_sequence"] if top else "",
                novel_9mers=score[0] if score else "",
                total_9mers=score[1] if score else "",
                self_fasta_sha256=sources.get("self_fasta_sha256", "")))


def run(args=None, *, prog=None):
    parsed = parser_for_program(parser, prog).parse_args(args)
    if parsed.output.suffix.lower() != ".pdf":
        parser.error("--output must end in .pdf")
    try:
        per_call = parsed.per_call or parsed.output.with_suffix(".tsv")
        evidence_json = parsed.evidence_json or parsed.output.with_suffix(".json")
        outputs = [path.resolve() for path in (parsed.output, per_call, evidence_json)]
        inputs = {Path(path).resolve() for path in
                  [*parsed.vcf, parsed.bam, parsed.self_fasta, parsed.nc_orf_catalog] if path}
        if len(set(outputs)) != len(outputs) or inputs.intersection(outputs):
            raise CommandInputError("Report outputs must be distinct and must not overwrite input files")
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
                variants, alignment, collector, creator, sample_id=parsed.sample_id,
                source=parsed.rna_source or parsed.bam, catalog=catalog)
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
                       dna_source=parsed.dna_source or ", ".join(parsed.vcf),
                       vcf_paths=parsed.vcf, alignment_path=parsed.bam,
                       rna_source=parsed.rna_source or parsed.bam,
                       self_fasta=str(parsed.self_fasta) if parsed.self_fasta else "",
                       self_fasta_sha256=self_fasta_sha256)
        parsed.output.parent.mkdir(parents=True, exist_ok=True)
        per_call.parent.mkdir(parents=True, exist_ok=True)
        write_per_call(per_call, calls, annotations, summary, sources, self_scores)
        evidence_json.parent.mkdir(parents=True, exist_ok=True)
        with evidence_json.open("w") as stream:
            json.dump(dict(schema="isovar.evidence_report.v1", sources=sources, calls=calls,
                           annotations=annotations, proteins=proteins), stream, sort_keys=True)
            stream.write("\n")
        from ..evidence_pdf import write_generic_report
        write_generic_report(parsed.output, calls, annotations, summary, sources,
                             self_scores, catalog_assessed=catalog is not None)
    except (OSError, ValueError, CommandInputError) as error:
        parser_for_program(parser, prog).error(str(error))
    print(parsed.output)
    print(per_call)
    print(evidence_json)
