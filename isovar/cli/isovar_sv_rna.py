"""Reconstruct exploratory RNA paths around one nominated SV from an indexed BAM."""
import argparse
import json
from pathlib import Path

from ..default_parameters import (
    FUSION_PEPTIDE_LENGTHS, SV_ASSEMBLE, SV_BREAKPOINT_WINDOW, SV_MAX_BREAKPOINT_SHIFT,
    SV_MAX_EXTENSION_SEGMENTS, SV_MAX_PATHS, SV_ANNOTATED_JUNCTION_TOLERANCE, SV_MAX_QUERIES, SV_MAX_RECORDS,
    SV_MIN_ALTERNATIVE_FRACTION, SV_MIN_ALTERNATIVE_FRAGMENTS, SV_MIN_LOCAL_VARIANT_FRACTION, SV_MIN_ANCHOR_BASES, SV_MIN_OVERLAP,
    SV_MIN_ORF_AMINO_ACIDS, SV_MAX_ORF_CANDIDATES,
)
from ..sv_rna import reconstruct_sv_rna, sv_rna_input_from_dict
from ..logging import configure_cli_logging
from ..sv_rna_orf_export import export_sv_rna_orfs, sv_rna_orf_output_paths, write_sv_rna_orfs
from .rna_args import add_rna_args, alignment_file_from_args, read_collector_from_args


def make_parser(prog="isovar sv-rna"):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    parser.add_argument("--input", required=True,
                        help="Nominated event, search regions, reference models and event provenance JSON")
    parser.add_argument("--output", required=True, help="Exploratory candidate paths JSON")
    parser.add_argument("--orf-output-prefix", help="Also export exploratory ORFs as JSON, TSV and protein/DNA FASTA")
    parser.add_argument("--predictions", type=Path,
                        help="Optional source-identified protein predictions JSON for sequence comparison")
    parser.add_argument("--source", help="Alignment source identity recorded in the output (default: --bam)")
    add_rna_args(parser)
    group = parser.add_argument_group("SV RNA reconstruction")
    group.add_argument("--min-anchor-bases", type=int, default=SV_MIN_ANCHOR_BASES,
                       help="Exact collinear CDS match needed to transfer a frame (default: %(default)s)")
    group.add_argument("--min-overlap", type=int, default=SV_MIN_OVERLAP,
                       help="Exact coordinate-anchored overlap needed to extend a path (default: %(default)s)")
    group.add_argument("--min-alternative-fragments", type=int, default=SV_MIN_ALTERNATIVE_FRAGMENTS,
                       help="Fragments for an unattributed junction to seed a path, or for a branch to "
                            "prune a weaker alternative (default: %(default)s)")
    group.add_argument("--min-alternative-fraction", type=float, default=SV_MIN_ALTERNATIVE_FRACTION,
                       help="Prune alternatives below this fraction of the best-supported one "
                            "(default: %(default)s)")
    group.add_argument("--min-local-variant-fraction", type=float, default=SV_MIN_LOCAL_VARIANT_FRACTION,
                       help="Base/small-indel alternatives also need this fraction of the best (default: %(default)s)")
    group.add_argument("--max-records", type=int, default=SV_MAX_RECORDS,
                       help="Alignment records examined; reaching the limit is reported (default: %(default)s)")
    group.add_argument("--max-queries", type=int, default=SV_MAX_QUERIES,
                       help="Region/mate/supplementary lookups; reaching the limit is reported (default: %(default)s)")
    group.add_argument("--max-paths", type=int, default=SV_MAX_PATHS,
                       help="Candidate RNA paths kept; reaching the limit is reported (default: %(default)s)")
    group.add_argument("--max-extension-segments", type=int, default=SV_MAX_EXTENSION_SEGMENTS,
                       help="Segments built to extend one path end, furthest-reaching first (default: %(default)s)")
    group.add_argument("--breakpoint-window", type=int, default=SV_BREAKPOINT_WINDOW,
                       help="Bases searched to either side of each breakpoint (default: %(default)s)")
    group.add_argument("--annotated-junction-tolerance", type=int, default=SV_ANNOTATED_JUNCTION_TOLERANCE,
                       help="Treat unannotated joins this close to annotated ones as wobble (default: %(default)s)")
    group.add_argument("--max-breakpoint-shift", type=int, default=SV_MAX_BREAKPOINT_SHIFT,
                       help="Junction bases an aligner may place past a breakpoint (default: %(default)s)")
    group.add_argument("--no-assembly", dest="assemble", action="store_false", default=SV_ASSEMBLE,
                       help="Use only reads which span each seed junction")
    group.add_argument("--peptide-lengths", type=int, nargs="+", default=FUSION_PEPTIDE_LENGTHS,
                       help="Candidate peptide lengths (default: %(default)s)")
    group.add_argument("--min-orf-amino-acids", type=int, default=SV_MIN_ORF_AMINO_ACIDS,
                       help="Minimum separate exploratory ATG ORF length (default: %(default)s)")
    group.add_argument("--max-orf-candidates", type=int, default=SV_MAX_ORF_CANDIDATES,
                       help="Per-path exploratory ORF cap, with explicit truncation (default: %(default)s)")
    return parser


def run(args=None, prog=None):
    parser = make_parser(prog or "isovar sv-rna")
    options = parser.parse_args(args)
    configure_cli_logging()
    try:
        if options.orf_output_prefix and Path(options.output).expanduser().resolve() in {
                path.resolve() for path in sv_rna_orf_output_paths(options.orf_output_prefix).values()}:
            raise ValueError("ORF output prefix collides with --output")
        inputs = sv_rna_input_from_dict(json.loads(Path(options.input).expanduser().read_text()))
        predictions = json.loads(options.predictions.expanduser().read_text()) if options.predictions else None
        with alignment_file_from_args(options) as bam:
            result = reconstruct_sv_rna(
                bam, source=options.source or options.bam, read_collector=read_collector_from_args(options),
                min_anchor_bases=options.min_anchor_bases, min_overlap=options.min_overlap,
                min_alternative_fragments=options.min_alternative_fragments,
                min_alternative_fraction=options.min_alternative_fraction,
                min_local_variant_fraction=options.min_local_variant_fraction, max_records=options.max_records,
                max_queries=options.max_queries, max_paths=options.max_paths, assemble=options.assemble,
                breakpoint_window=options.breakpoint_window, max_breakpoint_shift=options.max_breakpoint_shift,
                max_extension_segments=options.max_extension_segments,
                annotated_junction_tolerance=options.annotated_junction_tolerance,
                peptide_lengths=options.peptide_lengths, min_orf_amino_acids=options.min_orf_amino_acids,
                max_orf_candidates=options.max_orf_candidates, **inputs)
        if options.predictions:
            from ..sv_rna_comparison import compare_sv_rna_predictions
            result["prediction_comparison"] = compare_sv_rna_predictions(
                result, predictions)
        Path(options.output).expanduser().write_text(json.dumps(result, indent=2) + "\n")
        if options.orf_output_prefix:
            write_sv_rna_orfs(export_sv_rna_orfs(result), options.orf_output_prefix)
    except (OSError, ValueError, TypeError, KeyError) as error:
        parser.error(str(error))
