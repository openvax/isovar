"""Reconstruct exploratory RNA paths around one nominated SV from an indexed BAM."""
import argparse
import json
from pathlib import Path

from ..default_parameters import (
    FUSION_PEPTIDE_LENGTHS, SV_ASSEMBLE, SV_BREAKPOINT_WINDOW, SV_MAX_PATHS, SV_MAX_QUERIES, SV_MAX_RECORDS,
    SV_MIN_ALTERNATIVE_FRACTION, SV_MIN_ALTERNATIVE_FRAGMENTS, SV_MIN_ANCHOR_BASES, SV_MIN_OVERLAP,
)
from ..sv_rna import reconstruct_sv_rna, sv_rna_input_from_dict
from .rna_args import add_rna_args, alignment_file_from_args, read_collector_from_args


def make_parser(prog="isovar sv-rna"):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    parser.add_argument("--input", required=True,
                        help="Nominated event, search regions, reference models and event provenance JSON")
    parser.add_argument("--output", required=True, help="Exploratory candidate paths JSON")
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
    group.add_argument("--max-records", type=int, default=SV_MAX_RECORDS, help="(default: %(default)s)")
    group.add_argument("--max-queries", type=int, default=SV_MAX_QUERIES, help="(default: %(default)s)")
    group.add_argument("--max-paths", type=int, default=SV_MAX_PATHS, help="(default: %(default)s)")
    group.add_argument("--breakpoint-window", type=int, default=SV_BREAKPOINT_WINDOW,
                       help="Bases searched to either side of each breakpoint (default: %(default)s)")
    group.add_argument("--no-assembly", dest="assemble", action="store_false", default=SV_ASSEMBLE,
                       help="Use only reads which span each seed junction")
    group.add_argument("--peptide-lengths", type=int, nargs="+", default=FUSION_PEPTIDE_LENGTHS,
                       help="Candidate peptide lengths (default: %(default)s)")
    return parser


def run(args=None, prog=None):
    parser = make_parser(prog or "isovar sv-rna")
    options = parser.parse_args(args)
    try:
        inputs = sv_rna_input_from_dict(json.loads(Path(options.input).expanduser().read_text()))
        with alignment_file_from_args(options) as bam:
            result = reconstruct_sv_rna(
                bam, source=options.source or options.bam, read_collector=read_collector_from_args(options),
                min_anchor_bases=options.min_anchor_bases, min_overlap=options.min_overlap,
                min_alternative_fragments=options.min_alternative_fragments,
                min_alternative_fraction=options.min_alternative_fraction, max_records=options.max_records,
                max_queries=options.max_queries, max_paths=options.max_paths, assemble=options.assemble,
                breakpoint_window=options.breakpoint_window,
                peptide_lengths=options.peptide_lengths, **inputs)
        Path(options.output).expanduser().write_text(json.dumps(result, indent=2) + "\n")
    except (OSError, ValueError, TypeError, KeyError) as error:
        parser.error(str(error))
