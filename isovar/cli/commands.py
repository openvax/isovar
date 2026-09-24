"""Discoverable subcommands, delegating unchanged options to existing handlers."""

import argparse
from copy import copy
from importlib import import_module
import sys

from ..logging import configure_cli_logging, get_logger
from .output_args import write_dataframe
from .validation import CommandInputError, check_output_path

logger = get_logger(__name__)


COMMANDS = {
    "run": ("isovar_main", "Run the complete RNA-to-protein pipeline and result filters"),
    "protein-sequences": ("isovar_protein_sequences", "Export candidate protein sequences"),
    "protein-hypotheses": ("isovar_protein_hypotheses",
                           "Export every protein hypothesis and translation with scoped RNA evidence"),
    "translations": ("isovar_translations", "Export individual RNA translations"),
    "reference-contexts": ("isovar_reference_contexts", "Export reference sequence and reading-frame contexts"),
    "allele-reads": ("isovar_allele_reads", "Export reads overlapping each variant"),
    "allele-counts": ("isovar_allele_counts", "Count reference, alternate and other allele support"),
    "variant-reads": ("isovar_variant_reads", "Export alternate-allele supporting reads"),
    "variant-sequences": ("isovar_variant_sequences", "Export reconstructed cDNA sequences"),
    "plot": ("isovar_plot", "Draw protein, coverage, read-overlap and transcript figures"),
    "fusion": ("isovar_fusion", "Validate and translate supplied fusion RNA with junction evidence"),
    "sv-rna": ("isovar_sv_rna", "Reconstruct exploratory RNA paths around one nominated SV"),
}


def parser_for_program(parser, prog=None):
    """Give a call its own help/error name without mutating a shared parser."""
    parser = copy(parser)
    if prog is not None:
        parser.prog = prog
    return parser


def run_dataframe_command(parser, dataframe_from_args, args=None, prog=None):
    """Parse options, build one table from them and write it as CSV."""
    parser = parser_for_program(parser, prog)
    args = parser.parse_args(sys.argv[1:] if args is None else args)
    configure_cli_logging(args.log_level)
    logger.debug("Options: %s", args)
    try:
        check_output_path(args.output)
        df = dataframe_from_args(args)
        write_dataframe(df, args)
    except CommandInputError as error:
        parser.error(str(error))
    logger.info("Wrote %d rows to %s", len(df), args.output)


def run(args=None):
    from .. import __version__

    args = list(sys.argv[1:] if args is None else args)
    # Preserve the original option-first interface, including option values
    # which happen to equal a command name. Only the first token dispatches.
    if args and args[0].startswith("-") and args[0] not in {"-h", "--help", "--version"}:
        return import_module(".isovar_main", __package__).run(args, prog="isovar")
    parser = argparse.ArgumentParser(
        prog="isovar", description="Reconstruct mutant protein sequences from RNA evidence.",
        epilog=("Use 'isovar COMMAND --help' for options. 'isovar --vcf ... --bam ...' runs the "
                "pipeline, and table/plot commands also install as isovar-COMMAND scripts."))
    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)
    commands = parser.add_subparsers(dest="command", metavar="COMMAND")
    for name, (_, description) in COMMANDS.items():
        commands.add_parser(name, help=description, add_help=False)
    if not args:
        parser.print_help()
        return
    # Parse only the command; its existing parser owns every remaining option.
    selected = parser.parse_args(args[:1]).command
    module_name, _ = COMMANDS[selected]
    return import_module("." + module_name, __package__).run(args[1:], prog="isovar " + selected)
