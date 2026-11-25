# PYTHON_ARGCOMPLETE_OK
"""Discoverable subcommands, delegating unchanged options to existing handlers."""

import argparse
from copy import copy
from difflib import get_close_matches
from importlib import import_module
import shutil
import sys
import textwrap

import argcomplete

from ..logging import configure_cli_logging, get_logger
from .completion import completing_subcommand, is_completing
from .help import HelpArgumentParser
from .output_args import write_dataframe
from .validation import CommandInputError, check_outputs, input_error

logger = get_logger(__name__)


COMMANDS = {
    "rna-origin": ("isovar_rna_origin", "Compare mitochondrial and nuclear allele-bearing RNA contexts"),
    "rna-features": ("isovar_rna_features", "Count and normalize local RNA support for any antigen category"),
    "variant-rna-paths": ("isovar_variant_rna_paths", "Assemble original ALT/REF RNA paths before coding-model assessment"),
    "sv-source": ("isovar_sv_source", "Validate original SV source definitions and record-specific PCR evidence"),
    "normal-evidence": ("isovar_normal_evidence", "Assess source-pinned matched-normal genotypes and reference confidence"),
    "patient-background-pairs": ("isovar_patient_background_pairs", "Pair mutant RNA proteins with declared patient-background baselines"),
    "rna-proteome": ("isovar_proteome", "Adapt recorded RNA hypotheses to an ORFdata Proteome and FASTA"),
    "cohort-report": ("isovar_cohort", "Reanalyze source-separated RNA libraries and render cohort PDF reports"),
    "calibrate-phasing": ("isovar_calibrate_phasing", "Calibrate phasing error bounds from truth-labeled control BAMs"),
    "viral-rna": ("isovar_viral", "Recover observed viral RNA from combined-reference alignments without a GTF"),
    "isoforms": ("isovar_isoforms", "Reconstruct and quantify CTA or other targeted antigen isoforms"),
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
    "evidence-report": ("isovar_evidence_report",
                        "Build a visual variant/RNA evidence PDF and per-call TSV"),
    "antigen-report": ("isovar_antigen_report",
                       "Combine small-variant, SV, fusion, viral and antigen-isoform evidence"),
    "fusion-rna": ("isovar_fusion", "Validate and translate supplied fusion RNA with junction evidence"),
    "sv-rna": ("isovar_sv_rna", "Reconstruct exploratory RNA paths around one nominated SV"),
    "sv-audit": ("isovar_sv_audit", "Reconstruct and reconcile all nominated SVs in a source-pinned manifest"),
    "allele-interpretations": ("isovar_allele_interpretations",
                               "Compare competing interpretations of one locus with the RNA reads"),
}

# Backward-compatible aliases for commands renamed under the #654 naming
# convention (artifact/domain noun phrases, reserving bare verbs for top-level
# actions like `run`, `plot`, and `calibrate-phasing`).
COMMAND_ALIASES = {
    "fusion": "fusion-rna",
    "viral": "viral-rna",
    "proteome": "rna-proteome",
    "cohort": "cohort-report",
}


# Help sections, in display order. Every command appears in exactly one.
COMMAND_GROUPS = (
    ("Pipeline", ("run",)),
    ("Table exports", (
        "protein-sequences", "translations", "variant-sequences", "variant-reads",
        "allele-reads", "allele-counts", "reference-contexts")),
    ("RNA hypotheses and interpretation", (
        "protein-hypotheses", "variant-rna-paths", "allele-interpretations", "rna-features", "rna-origin")),
    ("Structural variants and fusions", ("sv-source", "sv-rna", "sv-audit", "fusion-rna")),
    ("Viral and isoform evidence", ("viral-rna", "isoforms")),
    ("Proteome and patient background", ("rna-proteome", "patient-background-pairs", "normal-evidence")),
    ("Reports and figures", ("plot", "evidence-report", "antigen-report", "cohort-report")),
    ("Calibration", ("calibrate-phasing",)),
)


def command_listing():
    """Return the command help grouped into sections, wrapped to the terminal."""
    columns = min(shutil.get_terminal_size((80, 24)).columns, 100)
    name_width = max(len(name) for name in COMMANDS) + 4
    lines = []
    for title, names in COMMAND_GROUPS:
        lines.append(title + ":")
        for name in names:
            wrapped = textwrap.wrap(
                COMMANDS[name][1], max(columns - name_width, 20), break_on_hyphens=False)
            lines.append("  %-*s%s" % (name_width - 2, name, wrapped[0]))
            lines.extend(" " * name_width + line for line in wrapped[1:])
        lines.append("")
    lines.append(textwrap.fill(
        "Use 'isovar COMMAND --help' for options. 'isovar --vcf ... --bam ...' runs the pipeline.", columns))
    return "\n".join(lines)


def parser_for_program(parser, prog=None):
    """Give a call its own help/error name without mutating a shared parser."""
    parser = copy(parser)
    if prog is not None:
        parser.prog = prog
    return parser


def load_subcommand_parser(name):
    """Load one subcommand's argument parser without importing other handlers.

    Parameters
    ----------
    name : str
        Subcommand name from :data:`COMMANDS` (or a compatibility alias in
        :data:`COMMAND_ALIASES`).

    Returns
    -------
    argparse.ArgumentParser
        The subcommand's argument parser.
    """
    canonical = COMMAND_ALIASES.get(name, name)
    module_name, _ = COMMANDS[canonical]
    module = import_module("." + module_name, __package__)
    parser = getattr(module, "parser", None)
    if parser is not None:
        return parser_for_program(parser, "isovar " + canonical)
    return module.make_parser("isovar " + canonical)


def build_parser(subcommand_parsers=None):
    """Build the top-level ``isovar`` argument parser without heavy imports.

    Parameters
    ----------
    subcommand_parsers : mapping of str to argparse.ArgumentParser, optional
        Parsers to attach for specific subcommands (used during shell tab
        completion so only the selected subcommand's module is imported).
        When ``None``, a flat ``COMMAND`` positional argument is registered
        and the grouped command listing is shown in the epilog.

    Returns
    -------
    HelpArgumentParser
        Top-level parser for ``isovar``.
    """
    from .. import __version__

    parser = HelpArgumentParser(
        prog="isovar",
        usage="%(prog)s [-h] [--version] COMMAND ...",
        description="Reconstruct mutant protein sequences from RNA evidence.",
        epilog=command_listing(),
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__,
        help="Show the isovar version number and exit",
    )
    if subcommand_parsers is None:
        # The grouped listing in the epilog replaces argparse's flat command list.
        parser.add_argument("command", metavar="COMMAND", help=argparse.SUPPRESS)
        return parser

    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")
    for name, (_, description) in COMMANDS.items():
        if name in subcommand_parsers:
            sub_parser = subcommand_parsers[name]
            has_nested_subparsers = any(
                isinstance(action, argparse._SubParsersAction)
                for action in sub_parser._actions
            )
            if has_nested_subparsers:
                # Directly attach nested-subcommand parsers (e.g. `isovar isoforms`)
                # so nested actions and their flags complete cleanly.
                subparsers._name_parser_map[name] = sub_parser
                subparsers._choices_actions.append(
                    subparsers._ChoicesPseudoAction(name, (), description)
                )
            else:
                subparsers.add_parser(
                    name,
                    help=description,
                    add_help=False,
                    parents=[sub_parser],
                )
        else:
            subparsers.add_parser(name, help=description, add_help=False)
    return parser


def build_completion_parser(environ=None):
    """Build the import-light parser that answers a shell completion request.

    When a known subcommand is already present on ``COMP_LINE``, only that
    subcommand's module is imported and attached so its flags complete.

    Parameters
    ----------
    environ : mapping, optional
        Environment containing argcomplete variables. Defaults to ``os.environ``.

    Returns
    -------
    HelpArgumentParser
        Completion parser ready for :func:`argcomplete.autocomplete`.
    """
    selected = completing_subcommand(COMMANDS, environ)
    if selected is not None:
        return build_parser({selected: load_subcommand_parser(selected)})
    return build_parser({})


def run_dataframe_command(parser, dataframe_from_args, args=None, prog=None):
    """Parse options, build one table from them and write it as CSV or TSV."""
    parser = parser_for_program(parser, prog)
    args = parser.parse_args(sys.argv[1:] if args is None else args)
    configure_cli_logging(args.log_level)
    logger.debug("Options: %s", args)
    try:
        check_outputs(args, "output", allow_stdout=True)
        df = dataframe_from_args(args)
        write_dataframe(df, args)
    except CommandInputError as error:
        input_error(parser, error)
    destination = "stdout" if str(args.output) == "-" else args.output
    logger.info("Wrote %d rows to %s", len(df), destination)


def run(args=None):
    if args is None and is_completing():
        argcomplete.autocomplete(build_completion_parser())

    args = list(sys.argv[1:] if args is None else args)
    if args[:1] == ["help"]:
        # `isovar help` and `isovar help COMMAND`, as in sibling tools.
        args = args[1:2] + ["--help"]
    # Preserve the original option-first interface, including option values
    # which happen to equal a command name. Only the first token dispatches.
    if args and args[0].startswith("-") and args[0] not in {"-h", "--help", "--version"}:
        main = import_module(".isovar_main", __package__)
        option = args[0].split("=", 1)[0]
        known = main.make_parser("isovar")._option_string_actions
        if option.startswith("--") and not any(name.startswith(option) for name in known):
            # A mistyped top-level flag (e.g. --verison) would otherwise report
            # the pipeline's missing --bam.
            close = get_close_matches(option, ["--help", "--version", *sorted(known)], n=1)
            hint = " (did you mean '%s'?)" % close[0] if close else ""
            sys.stderr.write("isovar: error: unrecognized option '%s'%s. Run 'isovar --help' for the "
                             "commands or 'isovar run --help' for the pipeline's options.\n" % (option, hint))
            sys.exit(2)
        return main.run(args, prog="isovar")
    parser = build_parser()
    if not args:
        parser.print_help()
        return
    # Parse only the command; its existing parser owns every remaining option.
    selected = parser.parse_args(args[:1]).command
    if selected in COMMAND_ALIASES:
        canonical = COMMAND_ALIASES[selected]
        sys.stderr.write(
            "isovar: warning: 'isovar %s' is deprecated; use 'isovar %s' instead.\n"
            % (selected, canonical)
        )
        selected = canonical
    if selected not in COMMANDS:
        # Not argparse's `choices`, which would print every command name on one line.
        close = get_close_matches(selected, COMMANDS, n=1)
        hint = " (did you mean '%s'?)" % close[0] if close else ""
        parser.error("argument COMMAND: invalid choice: '%s'%s. Run 'isovar --help' to list "
                     "the commands." % (selected, hint))
    module_name, _ = COMMANDS[selected]
    return import_module("." + module_name, __package__).run(args[1:], prog="isovar " + selected)
