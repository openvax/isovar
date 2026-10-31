"""Discoverable subcommands, delegating unchanged options to existing handlers."""

import argparse
from copy import copy
from difflib import get_close_matches
from importlib import import_module
import shutil
import sys
import textwrap

from ..logging import configure_cli_logging, get_logger
from .output_args import write_dataframe
from .validation import CommandInputError, check_output_path

logger = get_logger(__name__)


COMMANDS = {
    "variant-rna-paths": ("isovar_variant_rna_paths", "Assemble original ALT/REF RNA paths before coding-model assessment"),
    "sv-source": ("isovar_sv_source", "Validate original SV source definitions and record-specific PCR evidence"),
    "normal-evidence": ("isovar_normal_evidence", "Assess source-pinned matched-normal genotypes and reference confidence"),
    "patient-background-pairs": ("isovar_patient_background_pairs", "Pair mutant RNA proteins with declared patient-background baselines"),
    "proteome": ("isovar_proteome", "Adapt recorded RNA hypotheses to an ORFdata Proteome and FASTA"),
    "cohort": ("isovar_cohort", "Reanalyze source-separated RNA libraries and render cohort PDF reports"),
    "calibrate-phasing": ("isovar_calibrate_phasing", "Calibrate phasing error bounds from truth-labeled control BAMs"),
    "viral": ("isovar_viral", "Recover observed viral RNA from combined-reference alignments without a GTF"),
    "isoforms": ("isovar_isoforms", "Reconstruct and quantify CTA or other targeted antigen isoforms"),
    "uorf-evidence": ("isovar_uorf_evidence", "Query and export the offline human uORF evidence reference"),
    "audit-bundle": ("isovar_audit_bundle", "Create, verify, archive and replay portable offline audit inputs"),
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
    "fusion": ("isovar_fusion", "Validate and translate supplied fusion RNA with junction evidence"),
    "sv-rna": ("isovar_sv_rna", "Reconstruct exploratory RNA paths around one nominated SV"),
    "sv-audit": ("isovar_sv_audit", "Reconstruct and reconcile all nominated SVs in a source-pinned manifest"),
    "allele-interpretations": ("isovar_allele_interpretations",
                               "Compare competing interpretations of one locus with the RNA reads"),
}


# Help sections, in display order. Every command appears in exactly one.
COMMAND_GROUPS = (
    ("Pipeline", ("run",)),
    ("Table exports", (
        "protein-sequences", "translations", "variant-sequences", "variant-reads",
        "allele-reads", "allele-counts", "reference-contexts")),
    ("RNA hypotheses and interpretation", (
        "protein-hypotheses", "variant-rna-paths", "allele-interpretations")),
    ("Structural variants and fusions", ("sv-source", "sv-rna", "sv-audit", "fusion")),
    ("Viral, isoform and uORF evidence", ("viral", "isoforms", "uorf-evidence")),
    ("Proteome and patient background", ("proteome", "patient-background-pairs", "normal-evidence")),
    ("Reports and figures", ("plot", "evidence-report", "antigen-report", "cohort")),
    ("Audit and calibration", ("audit-bundle", "calibrate-phasing")),
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
        "Use 'isovar COMMAND --help' for options. 'isovar --vcf ... --bam ...' runs the "
        "pipeline, and table/plot commands also install as isovar-COMMAND scripts.", columns))
    return "\n".join(lines)


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
        prog="isovar", usage="%(prog)s [-h] [--version] COMMAND ...",
        description="Reconstruct mutant protein sequences from RNA evidence.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=command_listing())
    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)
    # The grouped listing in the epilog replaces argparse's flat command list.
    parser.add_argument("command", metavar="COMMAND", help=argparse.SUPPRESS)
    if not args:
        parser.print_help()
        return
    # Parse only the command; its existing parser owns every remaining option.
    selected = parser.parse_args(args[:1]).command
    if selected not in COMMANDS:
        # Not argparse's `choices`, which would print all 27 names on one line.
        close = get_close_matches(selected, COMMANDS, n=1)
        hint = " (did you mean '%s'?)" % close[0] if close else ""
        parser.error("argument COMMAND: invalid choice: '%s'%s. Run 'isovar --help' to list "
                     "the commands." % (selected, hint))
    module_name, _ = COMMANDS[selected]
    return import_module("." + module_name, __package__).run(args[1:], prog="isovar " + selected)
