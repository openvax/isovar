"""
Export every protein hypothesis for each variant, with the nucleotide
translations behind it and scoped RNA evidence, as JSON plus a TSV with one
row per translation. All proteins are kept unless
--max-protein-sequences-per-variant says otherwise.
"""

import json
import sys

from ..cell_reconstruction import MODES, reconstruct_cell_groups
from ..logging import configure_cli_logging
from ..protein_hypotheses import export_protein_hypotheses, write_protein_hypotheses
from .commands import parser_for_program
from .base_quality_args import add_base_quality_args, base_quality_policy_from_args
from .main_args import make_isovar_arg_parser, run_isovar_from_parsed_args
from .output_args import add_log_level_arg
from .rna_args import alignment_file_from_args, read_collector_from_args
from .validation import CommandInputError, check_output_path

parser = make_isovar_arg_parser(description=__doc__)
add_base_quality_args(parser)
parser.set_defaults(max_protein_sequences_per_variant=0)
add_log_level_arg(parser)
_group = parser.add_argument_group("Protein hypothesis export")
_group.add_argument("--sample-id", required=True, help="Sample the reads came from")
_group.add_argument(
    "--source",
    help="Identity of the read set, which scopes read IDs with --sample-id (default: --bam)")
_group.add_argument(
    "--cell-umi-labels", action="store_true",
    help="Also report cell barcode/UMI labels behind each allele and protein (single-cell data)")
_group.add_argument(
    "--output", default="isovar-protein-hypotheses.json",
    help="JSON output; the TSV is written beside it (default: %(default)s)")

_group.add_argument(
    "--cell-reconstruction", choices=MODES,
    help="Also reconstruct by this pool and attribute original reads to scoped cells (JSON)")
_group.add_argument("--cell-groups", metavar="JSON", help="Scoped cell ID to group-name mapping for group mode")
_group.add_argument("--partial-read-support", action="store_true",
                    help="Also score partial-read support as fractional fragments across compatible hypotheses")
_group.add_argument("--cell-max-edits", "--support-max-edits", type=int,
                    help="Edit tolerance for cell attribution/partial-read scoring (default: 1)")
_group.add_argument("--cell-min-overlap", "--support-min-overlap", type=int,
                    help="Minimum compared bases for cell attribution/partial-read scoring (default: 10)")


def run(args=None, *, prog=None):
    command = parser_for_program(parser, prog)
    args = command.parse_args(sys.argv[1:] if args is None else args)
    configure_cli_logging(args.log_level)
    try:
        check_output_path(args.output)
        quality_policy = base_quality_policy_from_args(args)
        if not args.cell_reconstruction and (args.cell_groups is not None or (
                not args.partial_read_support and any(value is not None for value in
                                                      (args.cell_max_edits, args.cell_min_overlap)))):
            raise ValueError("Cell attribution options require --cell-reconstruction or --partial-read-support")
        max_edits = 1 if args.cell_max_edits is None else args.cell_max_edits
        min_overlap = 10 if args.cell_min_overlap is None else args.cell_min_overlap
        if bool(args.cell_groups) != (args.cell_reconstruction == "group"):
            raise ValueError("--cell-groups is required exactly with --cell-reconstruction group")
        groups = None
        if args.cell_groups:
            with open(args.cell_groups) as handle:
                groups = json.load(handle)
        results = run_isovar_from_parsed_args(args)
        with alignment_file_from_args(args) as alignment_file:
            export = export_protein_hypotheses(
                results, sample_id=args.sample_id, source=args.source or args.bam,
                alignment_header=alignment_file.header.to_dict(),
                cell_umi_alignment_file=alignment_file if args.cell_umi_labels else None,
                read_collector=read_collector_from_args(args), base_quality_policy=quality_policy,
                partial_read_support=args.partial_read_support, support_max_edits=max_edits,
                support_min_overlap=min_overlap)
            if args.cell_reconstruction:
                export["cell_reconstruction"] = reconstruct_cell_groups(
                    results, alignment_file, sample_id=args.sample_id, source=args.source or args.bam,
                    mode=args.cell_reconstruction, cell_groups=groups,
                    max_edits=max_edits, min_overlap=min_overlap,
                    read_collector=read_collector_from_args(args), base_quality_policy=quality_policy,
                    partial_read_support=args.partial_read_support)
        write_protein_hypotheses(export, args.output)
    except (CommandInputError, ValueError, OSError) as error:
        command.error(str(error))
