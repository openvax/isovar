"""Validate and translate an explicitly supplied fusion transcript/evidence JSON."""
import argparse

from ..default_parameters import FUSION_PEPTIDE_LENGTHS, MIN_FUSION_FRAGMENTS, PLOT_DPI
from ..fusion import fusion_from_dict, reconstruct_fusion
from ..json_files import read_json, write_json
from ..logging import configure_cli_logging, get_logger
from .output_args import add_log_level_arg
from .validation import CommandInputError, check_outputs, dpi, input_error

logger = get_logger(__name__)


def make_parser(prog="isovar fusion"):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    parser.add_argument("--input", required=True,
                        help="Supplied fusion sequence, mapping, references and RNA evidence JSON (or .json.gz)")
    parser.add_argument("--output", required=True,
                        help="Evidence-backed fusion result JSON, gzip-compressed if it ends in .gz")
    parser.add_argument("--peptide-lengths", type=int, nargs="+", default=FUSION_PEPTIDE_LENGTHS,
                        help="Junction peptide lengths (default: %(default)s)")
    parser.add_argument("--min-fragments", type=int, default=MIN_FUSION_FRAGMENTS,
                        help="Minimum distinct directly junction-spanning fragments (default: %(default)s)")
    parser.add_argument("--plot-dir", help="Also write PNG/SVG/vector PDF panels in a new UTC-stamped directory")
    parser.add_argument("--dpi", type=dpi, default=PLOT_DPI, help="PNG resolution (default: %(default)s)")
    add_log_level_arg(parser)
    return parser


def run(args=None, prog=None):
    parser = make_parser(prog or "isovar fusion")
    options = parser.parse_args(args)
    configure_cli_logging(options.log_level)
    try:
        check_outputs(options, "output", directories=("plot_dir",))
        data = read_json(options.input)
        fusion, references, reads = fusion_from_dict(data)
        result = reconstruct_fusion(fusion, references, reads, peptide_lengths=options.peptide_lengths,
                                    min_fragments=options.min_fragments)
        # Write the result before plotting, so a plotting failure never loses it.
        write_json(options.output, result)
        logger.info("Wrote %s", options.output)
        if options.plot_dir:
            from ..fusion_visualization import save_fusion_figures
            from ..visualization import timestamped_run_directory

            print(save_fusion_figures(result, timestamped_run_directory(options.plot_dir), references,
                                      data.get("reference_names"), dpi=options.dpi))
    except (CommandInputError, OSError, ValueError, ImportError) as error:
        input_error(parser, error)
