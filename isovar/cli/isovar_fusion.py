"""Validate and translate an explicitly supplied fusion transcript/evidence JSON."""
import argparse
import json
from pathlib import Path

from ..default_parameters import FUSION_PEPTIDE_LENGTHS, MIN_FUSION_FRAGMENTS, PLOT_DPI
from ..fusion import fusion_from_dict, reconstruct_fusion
from ..logging import configure_cli_logging


def make_parser(prog="isovar fusion"):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    parser.add_argument("--input", required=True, help="Supplied fusion sequence, mapping, references and RNA evidence JSON")
    parser.add_argument("--output", required=True, help="Evidence-backed fusion result JSON")
    parser.add_argument("--peptide-lengths", type=int, nargs="+", default=FUSION_PEPTIDE_LENGTHS,
                        help="Junction peptide lengths (default: %(default)s)")
    parser.add_argument("--min-fragments", type=int, default=MIN_FUSION_FRAGMENTS,
                        help="Minimum distinct directly junction-spanning fragments (default: %(default)s)")
    parser.add_argument("--plot-dir", help="Also write PNG/SVG/vector PDF panels in a new UTC-stamped directory")
    parser.add_argument("--dpi", type=int, default=PLOT_DPI, help="PNG resolution (default: %(default)s)")
    return parser


def run(args=None, prog=None):
    parser = make_parser(prog or "isovar fusion")
    options = parser.parse_args(args)
    configure_cli_logging()
    try:
        data = json.loads(Path(options.input).expanduser().read_text())
        fusion, references, reads = fusion_from_dict(data)
        result = reconstruct_fusion(fusion, references, reads, peptide_lengths=options.peptide_lengths,
                                    min_fragments=options.min_fragments)
        # Write the result before plotting, so a plotting failure never loses it.
        Path(options.output).expanduser().write_text(json.dumps(result, indent=2) + "\n")
        if options.plot_dir:
            from ..fusion_visualization import save_fusion_figures
            from ..visualization import timestamped_run_directory

            print(save_fusion_figures(result, timestamped_run_directory(options.plot_dir), references,
                                      data.get("reference_names"), dpi=options.dpi))
    except (OSError, ValueError, TypeError, KeyError, ImportError) as error:
        parser.error(str(error))
