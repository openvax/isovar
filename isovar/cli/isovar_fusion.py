"""Validate and translate an explicitly supplied fusion transcript/evidence JSON."""
import argparse
import json
from pathlib import Path

from ..default_parameters import FUSION_PEPTIDE_LENGTHS, MIN_FUSION_FRAGMENTS
from ..fusion import fusion_from_dict, reconstruct_fusion


def make_parser(prog="isovar fusion"):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    parser.add_argument("--input", required=True, help="Supplied fusion sequence, mapping, references and RNA evidence JSON")
    parser.add_argument("--output", required=True, help="Evidence-backed fusion result JSON")
    parser.add_argument("--peptide-lengths", type=int, nargs="+", default=FUSION_PEPTIDE_LENGTHS,
                        help="Junction peptide lengths (default: %(default)s)")
    parser.add_argument("--min-fragments", type=int, default=MIN_FUSION_FRAGMENTS,
                        help="Minimum distinct directly junction-spanning fragments (default: %(default)s)")
    return parser


def run(args=None, prog=None):
    parser = make_parser(prog or "isovar fusion")
    options = parser.parse_args(args)
    try:
        data = json.loads(Path(options.input).expanduser().read_text())
        fusion, references, reads = fusion_from_dict(data)
        result = reconstruct_fusion(fusion, references, reads, peptide_lengths=options.peptide_lengths,
                                    min_fragments=options.min_fragments)
        Path(options.output).expanduser().write_text(json.dumps(result, indent=2) + "\n")
    except (OSError, ValueError, TypeError, KeyError) as error:
        parser.error(str(error))
