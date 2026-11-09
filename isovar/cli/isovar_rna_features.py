"""Measure explicit local RNA features from the original alignment records."""

import argparse
import json
from pathlib import Path

from .rna_args import add_rna_path_args, rna_path_context
from .validation import input_error
from .output_args import add_log_level_arg
from ..logging import configure_cli_logging
from ..audit_bundle import publish_json
from ..rna_feature_measurement import measure_rna_path_features
from ..rna_feature_support import RNAFeature, RNASamplingProfile


def make_parser(prog="isovar rna-features"):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    add_rna_path_args(parser)
    parser.add_argument("--genome-fasta", help="Reference FASTA when needed to read CRAM")
    parser.add_argument("--features", required=True,
                        help="JSON with paths (sequence, positions, features); see RNA_FEATURE_SUPPORT.md")
    parser.add_argument("--sample-id", required=True, help="Sample providing the RNA")
    parser.add_argument("--source-id", required=True, help="Stable identity of the original RNA library")
    parser.add_argument("--sampling-profile", help="JSON library geometry and whole-library denominators")
    parser.add_argument("--output", required=True, help="New JSON file with source-separated measurements")
    add_log_level_arg(parser)
    return parser


def run(args=None, prog=None):
    parser = make_parser(prog or "isovar rna-features")
    options = parser.parse_args(args)
    configure_cli_logging(options.log_level)
    try:
        output = Path(options.output)
        if output.exists():
            raise ValueError("Output already exists: %s" % output)
        request = json.loads(Path(options.features).read_text())
        if set(request) != {"paths"} or not request["paths"]:
            raise ValueError("Feature request must contain a nonempty paths list")
        profile = (RNASamplingProfile.from_dict(json.loads(Path(options.sampling_profile).read_text()))
                   if options.sampling_profile else None)
        records = []
        with rna_path_context(options) as (bam, collector, library):
            for path in request["paths"]:
                if set(path) != {"sequence", "positions", "features"}:
                    raise ValueError("Each path requires exactly sequence, positions and features")
                records.extend(measure_rna_path_features(
                    bam, sequence=path["sequence"], positions=path["positions"],
                    features=[RNAFeature(**f) for f in path["features"]],
                    sample_id=options.sample_id, source=options.source_id,
                    profile=profile, read_collector=collector, library_index=library))
        identities = [(r["feature"]["path_id"], r["feature"]["identifier"]) for r in records]
        if len(identities) != len(set(identities)):
            raise ValueError("Duplicate path/feature identity in request")
        publish_json(output, dict(schema="isovar.rna_feature_measurements.v1", measurements=records))
    except (OSError, ValueError, TypeError, KeyError) as error:
        input_error(parser, error)
