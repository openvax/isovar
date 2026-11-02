"""Inspect every SV source row before building an RNA audit inventory."""

import argparse
from typing import cast

from .output_args import add_log_level_arg
from ..logging import configure_cli_logging, get_logger
from ..sv_rna_audit import write_json
from ..sv_sources import SOURCE_PROFILES, SourceProfile, read_sv_source, sv_source_report
from .validation import check_outputs, input_error

logger = get_logger(__name__)


def run(args: list[str] | None = None, prog: str | None = None) -> None:
    parser = argparse.ArgumentParser(prog=prog or 'isovar sv-source', description=__doc__)
    parser.add_argument('--vcf', required=True, help='Original SV caller VCF (plain or gzip), unmodified')
    parser.add_argument('--sha256', required=True, help='Checksum of the unmodified source bytes')
    parser.add_argument('--reference-name', required=True,
                        help='Assembly the VCF coordinates refer to, such as GRCh38')
    parser.add_argument('--profile', choices=SOURCE_PROFILES, default='standard',
                        help='Interpretation profile; named profiles require their pinned original '
                             '(default: %(default)s)')
    parser.add_argument('--output', required=True,
                        help='Per-row source report JSON, gzip-compressed if it ends in .gz')
    add_log_level_arg(parser)
    options = parser.parse_args(args)
    configure_cli_logging(options.log_level)
    profile = cast(SourceProfile, options.profile)
    try:
        check_outputs(options, 'output')
        records = read_sv_source(options.vcf, sha256=options.sha256,
                                reference_name=options.reference_name, profile=profile)
        write_json(options.output, sv_source_report(records, sha256=options.sha256,
                                                    reference_name=options.reference_name, profile=profile))
        logger.info("Wrote %d source rows to %s", len(records), options.output)
    except (OSError, ValueError) as error:
        input_error(parser, error)
