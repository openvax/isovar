"""Inspect every SV source row before building an RNA audit inventory."""

import argparse
from typing import cast

from ..sv_rna_audit import write_json
from ..sv_sources import SOURCE_PROFILES, SourceProfile, read_sv_source, sv_source_report


def run(args: list[str] | None = None, prog: str | None = None) -> None:
    parser = argparse.ArgumentParser(prog=prog or 'isovar sv-source', description=__doc__)
    parser.add_argument('--vcf', required=True)
    parser.add_argument('--sha256', required=True, help='Checksum of the unmodified source bytes')
    parser.add_argument('--reference-name', required=True)
    parser.add_argument('--profile', choices=SOURCE_PROFILES, default='standard')
    parser.add_argument('--output', required=True)
    options = parser.parse_args(args)
    profile = cast(SourceProfile, options.profile)
    try:
        records = read_sv_source(options.vcf, sha256=options.sha256,
                                reference_name=options.reference_name, profile=profile)
        write_json(options.output, sv_source_report(records, sha256=options.sha256,
                                                    reference_name=options.reference_name, profile=profile))
    except (OSError, ValueError) as error:
        parser.error(str(error))
