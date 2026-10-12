"""Compatibility entry point for the human upstream-ORF reference."""

from orfdata.cli import run as _run
from orfdata import load_uorf_catalog


def run(argv=None, prog=None):
    return _run(argv, prog, catalog_loader=load_uorf_catalog)
