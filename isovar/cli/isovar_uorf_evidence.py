"""Query and export the offline human upstream-ORF reference."""

from orfdata import load_uorf_catalog
from orfdata.cli import query


def run(argv=None, prog=None):
    return query(argv, prog or "isovar-uorf-evidence", catalog_loader=load_uorf_catalog)
