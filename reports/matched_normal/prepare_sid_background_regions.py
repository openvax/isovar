"""Pin the exon/ncORF territory used by Sid's previous supplied background.

The library separately adds nomination flanks. This file captures distant
transcript regions and complete REF spans from previously supplied alleles, so
expanded normal assessment does not silently shrink the earlier scope.
"""

import argparse
from pathlib import Path
import sqlite3

import orfdata
from pyensembl import EnsemblRelease
from varcode import normalize_contig_name

from isovar.genomic_intervals import GenomicIntervalIndex
from isovar.sv_rna_audit import file_sha256, read_json, write_json


def catalog_regions(catalog_path: Path) -> list[tuple[str, int, int]]:
    """Convert ORFdata's 1-based inclusive blocks, including one-base blocks."""
    return [(record.contig, start - 1, end) for record in orfdata.load_catalog(catalog_path).records
            for start, end in record.coding_blocks]


def prepare_regions(catalog_path: Path, previous_background: Path, output: Path) -> dict:
    """Write merged half-open regions and a receipt pinning every input's bytes."""
    receipt_path = output.with_name(output.name + '.receipt.json')
    if output.exists() or receipt_path.exists():
        raise ValueError('Background scope already exists; choose a new output path')
    regions, sources = [], []
    for release in (87, 115):
        database = Path(EnsemblRelease(release).db.local_db_path)
        if not database.is_file():
            raise ValueError('Install the Ensembl annotation before preparing scope: ' + str(release))
        with sqlite3.connect(database.as_uri() + '?mode=ro', uri=True) as connection:
            exons = [(str(contig), start - 1, end) for contig, start, end in
                     connection.execute('SELECT seqname,start,end FROM exon')]
        regions.extend(exons)
        sources.append(dict(kind='ensembl_exons', release=release, path=str(database),
                            sha256=file_sha256(database), intervals=len(exons)))
    blocks = catalog_regions(catalog_path)
    regions.extend(blocks)
    sources.append(dict(kind='orfdata_coding_blocks', path=str(catalog_path), sha256=file_sha256(catalog_path), intervals=len(blocks)))
    previous = [(row['contig'], row['start'] - 1, row['start'] - 1 + max(1, len(row['ref'])))
                for row in read_json(previous_background)]
    regions.extend(previous)
    sources.append(dict(kind='previous_background_ref_spans', path=str(previous_background),
                        sha256=file_sha256(previous_background), intervals=len(previous)))
    names = {contig: normalize_contig_name(contig, reference_name='GRCh38') for contig in {r[0] for r in regions}}
    merged = GenomicIntervalIndex((names[contig], start, end)
                                   for contig, start, end in regions).intervals()
    output.parent.mkdir(parents=True, exist_ok=True)
    write_json(output, merged)
    receipt = dict(schema='isovar.background_region_preparation.v1', reference_name='GRCh38', sources=sources,
        output=dict(path=str(output), sha256=file_sha256(output), intervals=len(merged),
                    bases=sum(end - start for _, start, end in merged)),
        coordinate_system='zero-based half-open',
        interpretation='Requested background territory only; no genotype or reference confidence is inferred from annotation.')
    write_json(receipt_path, receipt)
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ncorf-catalog', type=Path, required=True)
    parser.add_argument('--previous-background', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    options = parser.parse_args()
    print(prepare_regions(options.ncorf_catalog.resolve(), options.previous_background.resolve(), options.output.resolve())['output'])
