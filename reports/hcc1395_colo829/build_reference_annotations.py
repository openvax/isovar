"""Annotate every HCC1395 and COLO829 atlas call against pinned reference ORFs.

Run with an environment containing Varcode, PyEnsembl, and orfdata 0.34.1.
The bundled orfdata catalog checksum is validated before annotation.
"""

import csv
import gzip
import json
import argparse
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path

import orfdata
from pyensembl import cached_release
from varcode import Variant


ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / 'reports/hcc1395_colo829'
OUTPUT = REPORT / 'reference_coding_orfdata_annotations.tsv.gz'
CATALOG = Path(orfdata.__file__).parent / 'data/orf-evidence/catalog.json.gz'
CATALOG_SHA256 = '06bc0bf3f81c25a7356f71a38581ff5704c05dcb75752454bc1ca75075ab8c09'
FIELDS = (
    'sample', 'dna_source', 'class', 'call_id', 'chrom', 'pos', 'ref', 'alt',
    'dna_filter', 'genomic_transcript_overlap', 'reference_coding_effect',
    'top_effect_type', 'top_effect_description', 'orfdata_coding_block_overlap',
    'orfdata_orf_count', 'orfdata_orf_ids', 'orfdata_biotypes',
    'orfdata_overlap_rule', 'annotation_error',
)


def read_calls():
    for filename in ('hcc1395_direct_rna_atlas_calls.tsv.gz',
                     'colo829_atlas_calls.tsv.gz'):
        with gzip.open(REPORT / filename, 'rt', newline='') as source:
            yield from csv.DictReader(source, delimiter='\t')


def affected_interval(call):
    """Return a zero-based half-open reference interval for an event.

    For insertions, query the VCF anchor base as a candidate overlap. That
    alone does not establish an insertion within a coding ORF.
    """
    pos = int(call['pos'])
    ref, alt = call['ref'], call['alt']
    if len(ref) > len(alt):
        return pos, pos + len(ref) - 1, 'deleted_reference_bases'
    if len(alt) > len(ref):
        return pos - 1, pos, 'insertion_anchor_candidate'
    return pos - 1, pos - 1 + len(ref), 'changed_reference_bases'


def catalog_index(records):
    """Index coding blocks with prefix maxima for reference-span queries."""
    blocks = defaultdict(list)
    for record in records:
        for first, last in record.coding_blocks:
            blocks[record.contig.removeprefix('chr')].append((first - 1, last, record.orf_id))
    result = {}
    for contig, intervals in blocks.items():
        intervals.sort()
        maximum = 0
        maxima = []
        for _, end, _ in intervals:
            maximum = max(maximum, end)
            maxima.append(maximum)
        result[contig] = ([x[0] for x in intervals], maxima, intervals)
    return result


def overlapping_orfs(index, contig, begin, end):
    """Return ORF IDs whose coding blocks intersect a half-open interval."""
    if contig not in index:
        return []
    starts, maxima, intervals = index[contig]
    i = bisect_left(starts, end) - 1
    ids = set()
    while i >= 0 and maxima[i] > begin:
        _, block_end, orf_id = intervals[i]
        if block_end > begin:
            ids.add(orf_id)
        i -= 1
    return sorted(ids)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    if orfdata.__version__ != '0.34.1':
        raise RuntimeError('This audit requires orfdata 0.34.1')
    catalog = orfdata.load_catalog(path=CATALOG, expected_sha256=CATALOG_SHA256)
    index = catalog_index(catalog.records)
    by_id = {r.orf_id: r for r in catalog.records}
    genome = cached_release(87)
    count = errors = 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.output, 'wt', newline='') as target:
        writer = csv.DictWriter(target, fieldnames=FIELDS, delimiter='\t',
                                lineterminator='\n')
        writer.writeheader()
        for call in read_calls():
            contig = call['chrom'].removeprefix('chr')
            begin, end, rule = affected_interval(call)
            ids = overlapping_orfs(index, contig, begin, end)
            result = {
                'sample': call['sample'], 'dna_source': call['dna_source'],
                'class': call['class'], 'call_id': call['call_id'],
                'chrom': call['chrom'], 'pos': call['pos'], 'ref': call['ref'],
                'alt': call['alt'], 'dna_filter': call['dna_filter'],
                'genomic_transcript_overlap': bool(call['genomic_overlapping_transcript_ids']),
                'reference_coding_effect': False,
                'top_effect_type': '', 'top_effect_description': '',
                'orfdata_coding_block_overlap': bool(ids),
                'orfdata_orf_count': len(ids), 'orfdata_orf_ids': ';'.join(ids),
                'orfdata_biotypes': ';'.join(sorted({by_id[i].biotype for i in ids})),
                'orfdata_overlap_rule': rule, 'annotation_error': '',
            }
            # An indel can affect an exon after its VCF anchor, even if the
            # anchor itself lies outside a transcript interval.
            if result['genomic_transcript_overlap'] or call['class'] == 'indel':
                try:
                    effects = Variant(contig, int(call['pos']), call['ref'],
                                      call['alt'], genome).effects()
                    top = effects.top_priority_effect()
                    result['reference_coding_effect'] = any(
                        effect.modifies_coding_sequence for effect in effects)
                    result['top_effect_type'] = type(top).__name__
                    result['top_effect_description'] = top.short_description
                except Exception as exc:
                    result['annotation_error'] = f'{type(exc).__name__}: {exc}'
                    errors += 1
            writer.writerow(result)
            count += 1
            if count % 5000 == 0:
                print('annotated', count, 'errors', errors, flush=True)
    metadata = {
        'records': count,
        'errors': errors,
        'orfdata_version': orfdata.__version__,
        'orfdata_catalog_records': len(catalog.records),
        'orfdata_catalog_version': catalog.metadata['dataset_version'],
        'orfdata_catalog_sha256': CATALOG_SHA256,
        'genome_annotation': 'Ensembl 87',
        'reference_build': 'GRCh38',
    }
    args.output.with_suffix('.metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(metadata, flush=True)
    if errors:
        raise RuntimeError(f'{errors} annotation errors')


if __name__ == '__main__':
    main()
