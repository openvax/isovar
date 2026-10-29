"""Reproduce source-file comparisons and exact reference checks, not RNA claims."""

import argparse
from collections import Counter
from pathlib import Path
from itertools import zip_longest

import pysam
from varcode import ContigResolver

from isovar.sv_rna_audit import file_sha256, write_json
from isovar.sv_sources import SOURCE_PROFILES, read_sv_source


def source_rows(path: Path) -> dict[str, list[str]]:
    rows = [line.split('\t') for line in path.read_text().splitlines() if line and not line.startswith('#')]
    result = {row[2]: row for row in rows}
    if len(result) != len(rows):
        raise ValueError('This source comparison requires unique source record IDs')
    return result


def compare(hcc: Path, author_hcc: Path, colo: Path, minda_colo: Path, reference: Path) -> dict:
    """Compare exact rows, corresponding breakend alleles and local REF matches."""
    hcc_records = read_sv_source(hcc, sha256=SOURCE_PROFILES['hcc1395_svscope_2026'],
                                 reference_name='GRCh38', profile='hcc1395_svscope_2026')
    colo_records = read_sv_source(colo, sha256=SOURCE_PROFILES['colo829_paulin_paz_v2'],
                                  reference_name='GRCh38', profile='colo829_paulin_paz_v2')
    current, previous = source_rows(hcc), source_rows(author_hcc)
    shared = current.keys() & previous.keys()
    changed_fields = Counter(tuple(index for index, (left, right) in enumerate(zip_longest(current[key], previous[key]))
                                   if left != right) for key in shared)
    reference_checks = []
    with pysam.FastaFile(str(reference)) as fasta:
        resolver = ContigResolver(fasta.references, reference_name='GRCh38')
        for record in hcc_records:
            fields = current[record.record_id]
            if set(fields[3] + fields[4]) - set('ACGTN'):
                continue
            resolution = resolver.resolve(fields[0])
            if resolution.resolved_name is None:
                raise ValueError(resolution.diagnostic)
            position, allele = int(fields[1]) - 1, fields[3]
            start = max(0, position - 500)
            end = min(fasta.get_reference_length(resolution.resolved_name), position + len(allele) + 500)
            sequence = fasta.fetch(resolution.resolved_name, start, end).upper()
            offsets, cursor = [], 0
            while (found := sequence.find(allele, cursor)) >= 0:
                offsets.append(start + found - position)
                cursor = found + 1
            reference_checks.append(dict(record_id=record.record_id, contig=resolution.resolved_name,
                vcf_position_1based=position + 1, reference_allele_length=len(allele),
                exact_reference_match_offsets_from_vcf_position=offsets))
    amended = list(source_rows(minda_colo).values())
    correspondences = []
    for record in colo_records:
        if not record.transformations:
            continue
        original = record.original_row.split('\t')
        matches = [row[2] for row in amended
                   if row[:2] == original[:2] and row[3:5] == [original[3], record.interpreted_alt]]
        correspondences.append(dict(original_record_id=record.record_id,
                                     amended_coordinate_and_allele_matches=matches))
    inputs = dict(hcc=hcc, author_hcc=author_hcc, colo=colo, minda_colo=minda_colo, reference=reference)
    return dict(schema='isovar.sv_source_comparison.v1',
        driver_sha256=file_sha256(__file__), inputs={key: dict(filename=path.name, sha256=file_sha256(path))
                                                   for key, path in inputs.items()},
        hcc_row_comparison=dict(current_rows=len(current), previous_rows=len(previous),
            shared_rows=len(shared), added_ids=sorted(current.keys() - previous.keys()),
            removed_ids=sorted(previous.keys() - current.keys()),
            changed_column_indexes=[dict(indexes=list(key), count=count) for key, count in sorted(changed_fields.items())]),
        hcc_literal_reference_checks=reference_checks, reference_search_flank_bases=500,
        colo_amended_correspondences=correspondences,
        interpretation='Exact source/reference comparison only; no coordinate repair or RNA evidence inferred')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('hcc', 'author-hcc', 'colo', 'minda-colo', 'reference', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    write_json(args.output, compare(args.hcc, args.author_hcc, args.colo, args.minda_colo, args.reference))
