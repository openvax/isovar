"""Render verified shared-library audit snapshots; do not reinterpret RNA."""

import argparse
import csv
import os
from collections import Counter
from pathlib import Path

from isovar.sv_rna_audit import file_sha256, read_json
from .build_unified_sv_manifest import OUTPUT, SV3


def snapshots(directory=OUTPUT / 'results', event_id=None):
    directory = Path(directory)
    report = read_json(directory / 'report.json')
    for row in report['runs']:
        if event_id is not None and row['event_id'] != event_id:
            continue
        path = directory / row['artifact']
        receipt = read_json(str(path) + '.receipt.json')
        if file_sha256(path) != receipt['sha256']:
            raise ValueError('Audit artifact checksum mismatch: ' + row['artifact'])
        snapshot = read_json(path)
        if snapshot['fingerprint'] != receipt['fingerprint']:
            raise ValueError('Audit receipt fingerprint mismatch')
        yield row, snapshot


def render(directory=OUTPUT / 'results', event_id=None, artifact_prefix='results'):
    values = list(snapshots(directory, event_id))
    rows = [row for row, _ in values]
    events = {row['event_id'] for row in rows}
    sources = {value['source_id']: value['source'] for _, value in values}
    with (OUTPUT.parent / 'report_locus_annotations.tsv').open() as stream:
        loci = {r['call_id']: r for r in csv.DictReader(stream, delimiter='\t')
                if r['sample'] == 'COLO829' and r['kind'] == 'sv'}
    lines = ['# ' + ('COLO829 SV3 shared-workflow assessment' if event_id == SV3 else 'COLO829 unified SV RNA audit'), '',
             f'{len(events)} source-defined events; {len(rows)} event/source/orientation runs.', '',
             'Every retained path is assessed by `isovar.assess_sv_rna_paths` and compared by '
             '`isovar.reconcile_sv_rna_hypotheses`. Complete results and original witness provenance '
             'are in the linked compressed snapshots. No report-specific frame inference is performed.', '',
             'Counts below describe retained computational paths, not independent molecules or resolved biological events. '
             'The two SRR8615617 subsets are the same library. Opposite orientations and paired VCF records '
             'are separate hypotheses, not replication. Exact-sequence compatibility does not validate base quality, '
             'establish translation or prove somatic causation.', '',
             'Gene names below overlap nominated DNA breakends under Ensembl 87; they do not assign a coding fusion. '
             'IGV loci are GRCh38 one-based inclusive, one locus per endpoint. '
             'The [source index](../SEQUENCING_SOURCES.md) lists every library and reference.', '',
             '| Source | Library | Culture | Role |', '|---|---|---|---|']
    for identifier, source in sources.items():
        lines.append(f"| {identifier} | {source['library_id']} | {source['culture']} | {source['role']} |")
    states = Counter(row['status'] for row in rows)
    lines += ['', 'Run outcomes: ' + ', '.join(f'{key}: {value}' for key, value in sorted(states.items())) + '.', '',
              'Missing original PacBio tumor/normal evidence remains open in [#574](https://github.com/openvax/isovar/issues/574). '
              'Malformed or sequence-unresolved source definitions remain open in [#575](https://github.com/openvax/isovar/issues/575). '
              'An unassessed or bounded no-candidate result is not evidence of absence.', '',
              '| Event / source / view | Status | Retained paths | Direct event join / compatible join paths | '
              'Event-linked partial / complete protein candidates | Limits |', '|---|---|---:|---:|---:|---|']
    for row, value in values:
        paths = value.get('assessment', {}).get('paths', [])
        direct = sum(any(j['relation'] == 'breakpoint_junction' and j['direct_support'].get('fragments', 0)
                         for j in p['junctions']) for p in paths)
        compatible = sum(any(j['relation'] == 'event_compatible_junction' and j['direct_support'].get('fragments', 0)
                             for j in p['junctions']) for p in paths)
        proteins = Counter(c['completeness'] for p in paths for c in p['proteins'] if c['event_linked'])
        limits = ', '.join(row['limitations']).replace('|', '\\|') or 'none recorded'
        identity = f"{row['event_id']} / {row['source_id']} / {row['view_id']}"
        call_id = row['event_id'].split(':', 1)[1]
        locus = loci[call_id]
        genes = locus['genes'].replace(';', ', ') or 'No overlapping gene in pinned annotation'
        location = locus['igv_loci'].replace(' ', '`<br>`')
        lines.append(f"| [{identity}]({artifact_prefix}/{row['artifact']})<br>**{genes}**<br>IGV: `{location}` | {row['status']} | {len(paths)} | "
                     f"{direct} / {compatible} | {proteins['partial']} / {proteins['complete']} | {limits} |")
    lines += ['', 'Protein candidates retain their annotated frame evidence and completeness. Partial candidates '
              'can assume an upstream annotated frame; exploratory ATGs are reported separately and never become '
              'annotated protein candidates merely because reads tile them. Multiple local hypotheses can coexist.', '',
              'Reproduce with `python -m isovar sv-audit --manifest '
              'reports/hcc1395_colo829/unified_sv_audit/manifest.json --output-dir NEW_DIRECTORY`, then '
              '`python -m reports.hcc1395_colo829.build_unified_sv_report --results NEW_DIRECTORY`.', '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=OUTPUT / 'results')
    parser.add_argument('--output', type=Path, default=OUTPUT / 'REPORT.md')
    parser.add_argument('--event-id')
    args = parser.parse_args()
    args.output.write_text(render(args.results, args.event_id, os.path.relpath(args.results, args.output.parent)))
    print(args.output)


if __name__ == '__main__':
    main()
