"""Measure reusable acquisition on the original libraries of a pinned SV audit.

This measures index preparation, validated cache reopening and complete-template
queries. It does not measure assembly or claim that a candidate template supports
an SV. Sources remain separate and temporary indexes are removed after each one.
"""

import argparse
from pathlib import Path
import resource
import shutil
import sys
from tempfile import TemporaryDirectory
from time import perf_counter

from isovar.read_collector import ReadCollector
from isovar.rna_library import prepare_rna_library
from isovar.sv_rna_audit import file_sha256, read_json, write_json


def benchmark(manifest_path: Path, output: Path, event_count: int, source_ids: list[str] | None = None) -> dict:
    """Write independent source metrics for the first declared inventory events."""
    if event_count < 1 or output.exists():
        raise ValueError('Use a positive event count and a new output path')
    manifest = read_json(manifest_path)
    collector = ReadCollector(**manifest.get('read_collection', {}))
    events = manifest['events'][:event_count]
    if source_ids is not None and set(source_ids) - {source['source_id'] for source in manifest['sources']}:
        raise ValueError('Unknown requested source')
    report = dict(status='incomplete', manifest_sha256=file_sha256(manifest_path),
                  selected_source_ids=source_ids,
                  event_ids=[event['event_id'] for event in events], sources=[],
                  scope='Acquisition only; conservative candidates are not SV support')
    output.parent.mkdir(parents=True, exist_ok=True)
    for source in manifest['sources']:
        if source.get('unavailable_reason') or (source_ids is not None and source['source_id'] not in source_ids):
            continue
        queries = [(event['event_id'], view) for event in events
                   if source['source_id'] in event['sources'] and not event.get('unassessed_reason')
                   for view in event['views']]
        references = {view['input']['reference_name'] for _, view in queries}
        if len(references) != 1:
            raise ValueError('Benchmark requires one explicit assembly per source')
        reference = references.pop()
        pin = manifest['files'][source['bam']]
        bam = Path(pin['path'])
        if not bam.is_absolute():
            bam = manifest_path.parent / bam
        with TemporaryDirectory(prefix='rna-index-benchmark-', dir=output.parent) as scratch:
            if shutil.disk_usage(scratch).free < 5 * 1024 ** 3:
                raise ValueError('Benchmark requires a 5 GiB free-space reserve before preparation')
            started = perf_counter()
            prepared = prepare_rna_library(bam, scratch, read_collector=collector, reference_name=reference)
            preparation_seconds = perf_counter() - started
            if prepared.identity['source_sha256'] != pin['sha256']:
                raise ValueError('Original library differs from the pinned manifest')
            print(source['source_id'], 'prepared', preparation_seconds, flush=True)
            started = perf_counter()
            reopened = prepare_rna_library(bam, scratch, read_collector=collector, reference_name=reference)
            cache_reopen_seconds = perf_counter() - started
            rows = []
            with reopened.open() as library:
                for event_id, view in queries:
                    started = perf_counter()
                    templates = records = 0
                    for template in library.templates(view['input']['regions']):
                        templates += 1
                        records += len(template)
                    rows.append(dict(event_id=event_id, view_id=view['view_id'],
                                     seconds=perf_counter() - started,
                                     candidate_templates=templates, distinct_eligible_records=records))
            report['sources'].append(dict(source_id=source['source_id'], identity=prepared.identity,
                statistics=prepared.statistics, preparation_seconds=preparation_seconds,
                validated_cache_reopen_seconds=cache_reopen_seconds,
                index_bytes=sum(path.stat().st_size for path in prepared.directory.iterdir()),
                queries=rows, total_query_seconds=sum(row['seconds'] for row in rows),
                process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                    * (1 if sys.platform == 'darwin' else 1024)))
        write_json(output, report)
        print(source['source_id'], 'queries', len(rows), report['sources'][-1]['total_query_seconds'], flush=True)
    report['status'] = 'complete'
    write_json(output, report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--event-count', type=int, default=12)
    parser.add_argument('--source-id', action='append', help='Limit to this source; repeat for additional libraries')
    options = parser.parse_args()
    benchmark(options.manifest.resolve(), options.output.resolve(), options.event_count, options.source_id)
