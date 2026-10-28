"""Acquire pinned regional normal VCF derivatives without retaining full downloads.

Run as a module from the repository root. The reviewed inventory selects complete
public source IDs, exact VCF columns and policies. It does not infer specimen
identity from filenames. Original HTTP bytes are hashed while HTSlib parses the
stream; all records overlapping the declared regions are retained, including
no-calls and rejected genotypes. Scientific interpretation belongs to the shared
``isovar normal-evidence`` command.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from hashlib import sha256
import os
from pathlib import Path
from threading import Thread
from urllib.request import Request, urlopen

import pysam
from varcode import ContigResolver, Variant

from isovar.genomic_intervals import GenomicIntervalIndex
from isovar.normal_sources import normal_query_interval
from isovar.sv_rna_audit import file_sha256, read_json, write_json


@contextmanager
def _hashed_vcf_stream(source):
    """Feed one conditional HTTP response to HTSlib and hash its original bytes."""
    provenance = source['provenance']
    request = Request(provenance['url'], headers={'If-Match': provenance['http_identity']['ETag']})
    receipt = dict(bytes=0, sha256=None)
    failures = []
    with urlopen(request, timeout=120) as response:
        if response.headers.get('ETag') != provenance['http_identity']['ETag']:
            raise ValueError('Source HTTP identity changed: ' + source['source_id'])
        read_fd, write_fd = os.pipe()

        def feed():
            digest = sha256()
            try:
                with os.fdopen(write_fd, 'wb') as destination:
                    while block := response.read(1024 * 1024):
                        digest.update(block)
                        receipt['bytes'] += len(block)
                        destination.write(block)
                receipt['sha256'] = digest.hexdigest()
            except Exception as error:
                failures.append(error)

        thread = Thread(target=feed, daemon=True)
        thread.start()
        try:
            with pysam.VariantFile(read_fd) as reader:
                yield reader, receipt
        finally:
            os.close(read_fd)
            thread.join()
        if failures:
            raise failures[0]
        if receipt['bytes'] != provenance['source_size']:
            raise ValueError('Source size changed or download incomplete: ' + source['source_id'])


def acquire(source, variants, flank, output):
    """Select record intervals with one resolver and one merged interval index."""
    destination = output / (source['source_id'] + '.vcf.gz')
    receipt_path = output / (source['source_id'] + '.acquisition.json.gz')
    if receipt_path.exists():
        result = read_json(receipt_path)
        if result['selection'] != source or result['flank'] != flank or result['nomination_sha256'] != variants['sha256']:
            raise ValueError('Acquisition scope changed; use a new directory')
        for name in ('path', 'index_path'):
            if file_sha256(result['source'][name]) != result['source']['sha256' if name == 'path' else 'index_sha256']:
                raise ValueError('Recorded acquisition bytes changed')
        return result
    if destination.exists():
        raise ValueError('Unreceipted derivative exists: ' + str(destination))
    temporary = destination.with_name(destination.name + '.partial')
    counts = dict(original_records=0, retained_records=0)
    with _hashed_vcf_stream(source) as (reader, original):
        if source['vcf_sample'] not in reader.header.samples:
            raise ValueError('Declared VCF sample is absent')
        resolver = ContigResolver(reader.header.contigs, reference_name=source['reference_name'])
        regions = []
        for variant in variants['values']:
            resolved = resolver.resolve(variant.original_contig).resolved_name
            if resolved is None:
                continue
            start, end = normal_query_interval(variant)
            length = reader.header.contigs[resolved].length
            regions.append((resolved, max(0, start - flank), min(length, end + flank) if length else end + flank))
        index = GenomicIntervalIndex(regions)
        with pysam.VariantFile(temporary, 'wz', header=reader.header) as writer:
            for record in reader:
                counts['original_records'] += 1
                if index.overlaps(record.contig, record.start, record.stop):
                    writer.write(record)
                    counts['retained_records'] += 1
    pysam.tabix_index(str(temporary), preset='vcf')
    temporary.replace(destination)
    Path(str(temporary) + '.tbi').replace(Path(str(destination) + '.tbi'))
    selected = dict(source, path=str(destination), sha256=file_sha256(destination),
        index_path=str(destination) + '.tbi', index_sha256=file_sha256(str(destination) + '.tbi'),
        acquired_regions=index.intervals(), provenance=dict(source['provenance'], original_file=original,
            acquisition='All original records overlapping nominated intervals plus declared flanks; no genotype filtering'))
    result = dict(schema='isovar.normal_acquisition.v1', selection=source, flank=flank,
                  nomination_sha256=variants['sha256'], source=selected, counts=counts)
    write_json(receipt_path, result)
    print(source['source_id'][:12], counts, flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory', type=Path, default=Path(__file__).with_name('sid_normal_inventory.json.gz'))
    parser.add_argument('--variants', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--flank', type=int, default=500)
    parser.add_argument('--workers', type=int, default=3)
    args = parser.parse_args()
    if args.flank < 0 or not 1 <= args.workers <= 4:
        parser.error('Use a nonnegative flank and 1-4 acquisition workers')
    inventory = read_json(args.inventory)
    nominations = read_json(args.variants)
    variants = dict(sha256=file_sha256(args.variants),
        values=[Variant(row['contig'], row['start'], row['ref'], row['alt'], genome='GRCh38') for row in nominations])
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    results = {}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs = [pool.submit(acquire, source, variants, args.flank, output)
                for source in inventory['sources'] if source['selected']]
        for job in as_completed(jobs):
            result = job.result()
            results[result['source']['source_id']] = result['source']
    manifest = dict(schema='isovar.normal_manifest.v1', sample_id='Sid', reference_name='GRCh38',
        snapshot=inventory['snapshot'], background_flank=args.flank,
        variants=dict(path=str(args.variants.resolve()), sha256=variants['sha256']),
        sources=[results.get(source['source_id'], source) for source in inventory['sources']])
    write_json(output / 'normal_manifest.json.gz', manifest)


if __name__ == '__main__':
    main()
