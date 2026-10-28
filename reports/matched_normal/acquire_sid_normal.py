"""Acquire pinned regional normal VCF derivatives without retaining full downloads.

Run as a module from the repository root. The reviewed inventory selects complete
public source IDs, exact VCF columns and policies. It does not infer specimen
identity from filenames. Original HTTP bytes are hashed into a temporary local file before HTSlib
parsing; all records overlapping the declared regions are retained, including
no-calls and rejected genotypes. Scientific interpretation belongs to the shared
``isovar normal-evidence`` command.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import Request, urlopen

import pysam
from varcode import Variant

from isovar.genomic_intervals import GenomicIntervalIndex
from isovar.normal_regions import normal_background_regions, resolve_normal_regions
from isovar.rna_evidence import canonical_json
from isovar.sv_rna_audit import file_sha256, read_json, write_json


@contextmanager
def _hashed_vcf_stream(source):
    """Hash a complete public download, parse locally, then discard raw bytes.

    Separating download from parsing avoids leaving a live HTTP response idle
    while regional indexing/filtering competes for CPU. The temporary original
    is removed on success or failure; only verified regional derivatives persist.
    """
    provenance = source['provenance']
    request = Request(provenance['url'], headers={'If-Match': provenance['http_identity']['ETag']})
    with TemporaryDirectory(prefix='isovar-public-normal-') as directory:
        original = Path(directory) / 'original.vcf.gz'
        digest, count = sha256(), 0
        with urlopen(request, timeout=120) as response, original.open('wb') as destination:
            if response.headers.get('ETag') != provenance['http_identity']['ETag']:
                raise ValueError('Source HTTP identity changed: ' + source['source_id'])
            while block := response.read(1024 * 1024):
                digest.update(block)
                count += len(block)
                destination.write(block)
        if count != provenance['source_size']:
            raise ValueError('Source size changed or download incomplete: ' + source['source_id'])
        with pysam.VariantFile(original) as reader:
            yield reader, dict(bytes=count, sha256=digest.hexdigest())


def _acquisition_identity(source):
    """Bind bytes and interval interpretation, independently of genotype policy."""
    metadata = source['provenance']
    return dict(source_id=source['source_id'], reference_name=source['reference_name'],
                url=metadata['url'], source_size=metadata['source_size'], http_identity=metadata['http_identity'],
                header_sha256=metadata.get('header_sha256'))


def acquire(source, variants, flank, output, *, additional_regions=()):
    """Select record intervals with one resolver and one merged interval index."""
    additional_regions = GenomicIntervalIndex(additional_regions).intervals()
    additional_regions_sha256 = sha256(canonical_json(additional_regions).encode()).hexdigest() if additional_regions else None
    destination = output / (source['source_id'] + '.vcf.gz')
    receipt_path = output / (source['source_id'] + '.acquisition.json.gz')
    if receipt_path.exists():
        result = read_json(receipt_path)
        if (canonical_json(_acquisition_identity(result['selection'])) != canonical_json(_acquisition_identity(source)) or result['flank'] != flank
                or result['nomination_sha256'] != variants['sha256']
                or result.get('additional_regions_sha256') != additional_regions_sha256):
            raise ValueError('Acquisition scope changed; use a new directory')
        for name in ('path', 'index_path'):
            if file_sha256(result['source'][name]) != result['source']['sha256' if name == 'path' else 'index_sha256']:
                raise ValueError('Recorded acquisition bytes changed')
        previous = result['source']
        # Acquisition retains every record and therefore does not depend on
        # confidence thresholds or sample assignment. A new assessment can use
        # the verified bytes with a new declared policy; the acquisition receipt
        # itself remains unchanged.
        updated = dict(source, **{key: previous[key] for key in (
            'path', 'sha256', 'index_path', 'index_sha256', 'acquired_regions')},
            provenance=dict(source['provenance'], original_file=previous['provenance']['original_file'],
                            acquisition=previous['provenance']['acquisition']))
        return dict(result, source=updated)
    if destination.exists():
        raise ValueError('Unreceipted derivative exists: ' + str(destination))
    print('Acquiring', source['source_id'][:12], flush=True)
    temporary = destination.with_name(destination.name + '.partial')
    counts = dict(original_records=0, retained_records=0)
    # Build the potentially large union before opening HTTP. Holding a live
    # response idle while merging millions of exon/ncORF intervals can cause
    # transport timeouts under CPU contention.
    requested = normal_background_regions(variants['values'], background_flank=flank,
                                           additional_regions=additional_regions)
    with _hashed_vcf_stream(source) as (reader, original):
        if source['vcf_sample'] not in reader.header.samples:
            raise ValueError('Declared VCF sample is absent')
        scope = resolve_normal_regions(requested, reader.header, reference_name=source['reference_name'])
        index = GenomicIntervalIndex(scope.regions)
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
            acquisition='All original records overlapping nominated intervals, declared flanks and explicit background regions; no genotype filtering'))
    result = dict(schema='isovar.normal_acquisition.v1', selection=source, flank=flank,
                  nomination_sha256=variants['sha256'], additional_regions_sha256=additional_regions_sha256,
                  unresolved_contigs=list(scope.unresolved_contigs), outside_contig_regions=list(scope.outside_contig_regions),
                  clipped_regions=list(scope.clipped_regions), source=selected, counts=counts)
    write_json(receipt_path, result)
    print(source['source_id'][:12], counts, flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory', type=Path, default=Path(__file__).with_name('sid_normal_inventory.json.gz'))
    parser.add_argument('--variants', type=Path, required=True)
    parser.add_argument('--additional-background-regions', type=Path,
                        help='JSON/JSON.gz array of [contig, zero-based start, half-open end] intervals')
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
    additional_regions = read_json(args.additional_background_regions) if args.additional_background_regions else ()
    additional_digest = file_sha256(args.additional_background_regions) if args.additional_background_regions else None
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    results = {}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs = [pool.submit(acquire, source, variants, args.flank, output, additional_regions=additional_regions)
                for source in inventory['sources'] if source['selected']]
        for job in as_completed(jobs):
            result = job.result()
            results[result['source']['source_id']] = result['source']
    manifest = dict(schema='isovar.normal_manifest.v1', sample_id='Sid', reference_name='GRCh38',
        snapshot=inventory['snapshot'], background_flank=args.flank,
        variants=dict(path=str(args.variants.resolve()), sha256=variants['sha256']),
        sources=[results.get(source['source_id'], source) for source in inventory['sources']])
    if args.additional_background_regions:
        manifest['additional_background_regions'] = dict(path=str(args.additional_background_regions.resolve()), sha256=additional_digest)
    digest = sha256(canonical_json(manifest).encode()).hexdigest()
    manifest_path = output / ('normal_manifest.' + digest[:16] + '.json.gz')
    if manifest_path.exists():
        if canonical_json(read_json(manifest_path)) != canonical_json(manifest):
            raise ValueError('Existing normal manifest content differs')
    else:
        write_json(manifest_path, manifest)
    print(manifest_path, flush=True)


if __name__ == '__main__':
    main()
