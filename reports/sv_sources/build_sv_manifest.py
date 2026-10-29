"""Prepare a source inventory for the shared RNA audit; never reconstruct RNA here."""

import argparse
from pathlib import Path

from pyensembl import EnsemblRelease
from varcode import ContigResolver

from isovar import __version__, prepare_sv_source_record, read_sv_source
from isovar.sv_rna_audit import file_sha256, read_json, write_json
from isovar.sv_sources import SOURCE_PROFILES


def build(args):
    library_path = args.libraries.resolve()
    libraries = read_json(library_path)
    files = {key: dict(value, path=str((library_path.parent / value['path']).resolve()))
             for key, value in libraries['files'].items()}
    sources = libraries['sources']
    if not sources:
        raise ValueError('Supply at least one separate RNA library')
    records = read_sv_source(args.vcf, sha256=args.sha256, reference_name=args.reference_name,
                             profile=args.profile)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if 'source-vcf' in files:
        raise ValueError('Library file IDs must not use reserved source-vcf')
    files['source-vcf'] = dict(path=str(args.vcf.resolve()), sha256=args.sha256,
                               description='Unmodified original source SV definitions')
    genome = EnsemblRelease(args.ensembl_release)
    resolver = ContigResolver(genome.contigs(), reference_name=genome.reference_name)
    if genome.reference_name != args.reference_name:
        raise ValueError('Source and annotation assemblies must agree')
    events, problems = [], []
    for ordinal, record in enumerate(records, 1):
        event_id = f'{args.sample_id}:{ordinal}:{record.record_id}'
        event = dict(event_id=event_id, sources=[source['source_id'] for source in sources],
                     provenance=dict(file_id='source-vcf', source_profile=args.profile, source_sha256=args.sha256,
                                     isovar_version=__version__, source_interpretation=record.to_dict()))
        events.append(event)
        if record.variant is None:
            event['unassessed_reason'] = '; '.join(record.reasons)
            problems.append(dict(event_id=event_id, reason=event['unassessed_reason']))
            continue
        prepared = prepare_sv_source_record(record, genome, resolver,
            annotation=f'Ensembl {args.ensembl_release}', region_flank=args.region_flank)
        reference_id, hypothesis_id = f'event-{ordinal}-refs', f'event-{ordinal}-hypotheses'
        for key, payload in ((reference_id, prepared.reference_payload()),
                             (hypothesis_id, prepared.hypotheses)):
            if key in files:
                raise ValueError('Library file ID conflicts with prepared event: ' + key)
            path = output / (key + '.json.gz')
            write_json(path, payload)
            files[key] = dict(path=str(path), sha256=file_sha256(path))
        event['views'] = [dict(view_id=view.view_id, references=reference_id, hypotheses=hypothesis_id,
            parameters=dict(exon_anchored_clips=bool(prepared.references)),
            input=view.to_dict(prepared.reference_name)) for view in prepared.views]
        event['provenance'].update(gene_names=prepared.gene_names, nominal_igv=prepared.nominal_igv)
        print(event_id, len(prepared.references), 'models', len(prepared.hypotheses), 'hypotheses', flush=True)
    manifest = dict(schema='isovar.sv_rna_audit_manifest.v1', files=files, sources=sources, events=events,
        description=f'All original {args.sample_id} source rows; explicit source interpretation, separate libraries.',
        parameters=dict(max_records=10000, max_queries=100, max_paths=20, max_extension_segments=300,
                        assemble=True, include_regional_candidates=False, dense_support=True, exon_anchored_clips=True),
        read_collection=dict(min_mapping_quality=0, use_soft_clipped_bases=True),
        assessment=dict(min_path_match=20, min_exon_overlap=15, junction_flank=12, min_overlap=30, max_chain_reads=30))
    write_json(output / 'manifest.json', manifest)
    write_json(output / 'input_problems.json', problems)
    print(len(events), 'source rows;', len(problems), 'unassessed definitions', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vcf', type=Path, required=True)
    parser.add_argument('--sha256', required=True)
    parser.add_argument('--profile', choices=SOURCE_PROFILES, default='standard')
    parser.add_argument('--reference-name', required=True)
    parser.add_argument('--sample-id', required=True)
    parser.add_argument('--libraries', type=Path, required=True,
                        help='JSON containing the audit schema files and sources for separate RNA libraries')
    parser.add_argument('--ensembl-release', type=int, default=87)
    parser.add_argument('--region-flank', type=int, default=10000)
    parser.add_argument('--output', type=Path, required=True)
    build(parser.parse_args())
