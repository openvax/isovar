"""Pin the complete COLO829 SV inventory for the shared library audit.

This prepares inputs, not RNA evidence. Varcode owns VCF adjacency semantics
and consequence prediction; the runner needs no annotation download.
"""

from collections import Counter
import json
from pathlib import Path

from pyensembl import EnsemblRelease
import varcode

from isovar.sv_rna_audit import file_sha256, write_json
from isovar.sv_sources import read_sv_source
from isovar.sv_source_preparation import prepare_sv_source_record

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
DATA = ROOT / 'tests/data/colo829_whole_run'
OUTPUT = HERE / 'unified_sv_audit'
SV3 = 'PaulinPaz_2024_v2:COLO829_SOMATIC_SV3_DEL'
SOURCE_URL = 'https://github.com/PacificBiosciences/DevNet/wiki/Melanoma--Cancer-Cell-Line-Iso-Seq-Data'


def build(output=OUTPUT):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output / 'inputs').mkdir(exist_ok=True)
    files = {}

    def pin(identifier, path, description):
        import os
        files[identifier] = dict(path=os.path.relpath(path, output), sha256=file_sha256(path), description=description)
        return identifier

    for key, name in [('all_sites_bam', 'SRR8615617.all_sites.grch38.bam'),
                      ('all_sites_index', 'SRR8615617.all_sites.grch38.bam.bai'),
                      ('cnih3_bam', 'SRR8615617.cnih3_source_candidates.grch38.bam'),
                      ('cnih3_index', 'SRR8615617.cnih3_source_candidates.grch38.bam.bai'),
                      ('sv3_raw', 'SRR8615617.sv3_exact_window_spots.fastq.gz')]:
        pin(key, DATA / name, 'Unmodified pinned original-read evidence; see tests/data/colo829_whole_run/README.md')
    sources = []
    for source_id, bam, selection in [
            ('SRR8615617-selected-sites', 'all_sites_bam', 'Full SRA run aligned; only nominated sites +/-500bp retained'),
            ('SRR8615617-CNIH3', 'cnih3_bam', 'Independent CNIH3-wide subset of the SAME library; not independent replication')]:
        sources.append(dict(source_id=source_id, sample_id='COLO829', library_id='SRR8615617',
                            culture='SRR8615617-separately-grown-from-benchmark-DNA', role='tumor', bam=bam,
                            provenance=dict(accession='SRR8615617', selection=selection,
                                            limitations=['selected-site evidence cannot exclude unsampled distant exons or placements',
                                                         'different processed subsets are not independent libraries'])))
    sources[1].update(raw_fastq='sv3_raw', fastq_name_format='sra_dot_mate')
    for sample, role in [('COLO829T', 'tumor'), ('COLO829BL', 'normal')]:
        sources.append(dict(source_id='PacBio2019-' + sample, sample_id=sample, library_id='PacBio2019-' + sample,
                            culture='PacBio2019-separate-ATCC-culture-' + sample, role=role,
                            provenance=dict(source_url=SOURCE_URL, issue='https://github.com/openvax/isovar/issues/574'),
                            unavailable_reason='Original CCS/FLNC RNA unavailable: documented PacBio host timed out 2026-10-26; alternate host 404. Polished representatives are not original molecules. See #574.'))
    genome = EnsemblRelease(87)
    annotation_contigs = varcode.ContigResolver(genome.contigs(), reference_name=genome.reference_name)
    events, problems = [], []
    inventories = [('PaulinPaz_2024_v2', 'colo829-sv-grch38.vcf'),
                   ('ValleInclan_2022_lifted', 'colo829-sv-2022-hg38lifted.vcf')]
    for inventory, name in inventories:
        pin(inventory, DATA / name, 'Original unmodified source VCF; distinct benchmark inventory')
        profile = 'colo829_paulin_paz_v2' if inventory == 'PaulinPaz_2024_v2' else 'standard'
        records = read_sv_source(DATA / name, sha256=files[inventory]['sha256'],
                                 reference_name='GRCh38', profile=profile)
        for record in records:
            event_id = inventory + ':' + record.record_id
            event = dict(event_id=event_id, sources=['SRR8615617-selected-sites'],
                         provenance=dict(inventory=inventory, file_id=inventory, line_number=record.line_number,
                                         original_vcf_row=record.original_row, assembly='GRCh38',
                                         source_interpretation=record.to_dict(),
                                         varcode_version=varcode.__version__))
            events.append(event)
            variant = record.variant
            if variant is None:
                event['unassessed_reason'] = '; '.join(record.reasons)
                problems.append(dict(event_id=event_id, reason=event['unassessed_reason']))
                continue
            prepared = prepare_sv_source_record(record, genome, annotation_contigs, annotation='Ensembl 87')
            models = prepared.reference_payload()
            key = 'event-%03d' % len(events)
            model_path = output / 'inputs' / (key + '.references.json.gz')
            write_json(model_path, models)
            ref = pin(key + '-refs', model_path, models['selection'])
            hypotheses = prepared.hypotheses
            hypothesis_path = output / 'inputs' / (key + '.varcode.json.gz')
            write_json(hypothesis_path, hypotheses)
            hyp = pin(key + '-hypotheses', hypothesis_path, 'All Varcode hypotheses, including null materialized cDNA')
            event['views'] = [dict(view_id=view.view_id, references=ref, hypotheses=hyp,
                parameters=dict(exon_anchored_clips=bool(prepared.references)),
                input=view.to_dict(prepared.reference_name)) for view in prepared.views]
            event['provenance'].update(gene_names=prepared.gene_names, nominal_igv=prepared.nominal_igv)
            if event_id == SV3:
                event['sources'] = [s['source_id'] for s in sources]
                # Match the established SV3 discovery budget. The exact same
                # library assessment/reconciliation runs on every retained path.
                for view in event['views']:
                    view['parameters'].update(max_records=20000, max_queries=200, max_paths=80, max_extension_segments=1200)
            print(event_id, len(models['models']), 'models', len(hypotheses), 'hypotheses', flush=True)
    manifest = dict(schema='isovar.sv_rna_audit_manifest.v1', description='All 188 COLO829 source VCF records; separate views, libraries and cultures; see README.md',
                    files=files, sources=sources, events=events,
                    parameters=dict(max_records=10000, max_queries=100, max_paths=20, max_extension_segments=300,
                                    assemble=True, include_regional_candidates=False, dense_support=True,
                                    exon_anchored_clips=True),
                    read_collection=dict(min_mapping_quality=0, use_soft_clipped_bases=True),
                    assessment=dict(min_path_match=20, min_exon_overlap=15, junction_flank=12, min_overlap=30, max_chain_reads=30))
    write_json(output / 'manifest.json', manifest)
    write_json(output / 'input_problems.json', problems)
    print(json.dumps(dict(events=len(events), by_inventory=Counter(e['provenance']['inventory'] for e in events),
                          unassessed_definitions=len(problems)), sort_keys=True))
    return manifest


if __name__ == '__main__':
    build()
