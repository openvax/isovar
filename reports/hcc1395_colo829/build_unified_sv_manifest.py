"""Pin the complete COLO829 SV inventory for the shared library audit.

This prepares inputs, not RNA evidence. Varcode owns VCF adjacency semantics
and consequence prediction; the runner needs no annotation download.
"""

from collections import Counter
from dataclasses import asdict
import json
from pathlib import Path

from pyensembl import EnsemblRelease
import varcode
from varcode.sv_allele_parser import parse_symbolic_alt

from isovar.fusion import FusionReference
from isovar.sv_rna_audit import file_sha256, write_json
from isovar.sv_rna_reconciliation import hypotheses_from_varcode

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
DATA = ROOT / 'tests/data/colo829_whole_run'
OUTPUT = HERE / 'unified_sv_audit'
SV3 = 'PaulinPaz_2024_v2:COLO829_SOMATIC_SV3_DEL'
SOURCE_URL = 'https://github.com/PacificBiosciences/DevNet/wiki/Melanoma--Cancer-Cell-Line-Iso-Seq-Data'


def reference_models(genome, junctions):
    transcripts = {}
    for junction in junctions:
        for end in junction:
            for gene in genome.genes_at_locus(end.contig.removeprefix('chr'), max(1, end.position - 10000), end.position + 10000):
                transcripts.update((t.id, t) for t in gene.transcripts)
    models, limitations = [], []
    for transcript in sorted(transcripts.values(), key=lambda t: t.id):
        sequence = transcript.sequence
        if not sequence or set(sequence.upper()) - set('ACGT'):
            limitations.append(transcript.id + ': missing_or_ambiguous_cdna')
            continue
        fields = dict(transcript_id=transcript.id, reference_name='GRCh38', annotation='Ensembl 87',
                      contig='chr' + transcript.contig.removeprefix('chr'), strand=transcript.strand,
                      exons=sorted((e.start - 1, e.end) for e in transcript.exons), sequence=sequence.upper())
        complete = (transcript.contains_start_codon and transcript.contains_stop_codon
                    and transcript.start_codon_complete and transcript.stop_codon_complete)
        if complete:
            try:
                reference = FusionReference(**fields, cds_start=transcript.first_start_codon_spliced_offset,
                                            cds_end=transcript.last_stop_codon_spliced_offset + 1)
            except ValueError as error:
                limitations.append(transcript.id + ': CDS validation: ' + str(error))
                reference = FusionReference(**fields)
        else:
            reference = FusionReference(**fields)
        models.append(asdict(reference))
    if not models:
        limitations.append('No reference transcripts in breakpoint windows; exon-anchored clip seeding unavailable')
    return dict(models=models, limitations=limitations, annotation='Ensembl 87', reference_name='GRCh38',
                selection='all transcripts of genes overlapping a 10kb breakpoint window')


def breakpoint(end, donor):
    """Convert Varcode retained 1-based bases to oriented interbase coordinates."""
    left = end.keeps == 'left'
    return dict(contig='chr' + end.contig.removeprefix('chr'), position=end.position if left else end.position - 1,
                strand=('+' if left else '-') if donor else ('-' if left else '+'))


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
    events, problems = [], []
    inventories = [('PaulinPaz_2024_v2', 'colo829-sv-grch38.vcf'),
                   ('ValleInclan_2022_lifted', 'colo829-sv-2022-hg38lifted.vcf')]
    for inventory, name in inventories:
        pin(inventory, DATA / name, 'Original unmodified source VCF; distinct benchmark inventory')
        for number, line in enumerate((DATA / name).read_text().splitlines(), 1):
            if line.startswith('#'):
                continue
            cols = line.split('\t')
            event_id = inventory + ':' + cols[2]
            event = dict(event_id=event_id, sources=['SRR8615617-selected-sites'],
                         provenance=dict(inventory=inventory, file_id=inventory, line_number=number,
                                         original_vcf_row=line, assembly='GRCh38', varcode_version=varcode.__version__))
            events.append(event)
            info = dict(part.split('=', 1) if '=' in part else (part, True) for part in cols[7].split(';'))
            for key in ('CIPOS', 'CIEND'):
                if key in info:
                    info[key] = [int(v) for v in info[key].split(',')]
            try:
                variant = parse_symbolic_alt(cols[0].removeprefix('chr'), int(cols[1]), cols[3], cols[4], info,
                                             genome=genome, convert_ucsc_contig_names=True)
                junctions = variant.junctions
            except ValueError as error:
                event['unassessed_reason'] = 'Original event definition cannot be decoded: ' + str(error)
                problems.append(dict(event_id=event_id, reason=event['unassessed_reason']))
                continue
            if not junctions or any(end.keeps not in ('left', 'right') for pair in junctions for end in pair):
                event['unassessed_reason'] = 'Nominated event has no resolved oriented adjacency; no orientation or inserted sequence is invented'
                problems.append(dict(event_id=event_id, reason=event['unassessed_reason']))
                continue
            models = reference_models(genome, junctions)
            key = 'event-%03d' % len(events)
            model_path = output / 'inputs' / (key + '.references.json.gz')
            write_json(model_path, models)
            ref = pin(key + '-refs', model_path, models['selection'])
            hypotheses = hypotheses_from_varcode(variant.effects())
            hypothesis_path = output / 'inputs' / (key + '.varcode.json.gz')
            write_json(hypothesis_path, hypotheses)
            hyp = pin(key + '-hypotheses', hypothesis_path, 'All Varcode hypotheses, including null materialized cDNA')
            event['views'] = []
            for index, (first, second) in enumerate(junctions):
                for orientation, (left, right) in [('forward', (first, second)), ('reverse', (second, first))]:
                    regions = [['chr' + end.contig.removeprefix('chr'), max(0, end.position - 501), end.position + 500]
                               for end in (left, right)]
                    event['views'].append(dict(view_id=f'{index + 1}-{orientation}', references=ref, hypotheses=hyp,
                        parameters=dict(exon_anchored_clips=bool(models['models'])),
                        input=dict(reference_name='GRCh38', donor=breakpoint(left, True),
                                   acceptor=breakpoint(right, False), regions=regions)))
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
