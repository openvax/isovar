"""Replay all ALT-supported ncORF overlaps through the shared RNA/ORF API."""

import gzip
import hashlib
import json
from pathlib import Path

import orfdata
import pysam
from pyensembl import EnsemblRelease
from varcode import Variant

from isovar import ProteinSequenceCreator, ReadCollector, reconstruct_variant_orfs
from .build_evidence_pdfs import ROOT, collect


HERE = Path(__file__).parent / 'ncorf_rna_audit'


def build():
    genome, catalog = EnsemblRelease(87), orfdata.load_catalog()
    results, manifest = [], []
    for sample, source, path, scope in (
            ('COLO829', 'SRR8615617', ROOT / 'tests/data/colo829_whole_run/SRR8615617.all_sites.grch38.bam',
             'Whole-run selected-site BAM; original ALT witnesses at every selected locus'),
            ('HCC1395', 'SRR23456709', HERE / 'inputs/HCC1395.original-witnesses.grch38.bam',
             'Recovered original atlas ALT witnesses; complete-genome realignment; no allele-frequency denominator')):
        data = collect(sample)
        manifest.append(dict(sample=sample, source=source, path=str(path.relative_to(ROOT)),
                             sha256=hashlib.sha256(path.read_bytes()).hexdigest(), scope=scope))
        with pysam.AlignmentFile(path) as alignment:
            for call in data['calls']:
                annotation = data['annotations'][call['call_id']]
                if annotation['orfdata_coding_block_overlap'] != 'True' or not int(call['alt_fragments']):
                    continue
                variant = Variant(call['chrom'].removeprefix('chr'), int(call['pos']), call['ref'], call['alt'], genome)
                references = [catalog.get(ident).bind(genome.transcript_by_id(catalog.get(ident).transcript_id))
                              for ident in annotation['orfdata_orf_ids'].split(';')]
                for transcript in {r.transcript.id: r.transcript for r in references}.values():
                    if transcript.is_protein_coding and transcript.complete:
                        references.append(orfdata.TranscriptORF.from_annotation(transcript))
                evidence = ReadCollector().read_evidence_for_variant(variant, alignment)
                for coverage in (2, 1):
                    result = reconstruct_variant_orfs(
                        variant, evidence, references, sample_id=sample, source=source,
                        protein_sequence_creator=ProteinSequenceCreator(
                            min_variant_sequence_coverage=coverage, max_protein_sequences_per_variant=0))
                    if result['allele_support']['fragments'] != int(call['alt_fragments']):
                        raise ValueError('Recovered original ALT support differs from atlas: ' + call['call_id'])
                    result.update(call_id=call['call_id'], input_scope=scope,
                                  atlas_alt_fragments=int(call['alt_fragments']),
                                  pass_name='default_coverage' if coverage == 2 else 'single_read_sensitivity',
                                  igv_locus=data['loci'][call['call_id']]['igv_loci'])
                    results.append(result)
                    print(sample, call['call_id'], result['pass_name'],
                          len(result['rna_candidates']), len(result['protein_candidates']), flush=True)
    with gzip.GzipFile(filename=str(HERE / 'results.json.gz'), mode='wb', mtime=0) as stream:
        stream.write((json.dumps(results, sort_keys=True, separators=(',', ':')) + '\n').encode())
    (HERE / 'manifest.json').write_text(json.dumps(dict(
        schema='isovar.ncorf_cohort_audit.v1', inputs=manifest, catalog_metadata=catalog.metadata,
        ensembl_release=87, passes=[dict(name='default_coverage', minimum_read_coverage=2),
                                  dict(name='single_read_sensitivity', minimum_read_coverage=1)],
        interpretation='Passes reuse the same original molecules; never pool their support. '
                       'Extra ORF models and annotated CDS competitors remain separate. '
                       'Local sequence candidates are not measured translation.',
        results_sha256=hashlib.sha256((HERE / 'results.json.gz').read_bytes()).hexdigest()), indent=2) + '\n')
    return results


if __name__ == '__main__':
    build()
