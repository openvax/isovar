"""Pin gene names and one-based IGV loci without changing scientific calls."""

import csv
from hashlib import sha256
import json
from pathlib import Path

from pyensembl import EnsemblRelease
from reports.hcc1395_colo829.build_evidence_pdfs import read_gzip, read_tsv

BASE = Path(__file__).resolve().parent


def main():
    annotation = EnsemblRelease(87)
    output = []
    for sample, prefix in [('COLO829', 'colo829'), ('HCC1395', 'hcc1395_direct_rna')]:
        protein_rows = read_gzip(BASE / f'{prefix}_atlas_proteins.tsv.gz')
        protein_calls = {r['call_id'] for r in protein_rows}
        protein_genes = {}
        for row in protein_rows:
            names = protein_genes.setdefault(row['call_id'], set())
            for tid in row['compatible_transcript_ids'].split(';'):
                if tid:
                    names.add(annotation.transcript_by_id(tid).gene_name)
        for row in read_gzip(BASE / f'{prefix}_atlas_calls.tsv.gz'):
            if row['call_id'] not in protein_calls and not int(row['alt_fragments']):
                continue
            start, end = int(row['pos']), int(row['pos']) + len(row['ref']) - 1
            genes = annotation.gene_names_at_locus(row['chrom'].removeprefix('chr'), start, end)
            output.append(dict(kind='small_variant', sample=sample, dna_source=row['dna_source'],
                               call_id=row['call_id'], genes=';'.join(sorted(set(genes))),
                               protein_genes=';'.join(sorted(protein_genes.get(row['call_id'], []))),
                               endpoint1_genes='', endpoint2_genes='',
                               igv_loci=f"{row['chrom']}:{start}-{end}",
                               genome_build='GRCh38', annotation='Ensembl 87'))
        filename = 'colo829_sv_rna_screen.tsv' if sample == 'COLO829' else 'hcc1395_direct_rna_sv_screen.tsv'
        for row in read_tsv(BASE / filename):
            genes, loci = [], []
            for number in (1, 2):
                chrom, pos = row[f'chrom{number}'], row[f'pos{number}']
                if chrom and chrom != 'NA' and pos.isdigit() and int(pos) > 0:
                    names = annotation.gene_names_at_locus(chrom.removeprefix('chr'), int(pos))
                    genes.append(sorted(set(names)))
                    loci.append(f'{chrom}:{pos}-{pos}')
                else:
                    genes.append([])
                    loci.append('unresolved')
            output.append(dict(kind='sv', sample=sample, dna_source=row['dna_source'], call_id=row['call_id'],
                               genes=';'.join(sorted(set(genes[0] + genes[1]))),
                               protein_genes='',
                               endpoint1_genes=';'.join(genes[0]), endpoint2_genes=';'.join(genes[1]),
                               igv_loci=' '.join(loci), genome_build='GRCh38', annotation='Ensembl 87'))
    with (BASE / 'report_locus_annotations.tsv').open('w') as stream:
        writer = csv.DictWriter(stream, output[0].keys(), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(output)
    transcript = annotation.transcript_by_id('ENST00000304494')
    (BASE / 'focused_proteins/cdkn2a_reference.json').write_text(json.dumps(dict(
        transcript_id=transcript.id, protein_id=transcript.protein_id, gene=transcript.gene_name,
        annotation='Ensembl 87', genome_build='GRCh38', sequence=transcript.protein_sequence,
        protein_sha256=sha256(transcript.protein_sequence.encode()).hexdigest()), indent=2) + '\n')
    print(f'Pinned {len(output)} source-specific locus annotations.')


if __name__ == '__main__':
    main()
