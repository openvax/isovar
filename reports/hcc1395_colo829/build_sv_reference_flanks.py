"""Pin small GRCh38 genomic controls beside the original-read SV fixtures."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import pysam

from isovar.dna import reverse_complement_dna


HERE = Path(__file__).parent
OUTPUT = HERE.parents[1] / 'tests/data/hcc1395_sv_junction/sv_reference_flanks.json'
SOURCE_URL = ('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/'
              'GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/'
              'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz')


def _arm(fasta, placement, offset, length):
    contig, base, strand = placement
    positions = [base + (offset + i) * (1 if strand == '+' else -1) for i in range(length)]
    if min(positions) < 0:
        raise ValueError('Reference flank outside chromosome')
    sequence = fasta.fetch(contig, min(positions), max(positions) + 1).upper()
    return sequence if strand == '+' else reverse_complement_dna(sequence)


def build(fasta_path, audit_path=HERE / 'sv_junction_audit.tsv'):
    """Extract 20 reference bases on either side of each observed RNA join."""
    with Path(fasta_path).open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    rows = {}
    with pysam.FastaFile(str(fasta_path)) as fasta, Path(audit_path).open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            if not row['junction_left']:
                continue
            left = row['junction_left'].split(':')
            right = row['junction_right'].split(':')
            left, right = (left[0], int(left[1]), left[2]), (right[0], int(right[1]), right[2])
            rows[row['call_id']] = dict(
                left=list(left), right=list(right),
                donor_before=_arm(fasta, left, -19, 20),
                donor_after=_arm(fasta, left, 1, 20),
                acceptor_before=_arm(fasta, right, -20, 20),
                acceptor_after=_arm(fasta, right, 0, 20))
    return dict(schema='isovar.sv_reference_flanks.v1', assembly='GRCh38',
                source_url=SOURCE_URL, source_fasta_sha256=digest, flank_bases=20, events=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-fasta', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    args.output.write_text(json.dumps(build(args.reference_fasta), sort_keys=True, indent=2) + '\n')


if __name__ == '__main__':
    main()
