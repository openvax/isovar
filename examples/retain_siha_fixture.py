"""Retain the existing SiHa regression slice, without downloading or realigning."""

import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import pysam

from isovar import reconstruct_viral_sequences


TARGET = 'isovar_viral_a22b9eba3846fcf4496e065f'
SOURCE = 'SRR12650977:SiHa-1:bounded-original-RNA'
ROOT = Path(__file__).resolve().parents[1]
FASTA = ROOT / 'tests/data/viral/hpv16-GDC.fa.gz'


def digest(value):
    return sha256(value).hexdigest()


def retain(bam_path, replay_path, output):
    bam_path, replay_path, output = Path(bam_path), Path(replay_path), Path(output)
    receipt_bytes = replay_path.read_bytes()
    receipt = json.loads(receipt_bytes)
    bam_sha = digest(bam_path.read_bytes())
    if receipt['input']['sha256'] != bam_sha:
        raise ValueError('Original BAM does not match the historical replay receipt')
    names = {row['fragment'][1] for row in receipt['support']['witnesses']
             if any(TARGET in segment['overlapping_sequences'] for segment in row['segments'])
             and not row['competing_placement_outside_targets']}
    output.mkdir(parents=True, exist_ok=False)
    with pysam.AlignmentFile(str(bam_path)) as source:
        header = source.header.to_dict()
        sam = [record.to_string() for record in source if record.query_name in names]
        (output / 'source-header.sam.gz').write_bytes(gzip.compress(str(source.header).encode(), mtime=0))
        header['SQ'] = [row for row in header['SQ'] if row['SN'] == 'HPV16']
        with pysam.AlignmentFile(str(output / 'original-rna.bam'), 'wb', header=header) as destination:
            for text in sam:
                # Parsing fails if a selected mate/alternative uses any other
                # contig. No record is silently dropped or its fields changed.
                record = pysam.AlignedSegment.fromstring(text, destination.header)
                if record.to_string() != text:
                    raise ValueError('A selected original SAM record changed')
                destination.write(record)
    (output / 'full-input-replay.json.gz').write_bytes(gzip.compress(receipt_bytes, mtime=0))
    with TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        fasta = tmp / 'reference.fa'
        fasta.write_bytes(gzip.decompress(FASTA.read_bytes()))
        with pysam.AlignmentFile(str(output / 'original-rna.bam')) as bam:
            result = reconstruct_viral_sequences(
                bam, genome_fasta=fasta, sample_id='SiHa-1', source=SOURCE,
                contigs=['HPV16'], reference_profile=None, max_discovery_reads=500,
                max_paths=32, reference_cache=tmp / 'reference-cache')
    target, = [row for row in result['sequences'] if row['id'] == TARGET]
    coding, = [row for row in target['coding']
               if any(ref['protein_id'] == 'NP_041326.1' for ref in row['reference_orfs'])]
    manifest = dict(
        source=dict(sample='GSM4785628', run='SRR12650977', project='PRJNA663801', series='GSE158033',
                    sample_name='SiHa-1', library_prefix_pairs=1000000,
                    url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4785628'),
        source_bam_sha256=bam_sha, source_result_sha256=digest(receipt_bytes),
        source_fasta=receipt['reference']['files']['genome_fasta'],
        selection='whole query-name groups with a segment overlapping the historical target; no recorded placement outside HPV16',
        selected_query_names=sorted(names), selected_fragments=len(names), selected_records=len(sam),
        sam_record_sha256=[digest(text.encode()) for text in sam], historical_target_id=TARGET,
        replay=dict(source=SOURCE, contigs=['HPV16'], reference_profile=None,
                    max_discovery_reads=500, max_paths=32, min_overlap=30),
        reference_fasta=dict(path='../hpv16-GDC.fa.gz', sha256=digest(FASTA.read_bytes()),
                             sequence_sha256=result['coding_reference']['contigs'][0]['sequence_sha256']),
        orfdata=dict(version=result['coding_reference']['version'],
                     catalogue_sha256=result['coding_reference']['contigs'][0]['reference']['catalogue_sha256'],
                     source=result['coding_reference']['contigs'][0]['reference']['source']),
        expected=dict(paths=len(result['sequences']), target_fragments=result['support']['target_fragments'],
                      target_rna_sha256=digest(target['sequence'].encode()),
                      target_rna_length=len(target['sequence']), target_protein_sha256=digest(coding['amino_acids'].encode()),
                      target_protein_length=len(coding['amino_acids']),
                      compatible_fragments=target['compatible_fragments'],
                      supporting_fragments=target['supporting_fragments'],
                      supporting_segments=target['supporting_segments']),
        files={path.name: digest(path.read_bytes()) for path in sorted(output.iterdir())})
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bam', required=True)
    parser.add_argument('--replay', required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    manifest = retain(args.bam, args.replay, args.output_dir)
    print(json.dumps(manifest['expected'], indent=2))
