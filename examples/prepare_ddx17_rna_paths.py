"""Prepare source-verified DDX17 inputs for the standard variant-rna-paths CLI.

This validates and projects the retained original local STAR alignments. It
performs no RNA assembly or protein inference; those belong to the shared CLI.
"""

import argparse
from dataclasses import asdict
from pathlib import Path

import pysam
from pyensembl import EnsemblRelease
from varcode import ContigResolver

from isovar.dna import reverse_complement_dna
from isovar.fusion import FusionReference
from isovar.sv_rna_audit import file_sha256, read_json, write_json


LOCAL_CONTIG = 'chr22_DDX17_38470001_38520000'
OFFSET = 38470000
LENGTH = 50000


def require(condition, message):
    if not condition:
        raise ValueError(message)


def references(audit, fasta, contig):
    """Validate pinned transcript sequences and their declared annotation CDS."""
    supplied = read_json(audit / 'models.json')
    genome = EnsemblRelease(supplied['ensembl_release'])
    result = []
    for row in supplied['models']:
        exons = tuple(sorted((start - 1, end) for start, end in row['exons']))
        sequence = ''.join(fasta.fetch(contig, start, end).upper() for start, end in exons)
        if row['strand'] == '-':
            sequence = reverse_complement_dna(sequence)
        require(sequence == row['sequence'], 'Transcript differs from pinned genomic sequence: ' + row['id'])
        start = end = None
        if row['complete']:
            transcript = genome.transcript_by_id(row['id'])
            require(transcript.first_start_codon_spliced_offset == row['cds_start'], 'Annotated CDS start differs')
            start, end = row['cds_start'], transcript.last_stop_codon_spliced_offset + 1
        result.append(FusionReference(row['id'], 'GRCh38', 'Ensembl:' + str(supplied['ensembl_release']),
                                     'chr22', row['strand'], exons, sequence, start, end))
    return result


def project(audit, fasta, contig, output):
    """Translate only verified coordinates, keeping all other SAM fields intact."""
    local = ''.join(line.strip() for line in (audit / 'reference.fa').read_text().splitlines() if not line.startswith('>'))
    require(local.upper() == fasta.fetch(contig, OFFSET, OFFSET + LENGTH).upper(), 'Local and target reference slices differ')
    target = output / 'DDX17.original.GRCh38.bam'
    count = 0
    with pysam.AlignmentFile(audit / 'expanded.star.bam') as source:
        require(source.references == (LOCAL_CONTIG,) and source.lengths == (LENGTH,), 'Unexpected local BAM dictionary')
        header = source.header.to_dict()
        header['SQ'] = [dict(SN='chr22', LN=fasta.get_reference_length(contig))]
        header.setdefault('CO', []).append('Sequence-verified projection +38470000 from a 50-kb local STAR reference. '
                                            'NH uniqueness remains LOCAL, not whole-genome uniqueness.')
        with pysam.AlignmentFile(target, 'wb', header=header) as dest:
            for read in source:
                require({tag for tag, _ in read.tags} <= {'NH', 'HI', 'AS', 'nM'},
                        'Unexpected tag; review coordinate-bearing metadata before projection')
                for identity, start in ((read.reference_id, read.reference_start),
                                         (read.next_reference_id, read.next_reference_start)):
                    require(identity == -1 or (identity == 0 and 0 <= start < LENGTH), 'Uninterpretable original coordinate')
                if read.reference_id >= 0:
                    read.reference_start += OFFSET
                if read.next_reference_id >= 0:
                    read.next_reference_start += OFFSET
                dest.write(read)
                count += 1
    pysam.index(str(target))
    return target, count


def prepare(audit, reference, output):
    """Pin source bytes, validate the projection, and write standard CLI inputs."""
    audit, reference, output = Path(audit).resolve(), Path(reference).resolve(), Path(output).resolve()
    manifest = read_json(audit / 'manifest.json')
    for name in ('expanded.star.bam', 'reference.fa', 'models.json'):
        require(file_sha256(audit / name) == manifest['files'][name], 'Original audit checksum mismatch: ' + name)
    require(manifest['local_reference_interval_0based'] == ['chr22', OFFSET, OFFSET + LENGTH], 'Unexpected source interval')
    output.mkdir(parents=True, exist_ok=False)
    with pysam.FastaFile(reference) as fasta:
        contig = ContigResolver(fasta.references, reference_name='GRCh38').resolve('chr22').resolved_name
        require(contig is not None, 'GRCh38 chr22 is unresolved in the supplied reference')
        models = references(audit, fasta, contig)
        bam, count = project(audit, fasta, contig, output)
    provenance = dict(source_manifest=manifest, records=count, projection_offset=OFFSET,
        scope='Recovered original SRR8615617 pairs aligned to a 50-kb DDX17 slice. '
              'NH and uniqueness remain local; source selection is not a whole-library RNA screen.',
        local_reference_matches_target_slice=True,
        files={name: dict(path=str(path), sha256=file_sha256(path)) for name, path in
               [('original_bam', audit / 'expanded.star.bam'), ('local_reference', audit / 'reference.fa'),
                ('target_reference', reference), ('projected_bam', bam), ('models', audit / 'models.json')]})
    write_json(output / 'projection.json', provenance)
    write_json(output / 'input.json', dict(reference_name='GRCh38', sample_id='COLO829',
        variant=dict(contig='chr22', start=38489856, ref='GA', alt='G'),
        references=[asdict(model) for model in models], provenance=provenance))
    files = {name: dict(path=str(path), sha256=file_sha256(path)) for name, path in
             [('original-rna', bam), ('ddx17-input', output / 'input.json')]}
    write_json(output / 'manifest.json', dict(schema='isovar.cohort_manifest.v1',
        samples=[dict(sample_id='COLO829')], files=files,
        sources=[dict(source_id='SRR8615617.recovered-local-STAR', sample_id='COLO829',
                      library_id='SRR8615617', bam='original-rna', scope=provenance['scope'], provenance=provenance)],
        variant_rna_paths=[dict(id='DDX17', title='DDX17 intronic deletion: original ALT and REF RNA',
            source_id='SRR8615617.recovered-local-STAR', sample_id='COLO829', input='ddx17-input',
            parameters=dict(max_paths=128, max_discovery_reads=10000, max_seed_segments=10000,
                            min_overlap=30, min_anchor_bases=12))]))
    return output / 'manifest.json'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original-audit', required=True)
    parser.add_argument('--genome-fasta', required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    try:
        print(prepare(args.original_audit, args.genome_fasta, args.output_dir))
    except (ValueError, OSError, KeyError) as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()
