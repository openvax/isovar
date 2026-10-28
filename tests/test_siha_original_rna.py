"""Offline original SiHa RNA replay, separate from its full-input receipt."""

from collections import Counter
import gzip
from hashlib import sha256
import json
from pathlib import Path

import pysam

from isovar.cli.commands import run


DATA = Path(__file__).parent / 'data/viral/siha'


def test_original_records_reference_pins_and_full_input_receipt():
    manifest = json.loads((DATA / 'manifest.json').read_text())
    for filename, digest in manifest['files'].items():
        assert sha256((DATA / filename).read_bytes()).hexdigest() == digest
    reference = manifest['reference_fasta']
    assert sha256((DATA / reference['path']).read_bytes()).hexdigest() == reference['sha256']
    with pysam.AlignmentFile(str(DATA / 'original-rna.bam')) as bam:
        rows = list(bam)
        assert bam.references == ('HPV16',)
        assert [sha256(row.to_string().encode()).hexdigest() for row in rows] == manifest['sam_record_sha256']
        assert Counter(row.query_name for row in rows) == {name: 2 for name in manifest['selected_query_names']}
        assert all({row.flag & 192 for row in rows if row.query_name == name} == {64, 128}
                   for name in manifest['selected_query_names'])
        assert all(row.query_qualities is not None and row.get_tag('RG') == 'SRR12650977' for row in rows)
        source_header = pysam.AlignmentHeader.from_text(gzip.decompress((DATA / 'source-header.sam.gz').read_bytes()).decode())
        assert set(bam.references) < set(source_header.references)
        assert bam.header.to_dict()['PG'] == source_header.to_dict()['PG']
    raw = gzip.decompress((DATA / 'full-input-replay.json.gz').read_bytes())
    assert sha256(raw).hexdigest() == manifest['source_result_sha256']
    full = json.loads(raw)
    assert full['input']['sha256'] == manifest['source_bam_sha256']
    assert full['reference']['files']['genome_fasta'] == manifest['source_fasta']
    assert len(full['sequences']) == 26 and full['support']['target_fragments'] == 71
    assert sum(w['competing_placement_outside_targets'] for w in full['support']['witnesses']
               if manifest['historical_target_id'] in w['compatible_sequences']) == 3


def test_cli_replays_original_slice_without_downloads_or_reference_padding(tmp_path):
    manifest = json.loads((DATA / 'manifest.json').read_text())
    fasta = tmp_path / 'reference.fa'
    fasta.write_bytes(gzip.decompress((DATA / manifest['reference_fasta']['path']).read_bytes()))
    output = tmp_path / 'result'
    run(['viral', '--bam', str(DATA / 'original-rna.bam'), '--genome-fasta', str(fasta),
         '--sample-id', 'SiHa-1', '--alignment-reference', 'custom', '--contig', 'HPV16',
         '--max-discovery-reads', '500', '--max-paths', '32',
         '--reference-cache', str(tmp_path / 'cache'), '--output-dir', str(output)])
    result = json.loads((output / 'viral.json').read_text())
    expected = manifest['expected']
    assert len(result['sequences']) == expected['paths'] == 7
    assert result['support']['target_fragments'] == expected['target_fragments'] == 23
    assert result['input']['sha256'] == manifest['files']['original-rna.bam']
    reference = result['coding_reference']['contigs'][0]
    assert reference['reference']['catalogue_sha256'] == manifest['orfdata']['catalogue_sha256']
    assert reference['reference']['source'] == manifest['orfdata']['source']
    assert reference['sequence_sha256'] == manifest['reference_fasta']['sequence_sha256']
    # Original header has no M5. The manifest pins its source reference; the
    # custom replay must not pretend its BAM carries a sequence checksum.
    assert result['reference']['bam_sequence_identity_unverified'] == ['HPV16']
    target, = [row for row in result['sequences'] if row['id'] == manifest['historical_target_id']]
    assert sha256(target['sequence'].encode()).hexdigest() == expected['target_rna_sha256']
    for field in ('compatible_fragments', 'supporting_fragments', 'supporting_segments'):
        assert target[field] == expected[field]
    coding, = [row for row in target['coding']
               if any(ref['protein_id'] == 'NP_041326.1' for ref in row['reference_orfs'])]
    assert sha256(coding['amino_acids'].encode()).hexdigest() == expected['target_protein_sha256']
    assert len(coding['amino_acids']) == 98 and coding['amino_acids'][27] == 'F'
    ref, = coding['reference_orfs']
    assert ref['reference_amino_acids'][27] == 'L' and not ref['protein_matches_reference']
    assert coding['status'] == 'complete_CDS' and coding['all_bases_observed']
    assert not coding['whole_cds_segments'] and not target['whole_path_segments']
    assert target['rna_strand'] == 'unknown' and target['transcript_completeness'] == 'unknown'
    assert coding['evidence'] == 'reference_anchored_prediction_from_observed_RNA'
    assert (output / 'proteins.fa').is_file() and (output / 'sequences.fa').is_file()
