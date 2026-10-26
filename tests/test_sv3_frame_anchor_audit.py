"""Original SV3 exon links remain distinct from an event-linked frame."""

from isovar.sv_rna_assessment import aligned_exon_clip_contains
from reports.hcc1395_colo829.build_sv3_frame_anchor_audit import build


def test_source_exon_clips_do_not_phase_the_sv3_junction():
    runs = build()
    value = next(r['assessment'] for r in runs if r['source_id'] == 'SRR8615617-CNIH3' and r['view_id'] == '1-forward')
    assert value['records_scanned'] == 52212
    assert value['fragments_scanned'] == 26976
    assert value['coding_exon_fragments'] == 673
    path = next(p for p in value['paths'] if p['path_id'] == '653d5d9c015d3e590fb149b3')
    assert len({tuple(link['fragment_id']) for link in path['exon_links']}) == 9
    assert len({tuple(link['fragment_id']) for link in path['exon_links']
                if any(m['path_interval'][1] - m['path_interval'][0] >= 30 for m in link['matches'])}) == 6
    join = next(j for j in path['junctions'] if j['relation'] == 'breakpoint_junction')
    assert not join['exon_linked_fragments']
    assert all(match['path_interval'][0] - join['query_interval'][1] == 280
               and match['same_end_exon_clip_contains_match'] == link['same_end_exon_alignment']
               for link in path['exon_links'] for match in link['matches'])
    assert any(link['opposite_mate_exon_alignment'] for link in path['exon_links'])
    assert value['coding_reference_starts'] == [dict(transcript_id='ENST00000272133', contig='chr1',
                                                   strand='+', position_0based=224617174)]
    assert not path['proteins']


def test_exon_clip_link_requires_the_exon_alignments_own_clipped_bases():
    record = dict(query_sequence='AACCGGTTGATCAGTA', query_aligned_interval=[0, 8])
    assert aligned_exon_clip_contains(record, 'TTGATCAG')
    assert not aligned_exon_clip_contains(record, 'AACCGGTT')
    assert not aligned_exon_clip_contains(dict(record, query_aligned_interval=[0, 16]), 'TTGATCAG')


def test_original_sv3_witnesses_reproduce_the_library_assessment():
    import json
    import pysam
    from isovar.sv_rna_assessment import assess_sv_rna_paths
    from isovar.sv_rna_audit import raw_fastq_reads
    from reports.hcc1395_colo829.build_unified_sv_manifest import DATA
    snapshot = next(r for r in build() if r['source_id'] == 'SRR8615617-CNIH3' and r['view_id'] == '1-forward')
    identifier = '653d5d9c015d3e590fb149b3'
    original = snapshot['reconstruction']
    path = next(p for p in original['paths'] if p['path_id'] == identifier)
    expected = next(p for p in snapshot['assessment']['paths'] if p['path_id'] == identifier)
    raw = tuple(raw_fastq_reads(DATA / 'SRR8615617.sv3_exact_window_spots.fastq.gz'))
    with pysam.AlignmentFile(str(DATA / 'SRR8615617.cnih3_source_candidates.grch38.bam')) as bam:
        replay = assess_sv_rna_paths(dict(original, paths=[path]), bam,
                                    snapshot['assessment']['reference_models'], raw_reads=raw,
                                    **snapshot['assessment']['parameters'])
    # Compare the public JSON representation: in-memory fragment tuples become
    # arrays on disk, without changing fragment identity or support.
    replay = json.loads(json.dumps(replay))
    assert replay['paths'] == [expected]
    assert {key: snapshot['assessment']['exon_alignments'][key] for key in replay['exon_alignments']} == replay['exon_alignments']
