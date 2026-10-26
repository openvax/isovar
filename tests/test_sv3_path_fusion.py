"""All-path SV3 audit retains alternatives and the missing normal/long reads."""

from reports.hcc1395_colo829.build_sv3_path_fusion import build


def test_source_reads_tile_exploratory_orfs_without_resolving_protein():
    runs = build()
    missing = [r for r in runs if r['status'] == 'unassessed']
    assert {r['source']['role'] for r in missing} == {'tumor', 'normal'}
    value = next(r['assessment'] for r in runs if r['source_id'] == 'SRR8615617-CNIH3' and r['view_id'] == '1-forward')
    assert len(value['paths']) > 1
    path = next(p for p in value['paths'] if p['path_id'] == '653d5d9c015d3e590fb149b3')
    join = next(j for j in path['junctions'] if j['relation'] == 'breakpoint_junction')
    assert join['direct_support']['fragments'] == 3
    assert len(join['exact_sequence_support']['complete_fragments']) == 5
    assert not path['exact_exon_anchors'] and path['frame_status'] == 'no_gene_anchor'
    assert not path['proteins'] and 'path_limit' in value['limitations']
    candidates = path['exploratory_orfs']
    assert sorted({len(c['candidate']['amino_acids']) for c in candidates}) == [37, 162]
    assert all(c['sequence_support']['overlap_chain'] for c in candidates)
    assert not any(c['initiation_observed'] or c['translation_observed'] for c in candidates)
