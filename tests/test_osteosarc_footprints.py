"""Pinned Sid RNA paths and explicit, non-conflated indel/deletion evidence."""

import pysam

from examples.osteosarc_footprint_figures import load, fusion_inputs
from isovar.fusion import fusion_from_dict, reconstruct_fusion
from tests.data.fusions.build_osteosarc import extract
from tests.data.osteosarc.expansion.references import translate
from tests.data.osteosarc.figure_comparisons.footprints import deletion_footprint, fusion_footprint, split_deletion_paths
from tests.test_chimeric_phasing import record, link


def test_all_nine_candidates_are_present_with_independent_indel_translation_checks():
    data=load()
    assert {e['gene'] for e in data['indels']}=={'GTF3C5','RNF213','GLIS3','KTN1'}
    assert {e['gene'] for e in data['deletions']}=={'DLG5','AFF3','KEAP1'}
    assert len(data['fusions'])==2
    for entry in data['indels']:
        assert len(entry['products'])==4
        for product in entry['products']:
            visualization=product['visualization']
            for mode in visualization['modes']:
                p=mode['protein']
                if p:
                    orf=p['witness']['orf']
                    actual,_=translate(orf['cdna'][orf['codon_offset']:])
                    assert actual[:len(p['amino_acids'])]==p['amino_acids']
            assert all(v['result'] is None or (v['result']['validation_status']=='ok' and v['result']['matches_expected'])
                       for v in visualization['provenance']['independent_validation'])


def test_real_sid_fusion_paths_and_junction_windows_recount_from_actual_records():
    data=load()
    header=pysam.AlignmentHeader.from_references(['chr'+str(i) for i in range(1,23)],[300000000]*22)
    for entry in data['fusions']:
        event=(entry['name'],*entry['breakpoints'])
        for product in entry['products']:
            records=[pysam.AlignedSegment.fromstring(s,header) for p in product['paths'] for s in p['original_records']]
            recounted=fusion_footprint(records,event)
            assert recounted['complete_paths']==product['complete_paths']
            assert recounted['cell_umi_labels']==product['cell_umi_labels']
        for supplied in fusion_inputs(entry,{s['id']:s for s in data['sources']}):
            fusion,refs,reads=fusion_from_dict(supplied)
            result=reconstruct_fusion(fusion,refs,reads)
            assert result['status']=='unresolved_frame' and result['translations']==[]
            for original in supplied['original_records']:
                records=[pysam.AlignedSegment.fromstring(original[k],header) for k in ('sam','partner_sam')]
                _,c1,p1,c2,p2,s1,s2=event
                row=extract(records[0],c1,p1,c2,p2,fusion.junction_start,records,(s1,s2))
                assert row['sequence']==fusion.sequence
                assert row['blocks']==supplied['fusion']['blocks']


def test_intronic_skip_is_not_a_cigar_deletion_or_coverage_of_deleted_bases():
    for op,label in [('N','splice'),('D','deletion')]:
        r=record(cigar='20M128'+op+'20M')
        result=deletion_footprint([r],dict(start=120,end=248),[dict(exons=[(100,120),(248,268)])])
        assert result['aligned_templates']==dict(left=1,inside=0,right=1)
        assert result['breakpoint_matching_templates'][label]==1
        assert result['breakpoint_matching_templates']['deletion' if op=='N' else 'splice']==0
        assert result['junctions'][0]['operation']==op


def test_dlg5_start_loss_is_not_relabelled_an_internal_inframe_deletion():
    data=load()
    dlg5=next(e for e in data['deletions'] if e['gene']=='DLG5')
    affected=[m for m in dlg5['models'] if any(dlg5['start']<p<=dlg5['end'] for p in m['start_codon_positions'])]
    assert affected and all(m['deleted_coding_bases'] for m in affected)
    for e in data['deletions']:
        if e['gene']!='DLG5':
            assert not any(m['deleted_coding_bases'] for m in e['models'])
        assert not any(p['complete_paths'] for p in e['supplementary_products'])


def test_split_deletion_requires_observed_supplementary_not_mate_or_sa_only():
    first,second=link(record(cigar='30M30S'),record(250,'30H30M',2048))
    assert len(split_deletion_paths([first,second],130,250))==1
    assert split_deletion_paths([first],130,250)==[]
    second.flag=129
    assert split_deletion_paths([first,second],130,250)==[]
