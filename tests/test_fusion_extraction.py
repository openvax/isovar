"""#288: junction bases come from observed full CIGARs, never SA summaries."""


import pysam
import pytest
from isovar import sid_data

from tests.data.fusions.build_osteosarc import extract, query_mapping, blocks
from tests.test_chimeric_phasing import link, record


def test_original_sid_insertion_is_not_relocated_to_sa_end():
    path = sid_data.path('chimeric/osteosarc-ont.sam')
    with pysam.AlignmentFile(path) as bam:
        reads = [r for r in bam if r.query_name == '46de198c-6f57-490f-9f0e-1679fb487660_0']
    partner, = [r for r in reads if r.is_supplementary]
    assert partner.cigarstring == '356H5M1I128M1H'
    mapping = query_mapping(partner)
    length = partner.infer_read_length()
    assert mapping[length-1-360] == ('chr17',63745253,'-')
    assert length-1-361 not in mapping  # Actual inserted query base.
    assert mapping[length-1-362] == ('chr17',63745254,'-')


@pytest.mark.parametrize('hard', [False, True])
def test_observed_partner_required_and_hard_clips_do_not_shift_bases(hard):
    first, second = link(record(cigar='30M30S'), record(200,'30H30M' if hard else '30S30M',2048))
    assert extract(first,'1',130,'1',201,20) is None
    assert extract(first,'1',130,'1',201,20,[first]) is None
    row=extract(first,'1',130,'1',201,20,[first,second])
    assert row['sequence']=='A'*40
    assert row['junction_start']==row['junction_end']==20
    assert row['blocks']==[
        dict(query_start=0,query_end=20,contig='1',reference_start=110,reference_end=130,strand='+'),
        dict(query_start=20,query_end=40,contig='1',reference_start=200,reference_end=220,strand='+')]


@pytest.mark.parametrize('change',['mate','read_group','secondary','missing_sa','sequence'])
def test_unobserved_or_conflicting_partner_cannot_supply_junction(change):
    first,second=link(record(),record(200,'30H30M',2048))
    if change=='mate':
        second.flag |= 129
    elif change=='read_group':
        second.set_tag('RG','different')
    elif change=='secondary':
        second.flag |= 256
    elif change=='missing_sa':
        second.set_tag('SA',None)
    else:
        second.query_sequence='C'*30
    assert extract(first,'1',130,'1',201,20,[first,second]) is None


def test_reverse_strand_blocks_remain_increasing_genomic_intervals():
    mapping={q:('1',100-q,'-') for q in range(4)}
    assert blocks(mapping,0,4)==[dict(query_start=0,query_end=4,contig='1',reference_start=97,reference_end=101,strand='-')]
