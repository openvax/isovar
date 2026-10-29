"""Antigen references use the same exact/alias/ambiguity contract as variants."""

import pytest
from varcode import ContigResolver

from isovar.antigen_reference import contig_in_reference
from isovar.contig_resolution import UnresolvedContigError


@pytest.mark.parametrize('query,names,expected', [
    ('MT', ['CHRM'], 'CHRM'), ('chrM', ['MT'], 'MT'),
    ('chr1', ['1'], '1'), ('chr1', ['1', 'chr1'], 'chr1'),
    ('HPV16', ['HPV16', 'HPV18'], 'HPV16'), ('chrEBV', ['chrEBV'], 'chrEBV'),
    ('hs38d1', ['hs38d1'], 'hs38d1'), ('0', ['0', 'HPV16'], '0'),
    ('GL000207.1', ['GL000207.1', 'GL000207.2'], 'GL000207.1'),
    ('chrUn_KI270442v1', ['chrUn_KI270442v1'], 'chrUn_KI270442v1'),
])
def test_antigen_reference_exact_and_unique_aliases(query, names, expected):
    assert contig_in_reference(query, names) == expected
    assert contig_in_reference(query, ContigResolver(names)) == expected


@pytest.mark.parametrize('query,names,status', [
    ('chr1', ['1', 'CHR1'], 'ambiguous'),
    ('HPV16', ['HPV18'], 'missing'), ('GL000207.1', ['GL000207.2'], 'missing'),
])
def test_unresolved_names_keep_structured_diagnostics(query, names, status):
    with pytest.raises(UnresolvedContigError) as caught:
        contig_in_reference(query, names)
    assert caught.value.resolution.status == status
    assert caught.value.to_dict()['original_name'] == query
    if status == 'ambiguous':
        assert set(caught.value.resolution.candidates) == set(names)
