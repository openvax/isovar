"""Acquisition coverage remains exact at boundaries and does not infer aliases."""

import pytest
from isovar.genomic_intervals import GenomicIntervalIndex


def test_merged_coverage_and_half_open_boundaries():
    index = GenomicIntervalIndex([('chr1', 20, 30), ('chr1', 10, 20), ('chr1', 15, 25), ('chr1', 40, 50)])
    assert index.intervals() == (('chr1', 10, 30), ('chr1', 40, 50))
    assert index.covers('chr1', 10, 30)
    assert not index.covers('chr1', 29, 41)
    assert index.overlaps('chr1', 29, 41)
    assert not index.overlaps('chr1', 30, 40)
    assert not index.covers('1', 10, 20)
    assert not index.overlaps('unknown', 10, 20)


@pytest.mark.parametrize('interval', [('chr1', -1, 3), ('chr1', 3, 3), ('', 1, 3), ('chr1', 1.5, 3)])
def test_invalid_regions_fail_at_construction(interval):
    with pytest.raises(ValueError, match='Genomic intervals'):
        GenomicIntervalIndex([interval])
