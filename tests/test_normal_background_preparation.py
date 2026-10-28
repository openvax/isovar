"""ORFdata block coordinates are inclusive, including one-base coding blocks."""

from types import SimpleNamespace

from isovar.genomic_intervals import GenomicIntervalIndex
from reports.matched_normal import prepare_sid_background_regions as preparation


def test_ncorf_background_keeps_first_base_and_single_base_blocks(monkeypatch, tmp_path):
    # Native gencode_phase1 c10norep34 / FAM107B block geometry.
    record = SimpleNamespace(contig='10', coding_blocks=((14530515, 14530515), (14553322, 14553383)))
    monkeypatch.setattr(preparation.orfdata, 'load_catalog', lambda path: SimpleNamespace(records=[record]))
    regions = preparation.catalog_regions(tmp_path / 'catalog.json')
    assert regions == [('10', 14530514, 14530515), ('10', 14553321, 14553383)]
    index = GenomicIntervalIndex(regions)
    assert index.covers('10', 14530514, 14530515)
    assert index.covers('10', 14553321, 14553322)
    assert sum(end - start for _, start, end in regions) == 63
