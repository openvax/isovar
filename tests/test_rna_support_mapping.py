"""Support recount cost must not rebuild a whole candidate for every read."""

from copy import deepcopy

import pysam
import pytest

from isovar.read_collector import ReadCollector
from isovar.rna_library import prepare_rna_library
import isovar.rna_reconstruction as reconstruction
from isovar.sv_rna import _Positions
from .test_antigen_isoforms import target_models
from .test_sv_rna import aligned, write_bam


@pytest.mark.parametrize('strand', ['+', '-'])
def test_candidate_mapping_is_reused_without_changing_support(tmp_path, monkeypatch, strand):
    targets, paths = target_models(strand)
    reads = [read for index in range(20) for read in
             aligned('original-%d' % index, targets[0].reference.sequence, paths[0])]
    bam = write_bam(tmp_path / 'original.bam', reads)
    collector = ReadCollector(use_secondary_alignments=True, min_mapping_quality=0, use_soft_clipped_bases=False)
    prepared = prepare_rna_library(bam, tmp_path / 'index', read_collector=collector, reference_name='test')
    candidates = [dict(id=str(index), sequence=target.reference.sequence, placements=tuple(mapping), coding=[])
                  for index, (target, mapping) in enumerate(zip(targets, paths))]
    before = deepcopy(candidates)
    candidate_positions = {tuple(path) for path in paths}
    builds = []
    runs = reconstruction._runs

    def counted_runs(positions):
        if tuple(positions) in candidate_positions and (not isinstance(positions, _Positions) or positions.runs is None):
            builds.append(len(positions))
        return runs(positions)

    monkeypatch.setattr(reconstruction, '_runs', counted_runs)
    regions = {'1': [(0, 2000, strand)]}
    with prepared.open() as library, pysam.AlignmentFile(bam) as reader:
        cached = reconstruction._support(library, candidates, regions, collector, reader.header.to_dict(), 'sample', 'original')
        # Original geometries may also be mapped once, but work is independent
        # of the 20 duplicate geometries / 40 matching attempts per candidate.
        assert len(builds) < 6
        builds.clear()
        monkeypatch.setattr(reconstruction, '_Positions', tuple)
        uncached_candidates = deepcopy(before)
        uncached = reconstruction._support(library, uncached_candidates, regions, collector,
                                           reader.header.to_dict(), 'sample', 'original')
    assert len(builds) >= 80
    assert cached == uncached
    assert candidates == uncached_candidates
    assert [c['placements'] for c in candidates] == [c['placements'] for c in before]
