"""Panel budgets and exact evidence remain independent of unrelated loci."""

from dataclasses import replace

import pytest

from isovar import IsoformTarget
from tests.test_antigen_isoforms import reconstruct, target_models
from tests.test_sv_rna import aligned, record


def shifted(target, *, gene_id, distance=1000):
    reference = target.reference
    return IsoformTarget(gene_id, gene_id, replace(
        reference, transcript_id=gene_id + '-transcript',
        exons=tuple((a + distance, b + distance) for a, b in reference.exons)))


def test_independent_gene_budgets_recover_later_sequence_without_capping_support(tmp_path):
    targets, positions = target_models()
    first = targets[0]
    second = shifted(first, gene_id='LATER')
    sequence = second.reference.sequence
    mutant = sequence[:75] + ('A' if sequence[75] != 'A' else 'T') + sequence[76:]
    other_positions = [(c, p + 1000, s) for c, p, s in positions[0]]
    reads = [r for i in range(3) for r in aligned('a%d' % i, sequence, positions[0])]
    reads += aligned('z', mutant, other_positions)
    panel = reconstruct(tmp_path, reads, [first, second], max_discovery_reads=1)
    genes = reconstruct(tmp_path, reads, [first, second], max_discovery_reads=1, discovery_partition='gene')
    assert not any(row['sequence'] == mutant for row in panel['isoforms'])
    assert any(row['sequence'] == mutant and row['supporting_fragments'] == 1 for row in genes['isoforms'])
    assert panel['support']['target_fragments'] == genes['support']['target_fragments'] == 4
    searches = {row['gene_id']: row for row in genes['discovery']['genes']}
    assert searches['G']['selected_segments'] == searches['LATER']['selected_segments'] == 1
    assert searches['G']['truncated'] and not searches['LATER']['truncated']
    assert 'selected_segments' not in genes['discovery']  # No sum masquerading as unique segments.


def test_shared_paralog_read_is_counted_once_across_gene_partitions(tmp_path):
    targets, positions = target_models()
    first = targets[0]
    second = shifted(first, gene_id='PARALOG')
    other = [(c, p + 1000, s) for c, p, s in positions[0]]
    reads = aligned('shared', first.reference.sequence, positions[0])
    reads += aligned('shared', first.reference.sequence, other, flag=256)
    result = reconstruct(tmp_path, reads, [first, second], discovery_partition='gene')
    assert result['support']['target_fragments'] == 1
    assert sum(row['fragments'] for row in result['support']['equivalence_classes']) == 1
    assert len(result['support']['witnesses'][0]['compatible_isoforms']) == 2


@pytest.mark.parametrize('strand', ['+', '-'])
def test_indexed_candidates_match_exhaustive_results_with_fewer_exact_comparisons(tmp_path, monkeypatch, strand):
    import isovar.antigen_isoforms as antigen
    import isovar.rna_reconstruction as shared

    targets, positions = target_models(strand)
    first = targets[0]
    targets += [shifted(first, gene_id='G%d' % i, distance=i * 1000) for i in range(1, 12)]
    # Splice variation plus shared/paralog placement exercises competing frames
    # and the negative-strand coordinate path, not merely an empty index.
    sequence = first.reference.sequence[:210] + first.reference.sequence[240:]
    mapping = positions[0][:210] + positions[0][240:]
    reads = aligned('observed', sequence, mapping)
    reads += aligned('observed', sequence, [(c, p + 1000, s) for c, p, s in mapping], flag=256)
    calls = []
    original = shared._matches
    def counted(*args, **kwargs):
        calls.append(None)
        return original(*args, **kwargs)
    monkeypatch.setattr(shared, '_matches', counted)
    indexed = reconstruct(tmp_path, reads, targets)
    indexed_calls = len(calls)

    class Exhaustive:
        def __init__(self, placements):
            self.keys = tuple(placements)

        def overlapping(self, positions, lo, hi):
            return {(key, 0) for key in self.keys}

    monkeypatch.setattr(shared, '_PlacementIndex', Exhaustive)
    monkeypatch.setattr(antigen, '_PlacementIndex', Exhaustive)
    calls.clear()
    exhaustive = reconstruct(tmp_path, reads, targets)
    assert indexed == exhaustive
    assert indexed_calls < len(calls) / 2


def test_invalid_partition_is_rejected(tmp_path):
    with pytest.raises(ValueError, match='discovery_partition'):
        reconstruct(tmp_path, [], discovery_partition='unknown')


def test_inserted_bases_retain_the_observed_coding_path(tmp_path):
    targets, _ = target_models()
    reference = targets[0].reference
    sequence = reference.sequence[:75] + 'GGG' + reference.sequence[75:]
    read = record('insertion', '1', 100, '75M3I105M120N90M110N180M', sequence)
    result = reconstruct(tmp_path, [read], targets, discovery_partition='gene')
    observed, = [row for row in result['isoforms'] if row['origin'] == 'observed_assembly']
    assert observed['sequence'] == sequence
    assert observed['whole_path_segments'] == [['a', 'insertion', 0]]
    assert observed['coding'][0]['status'] == 'complete_CDS'
    assert len(observed['coding'][0]['amino_acids']) == (reference.cds_end - reference.cds_start) // 3
