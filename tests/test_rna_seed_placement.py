"""A read's alternative placement is not an independent overlap witness."""

from dataclasses import replace

import pytest

from isovar.rna_reconstruction import _assemble
from isovar.sv_rna import RnaObservation, _mirror


def observations(independent=False, reverse=False):
    first = RnaObservation(
        "placement-a+", ("RG", "original", 0), ("sam-a",), "AAAA",
        tuple(("chr1", p, "+") for p in range(4)), (), False, (0, 4), False, False)
    second = replace(first, key="placement-b+", records=("sam-b",),
                     positions=tuple(("chr1", p, "+") for p in range(2, 6)), secondary=not independent,
                     identity=("RG", "independent" if independent else "original", 0))
    if reverse:
        first, second = [replace(observed, key=observed.key[:-1] + "-",
                                 sequence=_mirror(observed.sequence, observed.positions)[0],
                                 positions=_mirror(observed.sequence, observed.positions)[1],
                                 reverse_complement=True) for observed in (first, second)]
    return first, second


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("reverse_order", [False, True])
def test_alternative_placements_do_not_extend_their_own_seed(reverse, reverse_order):
    originals = observations(reverse=reverse)
    chosen = originals[::-1] if reverse_order else originals
    paths, notes = _assemble({o.comparison_key: o for o in chosen}, 10, 2)
    assert set(paths) == {o.comparison_key for o in originals}
    assert "incompatible_segment_placements_retained_as_alternatives" in notes


@pytest.mark.parametrize("reverse", [False, True])
def test_independent_original_read_still_extends(reverse):
    chosen = observations(independent=True, reverse=reverse)
    paths, notes = _assemble({o.comparison_key: o for o in chosen}, 10, 2)
    sequence, positions = paths[0]
    assert len(paths) == 1 and sequence == ("TTTTTT" if reverse else "AAAAAA")
    assert {position[1] for position in positions} == set(range(6))
    assert not notes


def test_shared_path_budget_retains_an_alternative_without_manufacturing_extension():
    chosen = observations()
    paths, notes = _assemble({o.comparison_key: o for o in chosen}, 1, 2)
    assert len(paths) == 1 and len(paths[0][0]) == 4
    assert "path_limit" in notes
