"""How joins relate to a nominated adjacency (#665)."""

import pytest

from isovar.fusion import FusionBreakpoint
from isovar.sv_rna_audit import _summary
from isovar.sv_rna_event import OrientedAdjacency, classify_join


def event(donor, acceptor, max_shift=5):
    return OrientedAdjacency(FusionBreakpoint("chr1", donor, "+"), FusionBreakpoint("chr1", acceptor, "+"),
                  max_shift, annotated=set(), tolerance=3)


def join(left, right):
    return ("chr1", left, "+"), ("chr1", right, "+")


# A tandem-duplication-like adjacency joins base 999 back to base 500.
DUPLICATION = event(1000, 500)
# A deletion-like adjacency joins base 499 forward to base 1000.
DELETION = event(500, 1000)


def test_duplication_breakpoint_and_backward_joins_are_event_linked():
    assert DUPLICATION.relation(*join(999, 500), 0)[0] == "breakpoint_junction"
    assert DUPLICATION.relation(*join(990, 520), 0)[0] == "event_compatible_junction"


def test_forward_joins_never_cross_a_duplication_like_adjacency():
    # One copy of the unrearranged reference makes every forward join: an
    # ordinary splice inside the duplicated span, or one spanning all of it.
    assert DUPLICATION.relation(*join(600, 700), 0) == (None, None)
    assert DUPLICATION.relation(*join(400, 1100), 0) == (None, None)


def test_forward_joins_across_a_deletion_like_adjacency_stay_splice_ambiguous():
    assert DELETION.relation(*join(400, 1100), 0)[0] == "splice_ambiguous_event_junction"


def test_small_deletions_relate_only_as_placements_of_the_adjacency():
    # The adjacency itself, as a CIGAR deletion: a short deletion SV.
    assert classify_join(DELETION, set(), *join(499, 1000), 0, "D") == ("event", "breakpoint_junction")
    # The same span, crossed far from the breakpoints: an alignment indel.
    assert classify_join(DELETION, set(), *join(400, 1100), 0, "D") == (None, None)
    assert DELETION.relation(*join(400, 1100), 0, kind="D") == (None, None)
    assert DELETION.relation(*join(497, 1002), 0, kind="D")[0] == "splice_ambiguous_event_junction"


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("kind", ["N", "D", "split"])
def test_duplication_direction_in_transcript_orientation(strand, kind):
    donor, acceptor = (1000, 500) if strand == "+" else (500, 1000)
    nominated = OrientedAdjacency(FusionBreakpoint("chr1", donor, strand), FusionBreakpoint("chr1", acceptor, strand),
                       5, annotated=set(), tolerance=3)
    left, right = (600, 700) if strand == "+" else (700, 600)
    assert nominated.relation(("chr1", left, strand), ("chr1", right, strand), 0, kind=kind) == (None, None)
    left, right = (990, 520) if strand == "+" else (520, 990)
    if kind != "D":
        assert nominated.relation(("chr1", left, strand), ("chr1", right, strand), 0, kind=kind)[0] == "event_compatible_junction"


def path(relations, linked):
    return dict(junctions=[dict(relation=r, direct_support=dict(fragments=3)) for r in relations],
                proteins=[dict(event_linked=linked)])


def test_audit_counts_breakpoint_and_event_compatible_junctions_separately():
    paths = [path(["breakpoint_junction"], True), path(["event_compatible_junction"], True),
             path(["splice_ambiguous_event_junction"], False)]
    row = _summary(dict(event_id="e", source_id="s", view_id="v", status="ok",
                        assessment=dict(paths=paths)), "artifact")
    assert row["paths_with_direct_junction"] == 1
    assert row["paths_with_event_compatible_junction"] == 1
    assert row["paths_with_event_linked_protein"] == 2
