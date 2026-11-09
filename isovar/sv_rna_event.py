"""SV adjacency classification and event-specific discovery seeds."""

from __future__ import annotations

from collections import defaultdict

from .rna_observation import RnaObservationStore, alignment_gaps
from .rna_path_assembly import SEQUENCE_BRIDGE_BASES
from .rna_path_geometry import reverse_position
from .sv_rna_relations import LINKAGE_STATUSES


def _distance(position, breakpoint, side):
    """Signed retained-partner bases between a base and a breakpoint.

    Negative when the base lies past the breakpoint, on the side the adjacency
    removes; ``None`` on another contig or strand.
    """
    contig, base, strand = position
    if (contig, strand) != (breakpoint.contig, breakpoint.strand):
        return None
    return breakpoint.position - base - 1 if (side == "donor") == (strand == "+") else base - breakpoint.position


def is_forward_splice(left, right):
    """Whether the unrearranged reference could make this join by cis-splicing.

    Parameters
    ----------
    left, right : tuple[str, int, str]
        The two flanking placed bases: contig, zero-based base and strand.

    Returns
    -------
    bool
        Whether the join has ordinary same-contig, same-strand forward geometry.
    """
    return left[::2] == right[::2] and (right[1] > left[1] if left[2] == "+" else right[1] < left[1])


class OrientedAdjacency:
    """The nominated oriented adjacency, compared with observed RNA joins.

    Parameters
    ----------
    donor, acceptor : FusionBreakpoint
        Oriented retained-partner boundaries of the nominated DNA adjacency.
    max_shift : int
        Declared bound for junction homology/breakpoint placement shifts.
    annotated : set of tuple
        Pairs of placed bases bounding annotated RNA junctions.
    tolerance : int
        Allowed displacement around annotated joins for seed classification.
    """

    def __init__(self, donor, acceptor, max_shift, annotated, tolerance):
        self.donor, self.acceptor, self.max_shift = donor, acceptor, max_shift
        self.donor_sites = {left for left, _ in annotated}
        self.acceptor_sites = {right for _, right in annotated}
        self.nearby = defaultdict(set)
        for (contig, low, strand), right in annotated:
            for shift in range(-tolerance, tolerance + 1):
                self.nearby[contig, low + shift, strand].add(right)
        self.tolerance = tolerance
        # The adjacency is read-through-ambiguous when any placement on its
        # diagonal joins an annotated exon end to an annotated exon start
        # downstream on the same strand: cis-splicing can make that RNA.
        self.read_through = any(
            is_forward_splice(left, right) and left in self.donor_sites and right in self.acceptor_sites
            for left, right in (self._at(d, -d) for d in range(-max_shift, max_shift + 1)))
        # A duplication-like adjacency joins a base back to one upstream of it
        # on the same strand. Forward joins happen within either copy, so they
        # never cross it.
        left, right = self._at(0, 0)
        self.backward = left[::2] == right[::2] and not is_forward_splice(left, right)

    def _at(self, d, a):
        """The donor and acceptor bases at these signed distances from the breakpoints."""
        donor, acceptor = self.donor, self.acceptor
        left = donor.position - 1 - d if donor.strand == "+" else donor.position + d
        right = acceptor.position + a if acceptor.strand == "+" else acceptor.position - 1 - a
        return (donor.contig, left, donor.strand), (acceptor.contig, right, acceptor.strand)

    def near_annotated(self, left, right):
        """Whether an unannotated join lies within ``tolerance`` bases of an annotated one.

        Long-read aligners commonly misplace annotated junctions by a few
        bases; such joins seed no ordinary-splicing or regional path.
        """
        return any(r[::2] == right[::2] and abs(r[1] - right[1]) <= self.tolerance for r in self.nearby.get(left, ()))

    def relation(self, left, right, unplaced, kind=None):
        """Relation of one join to the adjacency, and its breakpoint assignment.

        Unplaced junction bases (homology or untemplated insertion) may belong
        to either partner, and an aligner may place up to ``max_shift``
        homologous bases past a breakpoint: every such join on the adjacency's
        diagonal is the breakpoint junction. The homology is not verified
        against a genome sequence. Every placement of a read-through-ambiguous
        adjacency is splice-ambiguous. A join just off the diagonal, with one
        side past a breakpoint, is a (noisy) event join, not an unrelated one.
        Off the diagonal, a join must run in the adjacency's direction: a
        forward join (an ordinary splice or a small CIGAR deletion) cannot
        cross a duplication-like adjacency, because one copy can make it.
        """
        d, a = _distance(left, self.donor, "donor"), _distance(right, self.acceptor, "acceptor")
        if d is None or a is None:
            return None, None
        # CIGAR deletions are alignment indels unless they place the nominated
        # adjacency itself. Apply the same rule during seeding and export.
        if kind == "D" and not self.near_breakpoint(d, a, unplaced):
            return None, None
        if self.on_diagonal(d, a, unplaced):
            return ("splice_ambiguous_event_junction" if self.read_through else "breakpoint_junction"), [d, a]
        if self.backward and is_forward_splice(left, right):
            return None, None
        if min(d, a) < 0 and not self.near_breakpoint(d, a, unplaced):
            return None, None
        return ("splice_ambiguous_event_junction" if is_forward_splice(left, right)
                else "event_compatible_junction"), [d, a]

    def on_diagonal(self, d, a, unplaced):
        """Whether a join with this assignment is the adjacency itself.

        Junction bases may be assigned to either partner, but the RNA cannot
        contain reference bases from past both breakpoints (``d + a < 0``).
        """
        return 0 <= d + a <= unplaced and min(d, a) >= -self.max_shift

    def near_breakpoint(self, d, a, unplaced):
        """Within ``max_shift`` of the adjacency: a noisy placement of it, if weaker."""
        return max(d, a) <= unplaced + self.max_shift and min(d, a) >= -self.max_shift

    def clips(self, positions, length):
        """(i, j, assignment) where placed sequence stops exactly at a breakpoint."""
        found = []
        if positions[-1] is None:
            last = next(q for q in range(length - 1, -1, -1) if positions[q] is not None)
            if _distance(positions[last], self.donor, "donor") == 0:
                found.append((last, last + 1, [0, None]))
        if positions[0] is None:
            first = next(q for q, p in enumerate(positions) if p is not None)
            if _distance(positions[first], self.acceptor, "acceptor") == 0:
                found.append((first - 1, first, [None, 0]))
        return found


def canonical_join(left, right):
    """One key for a join and its reverse complement.

    Parameters
    ----------
    left, right : tuple[str, int, str]
        The two flanking placed bases: contig, zero-based base and strand.

    Returns
    -------
    tuple
        One deterministic key shared by a join and its reverse complement.
    """
    return min((left, right), (reverse_position(right), reverse_position(left)))


def classify_join(event, annotated, left, right, unplaced, kind):
    """How one join, seen in this orientation, is used.

    Parameters
    ----------
    event : OrientedAdjacency
        Nominated event and its explicit placement policy.
    annotated : set of tuple
        Known reference junctions.
    left, right : tuple[str, int, str]
        The two flanking placed bases: contig, zero-based base and strand.
    unplaced : int
        Observed unplaced query bases between the flanking bases.
    kind : str or None
        Observed CIGAR/split kind; None retains unknown kind.

    Returns
    -------
    tuple
        (Discovery class, event relation), preserving annotated/ambiguous/unrelated distinctions.

    Notes
    -----
    Returns ``("annotated", relation)``, ``("event", relation)`` for joins at
    or across the adjacency, ``("lazy", relation)`` for ordinary-splice joins
    across the event (seeded later, by support), ``("wobble", None)`` for
    unannotated joins beside annotated ones, or ``(None, relation)``.
    """
    relation, assignment = event.relation(left, right, unplaced, kind=kind)
    if unplaced == 0 and (left, right) in annotated:
        return "annotated", relation
    if relation and (relation != "splice_ambiguous_event_junction" or event.on_diagonal(*assignment, unplaced)):
        return "event", relation
    if kind == "D":
        return None, None  # A small deletion is an indel, not a splice across the event.
    if unplaced == 0 and event.near_annotated(left, right):
        return "wobble", None
    return ("lazy", relation) if relation else (None, None)


def has_breakpoint_clip(read, breakpoints):
    """Conservative clip screen using observed blocks, not terminal D/N spans.

    Parameters
    ----------
    read : pysam.AlignedSegment
        Original alignment record.
    breakpoints : sequence of tuple[str, int]
        Nominated contigs and boundary coordinates.

    Returns
    -------
    bool
        Whether a terminal clipped alignment could seed a breakpoint observation.
    """
    if not any(op == 4 for op, _ in read.cigartuples):
        return False
    blocks = read.get_blocks()
    return bool(blocks) and any(
        read.reference_name == contig and abs(edge - position) <= 1
        for contig, position in breakpoints for edge in (blocks[0][0], blocks[-1][1]))


def has_exon_clip(read, exonic, minimum):
    """Screen a terminal clip beside an exon before constructing read geometry.

    Parameters
    ----------
    read : pysam.AlignedSegment
        Original alignment record.
    exonic : callable
        Whether an oriented placed base lies in a supplied exon.
    minimum : int
        Required terminal clipped seed length.

    Returns
    -------
    bool
        Whether the record could supply an exon-anchored clipped seed.
    """
    if not read.cigartuples or not read.get_blocks():
        return False
    blocks = read.get_blocks()
    # A terminal soft clip can sit inside a terminal hard clip in valid SAM.
    operations = [item for item in read.cigartuples if item[0] != 5]
    edges = ((operations[0], blocks[0][0]),
             (operations[-1], blocks[-1][1] - 1))
    return any(operation == 4 and length >= minimum and any(
        exonic((read.reference_name, base, strand)) for strand in ('+', '-'))
        for (operation, length), base in edges)


def exon_clip_cores(observation, exonic, minimum):
    """Observed exon anchor and terminal unplaced query, in RNA orientation.

    Parameters
    ----------
    observation : RnaObservation
        Original processed observation in RNA orientation.
    exonic : callable
        Supplied exon-membership predicate.
    minimum : int
        Minimum observed unplaced terminal seed length.

    Returns
    -------
    sequence of tuple
        Original query intervals eligible as exon-anchored clip cores.
    """
    placed = [q for q, position in enumerate(observation.positions) if position is not None]
    if not placed:
        return ()
    first, last = placed[0], placed[-1]
    cores = []
    if first >= minimum and exonic(observation.positions[first]):
        cores.append((0, first, 'exon_anchored_clip_partner_unplaced', 'clip'))
    if len(observation.sequence) - last - 1 >= minimum and exonic(observation.positions[last]):
        cores.append((last, len(observation.sequence) - 1, 'exon_anchored_clip_partner_unplaced', 'clip'))
    return cores


def sv_path_seeds(observations, event, annotated, anchored, competing, key=None, deferred=None,
           exonic=None, exon_clipped=()):
    """Seeds at unannotated joins and breakpoint clips, strongest adjacency placement first.

    Parameters
    ----------
    observations : iterable of RnaObservation
        Original observations, retaining all selected placement alternatives.
    event : OrientedAdjacency
        Nominated adjacency and its relation policy.
    annotated : set of tuple
        Known annotated junctions.
    anchored : callable
        Whether a placed base lies in a supplied transcript span.
    competing : mapping
        Original fragments supporting competing annotated junctions.
    key : object, optional
        Deferred junction category to seed; None selects event seeds.
    deferred : mapping, optional
        Previously classified lazy junction candidates.
    exonic : callable, optional
        Exon-membership predicate for optional clipped seeds.
    exon_clipped : iterable, optional
        Original segment identities selected by the exon-clip screen.

    Returns
    -------
    list
        Ranked (sequence/placement seed, original evidence description) pairs.

    Notes
    -----
    Without ``key``, only joins at or across the event (and breakpoint clips)
    seed; other joins are added to ``deferred[relation][canonical join]`` so
    that all of their reads are judged together. With ``key``, only that
    (canonical) join seeds, in each orientation that is sense to a gene or
    crosses the event. Annotated joins across the event go to ``competing``.
    """
    seeds = {}
    for o in observations:
        cores = []
        for i, j, kind in o.breaks:
            left, right = o.positions[i], o.positions[j]
            use, relation = classify_join(event, annotated, left, right, j - i - 1, kind)
            if use == "annotated":
                if relation:
                    competing[left, right].add(o.fragment)
                continue
            if use == "event" and key is None:
                cores.append((i, j, relation, kind))
                continue
            regional = (use is None and relation is None and (anchored(left) or anchored(right))
                        and not event.relation(reverse_position(right), reverse_position(left), j - i - 1)[0] and kind != "D")
            if use != "lazy" and not regional:
                continue
            relation = relation or "regional_novel_junction"
            if key is None:
                if deferred is not None:
                    deferred[relation][canonical_join(left, right)].add(o.identity)
            elif canonical_join(left, right) == key:
                cores.append((i, j, relation, kind))
        if key is None:
            cores += [(i, j, "breakpoint_clip_partner_unplaced", "clip")
                      for i, j, _ in event.clips(o.positions, len(o.sequence))]
            if exonic is not None and o.identity in exon_clipped:
                cores += exon_clip_cores(o, exonic, SEQUENCE_BRIDGE_BASES)
        for i, j, relation, kind in cores:
            seed = seeds.setdefault((o.sequence[i:j + 1], o.positions[i:j + 1]), dict(
                relation=relation, fragments=set(), observations=[], kinds=set(),
                near_breakpoint=relation in ("breakpoint_junction", "event_compatible_junction",
                                             "splice_ambiguous_event_junction")
                and event.near_breakpoint(*event.relation(o.positions[i], o.positions[j], j - i - 1)[1], j - i - 1)))
            seed["fragments"].add(o.fragment)
            seed["observations"].append((o, i))
            seed["kinds"].add(kind)
    rank = {status: i for i, status in enumerate(LINKAGE_STATUSES)}
    # Placements of the adjacency are ordered by support alone, so the
    # strongest is the reference its weaker variants are pruned against.
    return sorted(seeds.items(), key=lambda item: (
        not item[1]["near_breakpoint"], 0 if item[1]["near_breakpoint"] else
        rank.get(item[1]["relation"], len(rank)),
        -len(item[1]["fragments"]), item[0][0], repr(item[0][1])))


def contains_seed(path, core):
    """Check whether an observed path contains an event seed.

    Parameters
    ----------
    path, core : tuple[str, sequence]
        Path and core sequence/placements; unplaced core bases do not assert
        a genomic coordinate.

    Returns
    -------
    bool
        Exact seed bases with agreement at every placed core coordinate.
    """
    sequence, positions = path
    start = sequence.find(core[0])
    while start >= 0:
        if all(x is None or x == y for x, y in zip(core[1], positions[start:start + len(core[0])])):
            return True
        start = sequence.find(core[0], start + 1)
    return False


class SvObservationStore(RnaObservationStore):
    """Lazy original observations with SV-specific join and seed classification.

    Parameters
    ----------
    *args, **kwargs : arguments
        RnaObservationStore construction arguments. The subclass adds
        SV join classification, gap kinds and exon-clipped seed identities.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.gaps = defaultdict(dict)
        self.gap_kinds = defaultdict(set)
        self.exon_clipped = set()

    def classify(self, event, annotated, anchored, exonic=None):
        """Sort segments by their CIGAR joins, without building single records.

        Returns segments to build now (multi-record paths, joins at or across
        the event, possible breakpoint clips); ``lazy[relation][canonical
        join]`` segments for splice-ambiguous and regional unannotated joins;
        and fragments of annotated joins crossing the event.
        """
        now, lazy, competing = set(), defaultdict(lambda: defaultdict(set)), defaultdict(set)
        breakpoints = [(b.contig, b.position) for b in (event.donor, event.acceptor)]
        for identity, records in self.groups.items():
            if len(records) > 1:
                now.add(identity)
                continue
            read, = records
            if read.has_tag("SA"):
                self.reasons["SA_declared_piece_not_linked"] += 1  # Its declared pieces were not collected.
            for left, right, kind, inserted in alignment_gaps(read, self.collector):
                self.gaps[left, right][identity] = inserted
                self.gap_kinds[left, right].add(kind)
                oriented = [(left, right), (reverse_position(right), reverse_position(left))]
                classes = [classify_join(event, annotated, *join, len(inserted), kind) for join in oriented]
                uses = [use for use, _ in classes]
                if "annotated" in uses:
                    for (use, relation), join in zip(classes, oriented):
                        if use == "annotated" and relation:
                            competing[join].add(identity[:2])
                elif "event" in uses:
                    now.add(identity)
                elif "wobble" in uses:
                    self.reasons["join_near_annotated_junction"] += 1
                elif "lazy" in uses:
                    lazy[next(r for use, r in classes if use == "lazy")][canonical_join(left, right)].add(identity)
                elif kind != "D" and all(r is None for _, r in classes) and any(
                        anchored(p) for p in (left, right, reverse_position(left), reverse_position(right))):
                    lazy["regional_novel_junction"][canonical_join(left, right)].add(identity)
            # SAM permits soft clips inside terminal hard clips.
            if self.collector.use_soft_clipped_bases and has_breakpoint_clip(read, breakpoints):
                now.add(identity)
            if (exonic is not None and self.collector.use_soft_clipped_bases
                    and has_exon_clip(read, exonic, SEQUENCE_BRIDGE_BASES)):
                now.add(identity)
                self.exon_clipped.add(identity)
        return now, lazy, competing

    def join_segments(self, left, right):
        """Segments whose own path makes this join: built ones plus never-built single records.

        A built segment without the join (a failed build, or a trimmed read
        end) does not count.
        """
        unbuilt = set(self.gaps.get((left, right), ())) | set(self.gaps.get((reverse_position(right), reverse_position(left)), ()))
        return {o.identity for o, *_ in self.joins.get((left, right), ())} | {i for i in unbuilt if i not in self.built}
