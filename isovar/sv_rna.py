"""Reconstruct RNA paths around one nominated structural variant.

Starting from an oriented DNA adjacency and annotated transcript models, this
module retrieves regional alignments, seeds candidate paths at unannotated RNA
junctions, extends them through exact, coordinate-anchored overlaps and
transfers reading frames from exact collinear CDS anchors.

Results are exploratory candidates, not the validated supplied-fusion schema.
RNA sequence support, reading frame and linkage to the nominated event are
independent evidence axes. A reference base never fills an unobserved RNA
base, an SA tag never substitutes for its record, and a branch that reads
cannot resolve is kept rather than decided by choosing one read.

Coordinates are 0-based. Each path base has one ``(contig, position, strand)``
placement in the path's 5'->3' orientation, or ``None`` when no single observed
placement exists (insertions, soft clips, junction homology).
"""

from bisect import bisect_right
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from hashlib import sha256
from itertools import combinations
from types import SimpleNamespace

from .chimeric_alignment import compatible_phasing_alignments, source_alignment_paths_from_pysam
from .default_parameters import (
    FUSION_PEPTIDE_LENGTHS, SV_ASSEMBLE, SV_BREAKPOINT_WINDOW, SV_MAX_PATHS, SV_MAX_QUERIES, SV_MAX_RECORDS,
    SV_MIN_ALTERNATIVE_FRACTION, SV_MIN_ALTERNATIVE_FRAGMENTS, SV_MIN_ANCHOR_BASES, SV_MIN_OVERLAP,
)
from .fusion import FusionBlock, FusionBreakpoint, FusionReference
from .genetic_code import standard_genetic_code
from .read_collector import ReadCollector
from .read_end_inference import reverse_complement
from .read_identity import source_alignments_from_pysam

_BIN = 64  # Genomic bin width of the observation index.
_HOP_WINDOW = 1000  # Nearby SA/mate targets share one indexed query.
_MAX_SEGMENT_PATHS = 64  # Alternative supplementary groupings per segment.

# Linkage of a path to the nominated adjacency, strongest first. The first
# three are event-linked; the last two may arise without the rearrangement.
LINKAGE_STATUSES = (
    "breakpoint_junction",  # The RNA join reproduces the DNA adjacency.
    "event_compatible_junction",  # Donor-side to acceptor-side join, e.g. spliced.
    "breakpoint_clip_partner_unplaced",  # Unaligned sequence at a breakpoint.
    "splice_ambiguous_event_junction",  # Also an ordinary forward splice.
    "regional_novel_junction",  # Unannotated join not crossing the event.
)
_SUPPORT_GATED = set(LINKAGE_STATUSES[3:])


def segment_identity(read):
    """(read group, QNAME, mate bits): one sequenced segment in one library."""
    return source_alignments_from_pysam(read, read.query_name)[0][0]


def _record_id(read):
    return sha256(read.to_string().encode()).hexdigest()


def _merge(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def _inside(intervals, position):
    i = bisect_right(intervals, (position, float("inf"))) - 1
    return i >= 0 and intervals[i][0] <= position < intervals[i][1]


def _hop_targets(read):
    """SA and mate locations: retrieval hints only, never observations."""
    targets, notes = [], set()
    if read.has_tag("SA"):
        try:
            for entry in read.get_tag("SA").rstrip(";").split(";"):
                contig, position = entry.split(",")[:2]
                targets.append((contig, int(position) - 1))
        except (AttributeError, ValueError):
            notes.add("unparseable_SA_tag")
    if read.is_paired:
        if read.next_reference_id >= 0 and read.next_reference_start >= 0:
            targets.append((read.next_reference_name, read.next_reference_start))
        else:
            notes.add("mate_location_unavailable")
    return targets, notes


def collect_sv_records(bam, regions, references, max_records=SV_MAX_RECORDS, max_queries=SV_MAX_QUERIES,
                       breakpoints=(), breakpoint_window=SV_BREAKPOINT_WINDOW):
    """Fetch regional records, then SA/mate records of the same fragments.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Indexed alignments.
    regions : iterable of (str, int, int)
        Explicit 0-based half-open search intervals. Exons of every reference
        model are searched too.
    references : iterable of FusionReference
    max_records, max_queries : int
        Resource limits. Reaching one is reported, never silent.
    breakpoints : iterable of FusionBreakpoint
        Each is searched ``breakpoint_window`` bases to either side.

    Returns
    -------
    records : list of pysam.AlignedSegment
        Unmodified records, including competing placements, in a stable order.
    acquisition : dict
        Searched intervals, hop queries and limitations. Records reachable only
        through unplaced mates or genome-wide realignment are not assessed.
    """
    lengths = dict(zip(bam.references, bam.lengths))
    spans = defaultdict(list)
    for breakpoint in breakpoints:
        if breakpoint.contig not in lengths:
            raise ValueError("Breakpoint contig %r is absent from the alignments" % breakpoint.contig)
        spans[breakpoint.contig].append((max(0, breakpoint.position - breakpoint_window),
                                         min(lengths[breakpoint.contig], breakpoint.position + breakpoint_window)))
    for contig, start, end in list(regions) + [(r.contig, a, b) for r in references for a, b in r.exons]:
        if (contig not in lengths or type(start) is not int or type(end) is not int
                or not 0 <= start < end <= lengths[contig]):
            raise ValueError("Invalid or unavailable SV search interval: %r" % ((contig, start, end),))
        spans[contig].append((start, end))
    searched = {contig: _merge(intervals) for contig, intervals in spans.items()}
    pending = [(contig, start, end, None) for contig in sorted(searched) for start, end in searched[contig]]
    records, hops, requested, limitations = {}, [], set(), set()
    queries = 0
    while pending:
        targets = defaultdict(set)
        for i, (contig, start, end, fragments) in enumerate(pending):
            if queries >= max_queries or "record_limit" in limitations:
                limitations.add("query_limit" if queries >= max_queries else "record_limit")
                hops.extend(dict(contig=c, start=s, end=e, fragments=len(f), fetched=False)
                            for c, s, e, f in pending[i:] if f is not None)
                break
            queries += 1
            if fragments is not None:
                hops.append(dict(contig=contig, start=start, end=end, fragments=len(fragments), fetched=True))
            for read in bam.fetch(contig, start, end):
                if read.query_name is None:
                    continue
                fragment = segment_identity(read)[:2]
                if fragments is not None and fragment not in fragments:
                    continue
                sam = read.to_string()
                if sam in records:
                    continue
                if len(records) >= max_records:
                    limitations.add("record_limit")
                    break
                records[sam] = read
                hop_targets, notes = _hop_targets(read)
                limitations.update(notes)
                for target_contig, position in hop_targets:
                    if target_contig not in lengths or not 0 <= position < lengths[target_contig]:
                        limitations.add("unavailable_SA_or_mate_contig")
                    elif (not _inside(searched.get(target_contig, ()), position)
                            and (target_contig, position, fragment) not in requested):
                        requested.add((target_contig, position, fragment))
                        targets[target_contig, position].add(fragment)
        if limitations & {"query_limit", "record_limit"}:
            break
        pending = []
        for contig, position in sorted(targets):
            if pending and pending[-1][0] == contig and position - pending[-1][2] < _HOP_WINDOW:
                previous = pending[-1]
                pending[-1] = (contig, previous[1], position + 1, previous[3] | targets[contig, position])
            else:
                pending.append((contig, position, position + 1, frozenset(targets[contig, position])))
    return sorted(records.values(), key=_record_id), dict(
        searched_regions=[[c, a, b] for c in sorted(searched) for a, b in searched[c]],
        hop_queries=hops, queries=queries, records=len(records), limitations=sorted(limitations),
        scope="indexed_regions_and_observed_SA_or_mate_locations",
        genome_wide_competing_placements_assessed=False)


@dataclass(frozen=True, eq=False)
class RnaObservation:
    """One sequenced segment's observed path, in one orientation.

    ``breaks`` lists consecutive placed bases ``(i, j, kind)`` which are not
    genomic neighbours; kind is the CIGAR ``N``/``D`` or ``split`` between
    supplementary records. ``reverse`` marks the reverse complement of the
    original sequencing orientation.
    """

    key: str
    identity: tuple
    records: tuple
    sequence: str
    positions: tuple
    breaks: tuple
    reverse: bool
    query_interval: tuple
    missing_qualities: bool
    secondary: bool

    @property
    def fragment(self):
        return self.identity[:2]


def _flip(position):
    contig, base, strand = position
    return contig, base, "-" if strand == "+" else "+"


def _mirror(sequence, positions):
    return reverse_complement(sequence), tuple(None if p is None else _flip(p) for p in reversed(positions))


def _successor(position):
    contig, base, strand = position
    return contig, base + 1 if strand == "+" else base - 1, strand


def _runs(positions):
    """Collinear placed runs: (query start, query end, contig, low, high, strand)."""
    runs = []
    for q, position in enumerate(positions):
        if position is None:
            continue
        if runs and runs[-1][1] == q and position == _successor(positions[q - 1]):
            start, _, contig, low, high, strand = runs[-1]
            runs[-1] = (start, q + 1, contig, min(low, position[1]), max(high, position[1] + 1), strand)
        else:
            runs.append((q, q + 1, position[0], position[1], position[1] + 1, position[2]))
    return runs


def _breaks(positions):
    placed = [q for q, p in enumerate(positions) if p is not None]
    return [(i, j) for i, j in zip(placed, placed[1:]) if positions[j] != _successor(positions[i])]


def _ineligible(read, collector):
    """The ReadCollector filters, applied to each record separately."""
    if read.is_unmapped:
        return "unmapped"
    if read.is_qcfail:
        return "qc_fail"
    if not read.query_sequence:
        return "sequence_unavailable"
    if read.is_duplicate and not collector.use_duplicate_reads:
        return "duplicate"
    if read.is_secondary and not collector.use_secondary_alignments:
        return "secondary"
    if read.mapping_quality < collector.min_mapping_quality:
        return "mapping_quality"
    if read.query_qualities is None and not collector.use_reads_without_base_qualities:
        return "qualities_unavailable"
    return None


def _linked_pieces(first, second, group):
    """Reuse the reciprocal SA-path validator on two actual records of one segment."""
    sources = [source_alignments_from_pysam(r, r.query_name) for r in (first, second)]
    if sources[0] == sources[1]:
        return False
    placements = defaultdict(set)
    for read in group:
        for segment, placement in source_alignments_from_pysam(read, read.query_name):
            placements[segment].add(placement)
    observations = []
    for read, source in zip((first, second), sources):
        pairs = read.get_aligned_pairs(matches_only=True)
        if not pairs:
            return False
        q, _ = pairs[len(pairs) // 2]
        observations.append(SimpleNamespace(
            source_alignments=source,
            source_alignment_paths=source_alignment_paths_from_pysam(read, source, (q, q + 1))))
    return compatible_phasing_alignments(dict(sources[0]), *observations, placements)


def _segment_groupings(group, reasons):
    """Maximal sets of mutually linked records; alternatives are never unioned."""
    links = {(i, j) for i, j in combinations(range(len(group)), 2) if _linked_pieces(group[i], group[j], group)}
    subsets = {frozenset([i]) for i in range(len(group))}
    for i in range(len(group)):
        for subset in sorted(subsets, key=sorted):
            if i not in subset and all(tuple(sorted((i, j))) in links for j in subset):
                subsets.add(subset | {i})
        if len(subsets) > _MAX_SEGMENT_PATHS:
            reasons["segment_grouping_limit"] += 1
            break
    maximal = [s for s in subsets if not any(s < other for other in subsets)]
    return [tuple(group[i] for i in sorted(s)) for s in sorted(maximal, key=sorted)]


def _record_placements(read, view):
    """Actual CIGAR placements and gap kinds, in original sequencing coordinates."""
    length = read.infer_read_length()
    hard = read.cigartuples[0][1] if read.cigartuples[0][0] == 5 else 0
    start, end = view.sequenced_interval(0, len(view.sequence))
    strand = "-" if read.is_reverse else "+"

    def original(q):
        return length - 1 - hard - q if read.is_reverse else hard + q

    placements, kinds = {}, {}
    q, r, previous, gap = 0, read.reference_start, None, None
    for operation, n in read.cigartuples:
        if operation in (0, 7, 8):
            for k in range(n):
                if start <= original(q + k) < end:
                    placements[original(q + k)] = (read.reference_name, r + k, strand)
            if previous is not None and gap is not None:
                kinds[tuple(sorted((original(previous), original(q))))] = gap
            previous, gap = q + n - 1, None
            q, r = q + n, r + n
        elif operation in (1, 4):
            q += n
        elif operation in (2, 3):
            gap = "N" if operation == 3 or gap == "N" else "D"
            r += n
    return placements, kinds


def _segment_path(records, collector):
    """Join one segment's linked records by original query offset."""
    bases, placements, sources, kinds = {}, defaultdict(set), defaultdict(set), {}
    for index, read in enumerate(records):
        view = collector.read_sequence_view(read)
        start, _ = view.sequenced_interval(0, len(view.sequence))
        sequence = (reverse_complement(view.sequence) if read.is_reverse else view.sequence).upper()
        for q, base in enumerate(sequence, start):
            if bases.setdefault(q, base) != base:
                return "conflicting_segment_sequence"
        mapping, gaps = _record_placements(read, view)
        for q, position in mapping.items():
            placements[q].add(position)
            sources[q].add(index)
        kinds.update(gaps)
    # Bases two pieces place differently (junction homology) have no single
    # observed placement; they stay in the sequence, unplaced.
    unique = {q: next(iter(p)) for q, p in placements.items() if len(p) == 1}
    if not unique:
        return "no_unique_placement"
    lo, hi = (min(bases), max(bases) + 1) if collector.use_soft_clipped_bases else (min(unique), max(unique) + 1)
    if any(q not in bases for q in range(lo, hi)):
        return "unavailable_query_bases"
    sequence = "".join(bases[q] for q in range(lo, hi))
    if set(sequence) - set("ACGT"):
        return "ambiguous_bases"
    positions = tuple(unique.get(q) for q in range(lo, hi))
    breaks = tuple(
        (i, j, kinds.get((i + lo, j + lo), "split" if sources[i + lo].isdisjoint(sources[j + lo]) else "gap"))
        for i, j in _breaks(positions))
    return sequence, positions, breaks, (lo, hi)


def segment_observations(records, collector):
    """Build both orientations of every observed segment path.

    Returns
    -------
    observations : list of RnaObservation
    excluded : dict
        Record ID to the ReadCollector filter which excluded it.
    reasons : collections.Counter
        Segment paths which could not be built, and single-record paths whose
        declared SA pieces were not observed and linked, by reason.
    """
    groups, excluded, reasons = defaultdict(list), {}, Counter()
    for read in records:
        reason = _ineligible(read, collector)
        if reason:
            excluded[_record_id(read)] = reason
        else:
            groups[segment_identity(read)].append(read)
    observations = []
    for identity in sorted(groups):
        group = sorted(groups[identity], key=_record_id)
        for members in _segment_groupings(group, reasons):
            if len(members) == 1 and members[0].has_tag("SA"):
                reasons["SA_declared_piece_not_linked"] += 1
            built = _segment_path(members, collector)
            if isinstance(built, str):
                reasons[built] += 1
                continue
            sequence, positions, breaks, interval = built
            record_ids = tuple(sorted(_record_id(r) for r in members))
            key = sha256("".join(record_ids).encode()).hexdigest()[:24]
            common = dict(identity=identity, records=record_ids, query_interval=interval,
                          missing_qualities=any(r.query_qualities is None for r in members),
                          secondary=any(r.is_secondary for r in members))
            n = len(sequence)
            observations.append(RnaObservation(key + "+", sequence=sequence, positions=positions,
                                               breaks=breaks, reverse=False, **common))
            sequence, positions = _mirror(sequence, positions)
            observations.append(RnaObservation(
                key + "-", sequence=sequence, positions=positions, reverse=True,
                breaks=tuple(sorted((n - 1 - j, n - 1 - i, kind) for i, j, kind in breaks)), **common))
    return observations, excluded, reasons


class _ObservationIndex:
    """Observations binned by their placed collinear runs, strand-specifically."""

    def __init__(self, observations):
        self.observations = {o.key: o for o in observations}
        self.bins = defaultdict(list)
        for o in observations:
            for start, _, contig, low, high, strand in _runs(o.positions):
                for b in range(low // _BIN, (high - 1) // _BIN + 1):
                    self.bins[contig, strand, b].append((o.key, start, low, high))

    def overlapping(self, positions, lo, hi):
        """(key, path offset of the observation start) sharing a placed base in positions[lo:hi]."""
        found = set()
        for start, _, contig, low, high, strand in _runs(positions[lo:hi]):
            for b in range(low // _BIN, (high - 1) // _BIN + 1):
                for key, other, other_low, other_high in self.bins.get((contig, strand, b), ()):
                    if other_low >= high or low >= other_high:
                        continue
                    if strand == "+":
                        shared = low if low > other_low else other_low
                        path_q, observed_q = start + shared - low, other + shared - other_low
                    else:
                        shared = (high if high < other_high else other_high) - 1
                        path_q, observed_q = start + high - 1 - shared, other + other_high - 1 - shared
                    found.add((key, lo + path_q - observed_q))
        return found


def _agrees(sequence, positions, observation, offset, lo, hi):
    """Identical bases and no conflicting placement over path[lo:hi]."""
    a, b = lo - offset, hi - offset
    if observation.sequence[a:b] != sequence[lo:hi]:
        return False
    mine, theirs = positions[lo:hi], observation.positions[a:b]
    return mine == theirs or all(x is None or y is None or x == y for x, y in zip(mine, theirs))


def _supported(support, min_fragments, min_fraction):
    """Alternatives which meet both thresholds; empty when none does."""
    best = max(support.values())
    return {key for key, count in support.items() if count >= min_fragments and count >= min_fraction * best}


def _votes(active, t):
    """Group reads by their (base, placement) at path offset ``t``.

    An unplaced base agrees with every placement of the same base, so a
    soft-clipped or inserted copy of a base does not split its placed support.
    """
    votes = defaultdict(list)
    for observation, offset in active:
        if t - offset < len(observation.sequence):
            votes[observation.sequence[t - offset], observation.positions[t - offset]].append((observation, offset))
    for base, position in [key for key in votes if key[1] is None]:
        placed = [key for key in votes if key[0] == base and key[1] is not None]
        if placed:
            for key in placed:
                votes[key].extend(votes[base, None])
            del votes[base, None]
    return votes


def _extend(sequence, positions, index, min_overlap, thresholds, budget, pruned, notes, mirrored):
    """Extend 3' through exact overlaps which share a placed base with the path.

    Continuations are added while every agreeing read agrees on the next base
    and its placement. At a disagreement, weakly supported alternatives are
    pruned (and reported) only when another alternative meets both
    thresholds; otherwise every alternative becomes its own path. Support is
    counted over every read overlapping the path end. A branch keeps the reads
    which chose it, so it can continue through unplaced sequence (e.g.
    soft-clipped tails) where no new read can be anchored.
    """
    done, stack = [], [(sequence, positions, {}, ())]
    while stack:
        sequence, positions, used, carried = stack.pop()
        n = len(sequence)
        need = min(min_overlap, n)
        extenders = []
        # Carried reads already agree with every base they voted for.
        known = {(o.key, offset) for o, offset in carried}
        for key, offset in sorted(index.overlapping(positions, n - need, n) | known):
            observation = index.observations[key]
            if (offset <= n - need and offset + len(observation.sequence) > n
                    and used.get(observation.identity, key) == key
                    and ((key, offset) in known
                         or _agrees(sequence, positions, observation, offset, max(0, offset), n))):
                extenders.append((observation, offset))
        if not extenders:
            done.append((sequence, positions))
            continue
        placed = {p for p in positions if p is not None}
        bases, added, active, forks, stopped, dropped = [], [], extenders, None, False, set()
        while True:
            votes = _votes(active, n + len(bases))
            if not votes or (len(votes) > 1 and bases):
                # Reads ending within this step leave fewer voters; decide a
                # disagreement only after re-collecting every overlapping read.
                break
            if len(votes) > 1:
                support = {key: len({o.fragment for o, _ in v}) for key, v in votes.items()}
                strong = _supported(support, *thresholds)
                if strong and len(strong) < len(votes):
                    anchor = next((p for p in reversed(positions + tuple(added)) if p is not None), None)
                    for key in sorted(set(votes) - strong, key=repr):
                        base, position = key
                        pruned.append(dict(
                            direction="5prime" if mirrored else "3prime",
                            after=None if anchor is None else list(_flip(anchor) if mirrored else anchor),
                            base=reverse_complement(base) if mirrored else base,
                            placement=None if position is None else list(_flip(position) if mirrored else position),
                            fragments=support[key], observations=sorted({o.key for o, _ in votes[key]})))
                        dropped.update(o.key for o, _ in votes[key])
                    votes = {key: votes[key] for key in strong}
                if len(votes) > 1:
                    forks = votes
                    break
            ((base, position), active), = votes.items()
            if position is not None:
                if position in placed:
                    # A path may not revisit a placement (e.g. around a tandem
                    # duplication); copy number cannot be read from overlaps.
                    notes.add("repeated_genomic_position")
                    stopped = True
                    break
                placed.add(position)
            bases.append(base)
            added.append(position)
        sequence, positions = sequence + "".join(bases), positions + tuple(added)
        agreeing = {o.key: o for o, _ in extenders if o.key not in dropped}
        if stopped or (forks and len(done) + len(stack) + len(forks) > budget):
            if forks:
                notes.add("path_limit")
            done.append((sequence, positions))
        elif forks:
            for key in sorted(forks, key=repr):
                chosen = {o.key for o, _ in forks[key]}
                rejected = {o.key for other, v in forks.items() if other != key for o, _ in v} - chosen
                branch = dict(used)
                for o in agreeing.values():
                    if o.key not in rejected:
                        branch.setdefault(o.identity, o.key)
                stack.append((sequence + key[0], positions + (key[1],), branch, forks[key]))
        else:
            used = dict(used)
            for o in agreeing.values():
                used.setdefault(o.identity, o.key)
            stack.append((sequence, positions, used, tuple(active) if bases else ()))
    return done


def _donor_distance(position, breakpoint):
    """Retained donor bases between this base and the breakpoint, if retained."""
    contig, base, strand = position
    if (contig, strand) != (breakpoint.contig, breakpoint.strand):
        return None
    distance = breakpoint.position - base - 1 if strand == "+" else base - breakpoint.position
    return distance if distance >= 0 else None


def _acceptor_distance(position, breakpoint):
    contig, base, strand = position
    if (contig, strand) != (breakpoint.contig, breakpoint.strand):
        return None
    distance = base - breakpoint.position if strand == "+" else breakpoint.position - base - 1
    return distance if distance >= 0 else None


def _forward_splice(left, right):
    """Whether the unrearranged reference could make this join by cis-splicing."""
    return left[::2] == right[::2] and (right[1] > left[1] if left[2] == "+" else right[1] < left[1])


def _event_relation(left, right, unplaced, donor, acceptor):
    """Compare one RNA join with the oriented DNA adjacency.

    A breakpoint junction may assign unplaced junction bases (homology or
    untemplated insertion) to either partner; the assignment is reported and
    is not verified against a genome sequence. A non-breakpoint join crossing
    from the donor to the acceptor side is splice-ambiguous when ordinary
    forward splicing (or read-through) of the reference could also make it.
    """
    d, a = _donor_distance(left, donor), _acceptor_distance(right, acceptor)
    if d is None or a is None:
        return None, None
    if d + a <= unplaced:
        return "breakpoint_junction", [d, a]
    if _forward_splice(left, right):
        return "splice_ambiguous_event_junction", [d, a]
    return "event_compatible_junction", [d, a]


class _Model:
    """A reference transcript's per-base placements, offsets and protein."""

    def __init__(self, reference):
        self.reference = reference
        positions = [(reference.contig, p, reference.strand) for a, b in reference.exons for p in range(a, b)]
        self.positions = positions[::-1] if reference.strand == "-" else positions
        self.offsets = {p: i for i, p in enumerate(self.positions)}
        self.protein = None
        if reference.cds_start is not None:
            self.protein = standard_genetic_code.translate(
                reference.sequence[reference.cds_start:reference.cds_end], first_codon_is_start=True)[0]

    def junctions(self):
        pairs = set()
        for a, b in zip(self.positions, self.positions[1:]):
            if b != _successor(a):
                pairs.update({(a, b), (_flip(b), _flip(a))})
        return pairs


def _anchor_runs(sequence, positions, model):
    """Exact, collinear matches to one transcript: (query start, end, transcript offset)."""
    runs = []
    for q, (base, position) in enumerate(zip(sequence, positions)):
        offset = model.offsets.get(position)
        if offset is None or model.reference.sequence[offset] != base:
            continue
        if runs and runs[-1][1] == q and runs[-1][2] + q - runs[-1][0] == offset:
            runs[-1] = (runs[-1][0], q + 1, runs[-1][2])
        else:
            runs.append((q, q + 1, offset))
    return runs


def _translations(sequence, positions, model, min_anchor, departures):
    """Readings from each coding anchor which leave the transcript at a junction.

    The frame is transferred downstream only; sequence after the anchor need
    not match any annotation. A base substitution keeps the frame and does not
    end the collinear match. A reading which first departs from the model
    anywhere other than ``departures`` (a CIGAR indel, sequencing error or
    another model's splice junction) is counted, not translated. Readings
    ending at the same codon share a frame; the longest is kept.
    """
    reference = model.reference
    runs = [r for r in _anchor_runs(sequence, positions, model) if r[1] - r[0] >= min_anchor]
    if not runs:
        return None
    if reference.cds_start is None:
        return dict(reason="CDS_unavailable", anchors=[list(r) for r in runs])
    readings, reference_reading, coding, other_departures = {}, False, [], 0
    for q0, q1, o0 in runs:
        lo = max(q0, q0 + reference.cds_start - o0)
        hi = min(q1, q0 + reference.cds_end - 3 - o0)
        if hi - lo < min_anchor:
            continue
        coding.append([q0, q1, o0])
        start = lo + (reference.cds_start - (o0 + lo - q0)) % 3
        offset = o0 + start - q0
        protein, stop = standard_genetic_code.translate(
            sequence[start:], first_codon_is_start=offset == reference.cds_start)
        end = start + 3 * (len(protein) + stop)
        departure, substitutions = None, []
        for q in range(start, end):
            expected = offset + q - start
            if expected >= len(model.positions) or positions[q] not in (None, model.positions[expected]):
                departure = q
            elif sequence[q] != reference.sequence[expected]:
                if positions[q] is None:
                    departure = q
                else:
                    substitutions.append(q)
            if departure is not None:
                break
        if departure is None:
            reference_reading = True
            continue
        if departure not in departures:
            other_departures += 1
            continue
        reading = readings.setdefault(end, dict(
            translation_start=start, translation_end=end, amino_acids=protein, ends_with_stop_codon=stop,
            complete_5prime=offset == reference.cds_start, departure=departure,
            substitutions=substitutions, anchors=[]))
        if start < reading["translation_start"]:
            reading.update(translation_start=start, amino_acids=protein, substitutions=substitutions,
                           complete_5prime=offset == reference.cds_start)
        reading["anchors"].append([q0, q1, o0])
    if not coding:
        return dict(reason="no_coding_anchor", anchors=[list(r) for r in runs])
    return dict(readings=sorted(readings.values(), key=lambda r: r["translation_start"]),
                reference_reading=reference_reading, other_departures=other_departures,
                anchors=[list(r) for r in runs])


def _frame_evidence(sequence, positions, departures, models, min_anchor, reference_peptides, lengths):
    """Group translations across all anchored models; keep every alternative."""
    groups, unresolved, reference_readings, other_departures = {}, [], [], 0
    for model in models:
        result = _translations(sequence, positions, model, min_anchor, departures)
        if result is None:
            continue
        reference = model.reference
        if "reason" in result:
            unresolved.append(dict(transcript_id=reference.transcript_id, annotation=reference.annotation,
                                   reason=result["reason"], anchors=result["anchors"]))
            continue
        if result["reference_reading"]:
            reference_readings.append(reference.transcript_id)
        other_departures += result["other_departures"]
        for reading in result["readings"]:
            key = (reading["translation_start"], reading["amino_acids"], reading["ends_with_stop_codon"])
            if key not in groups:
                protein = reading["amino_acids"]
                groups[key] = dict(
                    translation_start=reading["translation_start"], translation_end=reading["translation_end"],
                    amino_acids=protein, ends_with_stop_codon=reading["ends_with_stop_codon"],
                    transcript_ids=[], frame_evidence=[], competing_noncoding_models=[],
                    candidate_peptides=[dict(sequence=protein[i:i + k], protein_interval=[i, i + k])
                                        for k in lengths for i in range(len(protein) - k + 1)
                                        if protein[i:i + k] not in reference_peptides],
                    peptide_novelty="absent_from_supplied_reference_proteins", translation_observed=False)
            groups[key]["transcript_ids"].append(reference.transcript_id)
            groups[key]["frame_evidence"].append(dict(
                transcript_id=reference.transcript_id, annotation=reference.annotation,
                basis="exact_collinear_CDS_anchor", anchors=reading["anchors"],
                departure=reading["departure"], departure_relation=departures[reading["departure"]],
                frame_preserving_substitutions=reading["substitutions"],
                complete_5prime=reading["complete_5prime"],
                upstream_frame_assumed=not reading["complete_5prime"]))
    translations = sorted(groups.values(), key=lambda g: (g["translation_start"], g["amino_acids"]))
    for translation in translations:
        translation["departure_relations"] = sorted({e["departure_relation"] for e in translation["frame_evidence"]})
        # A noncoding model matching the same frame-anchor sequence is an
        # alternative explanation, not permission to drop the protein.
        spans = [a for e in translation["frame_evidence"] for a in e["anchors"]]
        translation["competing_noncoding_models"] = sorted(
            u["transcript_id"] for u in unresolved
            if any(a[0] < s[1] and s[0] < a[1] for a in u["anchors"] for s in spans))
    if not translations:
        status = ("reference_protein_only" if reference_readings else "unresolved" if unresolved
                  else "no_gene_anchor")
    elif len(translations) == 1 and not translations[0]["competing_noncoding_models"]:
        status = "translated"
    else:
        status = "ambiguous"
    return dict(frame_status=status, translations=translations, unresolved_models=unresolved,
                reference_readings=sorted(reference_readings), readings_departing_elsewhere=other_departures)


def _seeds(observations, models, annotated, donor, acceptor, min_fragments):
    """Unannotated junctions (and breakpoint clips) at which paths start."""
    spans = defaultdict(list)
    for model in models:
        r = model.reference
        spans[r.contig, r.strand].append((r.exons[0][0], r.exons[-1][1]))

    def anchored(position):
        return any(a <= position[1] < b for a, b in spans.get((position[0], position[2]), ()))

    seeds, competing = {}, defaultdict(set)
    for o in observations:
        cores = []
        for i, j, kind in o.breaks:
            left, right = o.positions[i], o.positions[j]
            relation, assignment = _event_relation(left, right, j - i - 1, donor, acceptor)
            if j == i + 1 and (left, right) in annotated:
                if relation:
                    competing[left, right].add(o.fragment)
                continue
            if relation is None:
                if kind == "D" or not (anchored(left) or anchored(right)):
                    continue
                relation = "regional_novel_junction"
            cores.append((i, j, relation, kind))
        placed = [q for q, p in enumerate(o.positions) if p is not None]
        if placed[-1] < len(o.sequence) - 1 and _donor_distance(o.positions[placed[-1]], donor) == 0:
            cores.append((placed[-1], placed[-1] + 1, "breakpoint_clip_partner_unplaced", "clip"))
        if placed[0] > 0 and _acceptor_distance(o.positions[placed[0]], acceptor) == 0:
            cores.append((placed[0] - 1, placed[0], "breakpoint_clip_partner_unplaced", "clip"))
        for i, j, relation, kind in cores:
            seed = seeds.setdefault((o.sequence[i:j + 1], o.positions[i:j + 1]), dict(
                relation=relation, fragments=set(), observations=[], kinds=set()))
            seed["fragments"].add(o.fragment)
            seed["observations"].append((o, i))
            seed["kinds"].add(kind)
    rank = {status: i for i, status in enumerate(LINKAGE_STATUSES)}
    ordered = sorted(
        ((key, seed) for key, seed in seeds.items()
         if seed["relation"] not in _SUPPORT_GATED or len(seed["fragments"]) >= min_fragments),
        key=lambda item: (rank[item[1]["relation"]], -len(item[1]["fragments"]), item[0][0], repr(item[0][1])))
    return ordered, [dict(left=list(left), right=list(right), fragments=len(fragments))
                     for (left, right), fragments in sorted(competing.items())]


def _contains(path, core):
    sequence, positions = path
    start = sequence.find(core[0])
    while start >= 0:
        if all(x is None or x == y for x, y in zip(core[1], positions[start:start + len(core[0])])):
            return True
        start = sequence.find(core[0], start + 1)
    return False


def _support(sequence, positions, index):
    """Observations contained in, and agreeing with, the whole path."""
    result = []
    for key, offset in sorted(index.overlapping(positions, 0, len(positions))):
        observation = index.observations[key]
        if (0 <= offset and offset + len(observation.sequence) <= len(sequence)
                and _agrees(sequence, positions, observation, offset, offset, offset + len(observation.sequence))):
            result.append((observation, offset))
    return result


def _path_result(sequence, positions, index, donor, acceptor, annotated, frame_evidence):
    """Serialize one path with junction, sequence and frame evidence.

    Returns the result and the observations reported as direct evidence:
    those spanning a reported junction or breakpoint clip, or the whole path.
    """
    support = _support(sequence, positions, index)
    kinds = defaultdict(set)
    for observation, offset in support:
        for i, j, kind in observation.breaks:
            kinds[i + offset, j + offset].add(kind)
    rows = []
    for i, j in _breaks(positions):
        left, right = positions[i], positions[j]
        relation, assignment = _event_relation(left, right, j - i - 1, donor, acceptor)
        is_annotated = j == i + 1 and (left, right) in annotated
        # Annotated splicing and CIGAR deletions are reported only when they
        # also cross the nominated event.
        if relation is None and (is_annotated or kinds[i, j] <= {"D"}):
            continue
        rows.append((i, j, relation or "regional_novel_junction", assignment, is_annotated))
    placed = [q for q, p in enumerate(positions) if p is not None]
    if placed[-1] < len(sequence) - 1 and _donor_distance(positions[placed[-1]], donor) == 0:
        rows.append((placed[-1], placed[-1] + 1, "breakpoint_clip_partner_unplaced", [0, None], False))
    if placed[0] > 0 and _acceptor_distance(positions[placed[0]], acceptor) == 0:
        rows.append((placed[0] - 1, placed[0], "breakpoint_clip_partner_unplaced", [None, 0], False))
    junctions, evidence = [], {}
    for i, j, relation, assignment, is_annotated in sorted(rows, key=lambda row: row[:2]):
        spanning = [(o, offset) for o, offset in support if offset <= i and j < offset + len(o.sequence)]
        evidence.update((o.key, o) for o, _ in spanning)
        junctions.append(dict(
            query_interval=[i, j], left=positions[i] and list(positions[i]),
            right=positions[j] and list(positions[j]), unplaced_bases=sequence[i + 1:j],
            forward_splice_geometry=bool(positions[i] and positions[j] and _forward_splice(positions[i], positions[j])),
            kinds=sorted(kinds[i, j]) or (["clip"] if relation == "breakpoint_clip_partner_unplaced" else []),
            annotated=is_annotated, relation=relation, breakpoint_assignment=assignment,
            spanning_segments=len({o.identity for o, _ in spanning}),
            spanning_fragments=len({o.fragment for o, _ in spanning}),
            spanning_observations=sorted({o.key for o, _ in spanning}),
            linked_interval=([min(offset for _, offset in spanning),
                              max(offset + len(o.sequence) for o, offset in spanning)] if spanning else None)))
    departures = {q: j["relation"] for j in junctions if not j["annotated"]
                  for q in range(j["query_interval"][0] + 1,
                                 len(sequence) if j["right"] is None else j["query_interval"][1] + 1)}
    frame = frame_evidence(sequence, positions, departures)
    relations = {j["relation"] for j in junctions if not j["annotated"]}
    status = next((s for s in LINKAGE_STATUSES if s in relations), "no_novel_junction")
    whole = sorted({o.key for o, offset in support if offset == 0 and len(o.sequence) == len(sequence)})
    evidence.update((key, index.observations[key]) for key in whole)
    segments = {o.identity for o, _ in support}
    return dict(
        path_id=sha256(repr((sequence, positions)).encode()).hexdigest()[:24],
        sequence=sequence,
        blocks=[asdict(FusionBlock(a, b, contig, low, high, strand)) for a, b, contig, low, high, strand in _runs(positions)],
        unplaced_intervals=_intervals(q for q, p in enumerate(positions) if p is None),
        junctions=junctions,
        event_linkage=dict(status=status, junctions=[k for k, j in enumerate(junctions)
                                                     if j["relation"] == status and not j["annotated"]],
                           somatic_causation_proven=False),
        sequence_evidence=dict(
            contained_segments=len(segments), contained_fragments=len({s[:2] for s in segments}),
            whole_path_observations=whole,
            missing_quality_segments=len({o.identity for o, _ in support if o.missing_qualities}),
            secondary_segments=len({o.identity for o, _ in support if o.secondary}),
            assembly_phase="hypothesis_not_proven_long_range_phase"),
        **frame), evidence


def _intervals(indices):
    result = []
    for q in indices:
        if result and result[-1][1] == q:
            result[-1][1] = q + 1
        else:
            result.append([q, q + 1])
    return result


def reconstruct_sv_rna(bam, *, event_id, reference_name, donor, acceptor, regions, references,
                       sample_id, source, event_provenance, read_collector=None,
                       min_anchor_bases=SV_MIN_ANCHOR_BASES, min_overlap=SV_MIN_OVERLAP,
                       min_alternative_fragments=SV_MIN_ALTERNATIVE_FRAGMENTS,
                       min_alternative_fraction=SV_MIN_ALTERNATIVE_FRACTION, max_records=SV_MAX_RECORDS, max_queries=SV_MAX_QUERIES,
                       max_paths=SV_MAX_PATHS, assemble=SV_ASSEMBLE,
                       breakpoint_window=SV_BREAKPOINT_WINDOW, peptide_lengths=FUSION_PEPTIDE_LENGTHS):
    """Reconstruct RNA paths for one nominated, oriented DNA adjacency.

    Paths start at unannotated RNA junctions (and, when soft-clipped bases are
    used, unaligned sequence at a breakpoint) and extend through exact overlaps
    anchored by shared placements. With ``assemble=False`` only observations
    containing the seed junction contribute, so every base is observed in a
    read which also spans that junction. Frames are transferred from exact,
    collinear annotated CDS anchors into downstream observed sequence.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Indexed alignments from one sample/library source. Missing QUAL is kept.
    donor, acceptor : FusionBreakpoint
        Nominated oriented adjacency, retained-partner boundaries as in
        ``FusionTranscript``. It need not be an RNA junction.
    regions : iterable of (str, int, int)
        Additional genomic search intervals. Annotated exons and
        ``breakpoint_window`` bases around each breakpoint are always searched.
    references : iterable of FusionReference
        All candidate transcript models. Peptide novelty is relative to them.
    sample_id, source : str
        Sample and alignment-source identity.
    event_provenance : dict
        Caller-supplied DNA event evidence, retained verbatim.
    read_collector : ReadCollector, optional
        Read filters, soft-clip use and read-end handling.
    min_anchor_bases : int
        Exact collinear transcript match required to transfer a frame.
    min_overlap : int
        Overlap required to extend a path with another observation.
    min_alternative_fragments : int
        Fragments needed for an unattributed junction to seed a path, and for
        an extension branch to outweigh (and prune) a weaker alternative.
    min_alternative_fraction : float
        An alternative extension branch, or a sequence variant of an already
        seeded junction, is also pruned when its fragments are below this
        fraction of the best-supported alternative. Pruning is reported.
    max_records, max_queries, max_paths : int
        Resource limits; reaching one is reported in ``limitations``.
    assemble : bool
        Extend through overlapping observations which do not span the seed.
    breakpoint_window : int
        Bases searched to either side of each nominated breakpoint.
    peptide_lengths : iterable of int

    Returns
    -------
    dict
        JSON-serializable exploratory candidates with independent sequence,
        frame and event-linkage evidence, original records and limitations.
        This is not the validated supplied-fusion result consumed by Vaxrank.
    """
    if not all((event_id, reference_name, sample_id, source)) or not isinstance(event_provenance, dict) \
            or not event_provenance:
        raise ValueError("SV reconstruction requires event, reference, sample, source and event provenance")
    if not isinstance(donor, FusionBreakpoint) or not isinstance(acceptor, FusionBreakpoint):
        raise TypeError("Expected oriented FusionBreakpoint inputs")
    for value in (min_anchor_bases, min_overlap, min_alternative_fragments, max_records, max_queries, max_paths,
                  breakpoint_window):
        if type(value) is not int or value < 1:
            raise ValueError("SV evidence thresholds and search limits must be positive integers")
    if isinstance(min_alternative_fraction, bool) or not 0 <= min_alternative_fraction <= 1:
        raise ValueError("min_alternative_fraction must be between 0 and 1")
    if min_anchor_bases < 3:
        raise ValueError("A frame anchor needs at least one codon")
    lengths = tuple(sorted(set(peptide_lengths)))
    if not lengths or any(type(n) is not int or n < 1 for n in lengths):
        raise ValueError("Peptide lengths must be positive integers")
    references = tuple(references)
    if not references or any(not isinstance(r, FusionReference) or r.reference_name != reference_name
                             for r in references):
        raise ValueError("Matching, explicitly versioned reference models are required")
    if len({(r.annotation, r.transcript_id) for r in references}) != len(references):
        raise ValueError("Duplicate reference transcript identity")
    collector = ReadCollector() if read_collector is None else read_collector
    models = [_Model(r) for r in references]
    annotated = set().union(*(m.junctions() for m in models))
    reference_peptides = {m.protein[i:i + k] for m in models if m.protein
                          for k in lengths for i in range(len(m.protein) - k + 1)}

    records, acquisition = collect_sv_records(bam, regions, references, max_records, max_queries,
                                              (donor, acceptor), breakpoint_window)
    observations, excluded, reasons = segment_observations(records, collector)
    index = _ObservationIndex(observations)
    seeds, competing = _seeds(observations, models, annotated, donor, acceptor, min_alternative_fragments)

    thresholds = (min_alternative_fragments, min_alternative_fraction)
    pairs = list(zip(observations[::2], observations[1::2]))  # Both orientations of each segment path.
    mirrors = {**{a.key: b for a, b in pairs}, **{b.key: a for a, b in pairs}}
    paths, pruned, notes, seed_rows, strongest = {}, [], set(acquisition["limitations"]), [], {}
    for (core_sequence, core_positions), seed in seeds:
        row = dict(relation=seed["relation"], kinds=sorted(seed["kinds"]),
                   left=list(core_positions[0]) if core_positions[0] else None,
                   right=list(core_positions[-1]) if core_positions[-1] else None,
                   unplaced_bases=core_sequence[1:-1], fragments=len(seed["fragments"]),
                   segments=len({o.identity for o, _ in seed["observations"]}), paths=[])
        seed_rows.append(row)
        # Seeds are ordered by support: a weak base variant of an already
        # seeded junction is pruned like any other weak extension branch.
        if core_positions in strongest:
            strong = _supported({"best": strongest[core_positions], "this": row["fragments"]}, *thresholds)
            if "best" in strong and "this" not in strong:
                row["minor_variant_of_seeded_junction"] = True
                continue
        strongest.setdefault(core_positions, row["fragments"])
        if any(_contains(path, (core_sequence, core_positions)) for path in paths):
            row["represented_by_earlier_path"] = True
            continue
        if len(paths) >= max_paths:
            notes.add("path_limit")
            row["unexplored"] = True
            continue
        pool = index
        if not assemble:
            spanning = [o for o, _ in seed["observations"]]
            pool = _ObservationIndex(spanning + [mirrors[o.key] for o in spanning])
        for sequence, positions in _extend(core_sequence, core_positions, pool, min_overlap,
                                           thresholds, max_paths - len(paths), pruned, notes, False):
            for mirror in _extend(*_mirror(sequence, positions), pool, min_overlap, thresholds,
                                  max(1, max_paths - len(paths)), pruned, notes, True):
                path = _mirror(*mirror)
                paths.setdefault(path, [])
                if len(paths) > max_paths:
                    del paths[path]
                    notes.add("path_limit")
                    break
                paths[path].append(row)
    results, used = [], {}
    for (sequence, positions), rows in paths.items():
        result, evidence = _path_result(
            sequence, positions, index, donor, acceptor, annotated,
            lambda *path: _frame_evidence(*path, models, min_anchor_bases, reference_peptides, lengths))
        results.append(result)
        for row in rows:
            if result["path_id"] not in row["paths"]:
                row["paths"].append(result["path_id"])
        used.update(evidence)
    rank = {status: i for i, status in enumerate(LINKAGE_STATUSES)}
    results.sort(key=lambda r: (rank.get(r["event_linkage"]["status"], len(rank)),
                                -max([j["spanning_fragments"] for j in r["junctions"]], default=0),
                                r["sequence"], r["path_id"]))
    record_ids = {rid for o in used.values() for rid in o.records}
    statuses = {r["event_linkage"]["status"] for r in results}
    return dict(
        schema="isovar.sv_rna_candidates.v1",
        status=("event_linked_candidates" if statuses - _SUPPORT_GATED - {"no_novel_junction"}
                else "regional_candidates_only" if results else "no_candidate_paths"),
        event_id=event_id, reference_name=reference_name, sample_id=sample_id, source=source,
        event_provenance=event_provenance, donor=asdict(donor), acceptor=asdict(acceptor),
        reference_models=[dict(transcript_id=r.transcript_id, annotation=r.annotation, contig=r.contig,
                               strand=r.strand, sequence_sha256=sha256(r.sequence.encode()).hexdigest())
                          for r in references],
        acquisition=acquisition, limitations=sorted(notes),
        parameters=dict(min_anchor_bases=min_anchor_bases, min_overlap=min_overlap,
                        min_alternative_fragments=min_alternative_fragments,
                        min_alternative_fraction=min_alternative_fraction, max_records=max_records,
                        max_queries=max_queries, max_paths=max_paths, assemble=assemble,
                        breakpoint_window=breakpoint_window, peptide_lengths=lengths,
                        genetic_code=1, use_soft_clipped_bases=collector.use_soft_clipped_bases,
                        min_mapping_quality=collector.min_mapping_quality),
        seeds=seed_rows, competing_annotated_junctions=competing, pruned_branches=pruned, paths=results,
        observations={key: dict(identity=list(o.identity), records=list(o.records), reverse_complement=o.reverse,
                                original_query_interval=list(o.query_interval),
                                missing_qualities=o.missing_qualities, secondary=o.secondary)
                      for key, o in sorted(used.items())},
        excluded_records=dict(Counter(excluded.values())), segment_path_notes=dict(reasons),
        original_records={_record_id(r): r.to_string() for r in records if _record_id(r) in record_ids})


def sv_rna_input_from_dict(data):
    """Decode the explicit JSON input used by ``isovar sv-rna``.

    Returns keyword arguments for ``reconstruct_sv_rna`` (other than the BAM,
    source and evidence parameters).
    """
    return dict(
        event_id=data["event_id"], reference_name=data["reference_name"], sample_id=data["sample_id"],
        donor=FusionBreakpoint(**data["donor"]), acceptor=FusionBreakpoint(**data["acceptor"]),
        regions=[tuple(region) for region in data.get("regions", [])],
        references=[FusionReference(**r) for r in data["references"]],
        event_provenance=data["event_provenance"])
