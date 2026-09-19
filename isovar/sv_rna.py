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
from heapq import nsmallest
from dataclasses import asdict, dataclass
from hashlib import sha256
from itertools import combinations
from types import SimpleNamespace

from .chimeric_alignment import compatible_phasing_alignments, source_alignment_paths_from_pysam
from .default_parameters import (
    FUSION_PEPTIDE_LENGTHS, SV_ASSEMBLE, SV_BREAKPOINT_WINDOW, SV_MAX_BREAKPOINT_SHIFT,
    SV_MAX_EXTENSION_SEGMENTS, SV_MAX_PATHS, SV_MAX_QUERIES, SV_MAX_RECORDS, SV_ANNOTATED_JUNCTION_TOLERANCE,
    SV_MIN_ALTERNATIVE_FRACTION, SV_MIN_ALTERNATIVE_FRAGMENTS, SV_MIN_LOCAL_VARIANT_FRACTION, SV_MIN_ANCHOR_BASES, SV_MIN_OVERLAP,
)
from .fusion import FusionBlock, FusionBreakpoint, FusionReference
from .genetic_code import standard_genetic_code
from .read_collector import ReadCollector
from .read_end_inference import reverse_complement
from .read_identity import source_alignments_from_pysam

_BIN = 64  # Genomic bin width of the observation index.
_RECORD_BIN = 1024  # Genomic bin width of the unbuilt-record index.
_MIN_INTRON = 21  # Shorter CIGAR N gaps are deletions (STAR's alignIntronMin default).
_LOCAL = 20  # Placements this close are local (base/indel) alternatives, not splice choices.
_HOP_WINDOW = 1000  # Nearby SA/mate targets share one indexed query.
_MAX_SEGMENT_PATHS = 64  # Alternative supplementary groupings per segment.

# Linkage of a path to the nominated adjacency, strongest first. The first
# three are event-linked; the last two may arise without the rearrangement.
LINKAGE_STATUSES = (
    "breakpoint_junction",  # The RNA join reproduces the DNA adjacency.
    "event_compatible_junction",  # Donor-side to acceptor-side join, e.g. spliced.
    "breakpoint_clip_partner_unplaced",  # Unaligned sequence at a breakpoint.
    # Crosses the event, but ordinary forward splicing (or read-through) of
    # the unrearranged reference could make it: any such non-breakpoint join,
    # and a breakpoint join between annotated splice sites.
    "splice_ambiguous_event_junction",
    "regional_novel_junction",  # Unannotated join not crossing the event.
)


def segment_identity(read):
    """(read group, QNAME, mate bits): one sequenced segment in one library."""
    return source_alignments_from_pysam(read, read.query_name)[0][0]


def _record_id(read):
    """Stable record ID; hashing the key avoids hashing long reads' SEQ/QUAL."""
    return sha256(repr(_record_key(read)).encode()).hexdigest()


def _record_key(read):
    """Cheap identity of one SAM record, for deduplicating repeated fetches."""
    return (read.query_name, read.flag, read.reference_id, read.reference_start, read.cigarstring or "",
            read.get_tag("RG") if read.has_tag("RG") else "")


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
                key = _record_key(read)
                if key in records:
                    continue
                if len(records) >= max_records:
                    limitations.add("record_limit")
                    break
                records[key] = read
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
    return [records[key] for key in sorted(records)], dict(
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


def _record_placements(read, view, intern):
    """Actual CIGAR placements and gap kinds, in original sequencing coordinates."""
    length = read.infer_read_length()
    hard = read.cigartuples[0][1] if read.cigartuples[0][0] == 5 else 0
    start, end = view.sequenced_interval(0, len(view.sequence))
    strand = "-" if read.is_reverse else "+"

    def original(q):
        return length - 1 - hard - q if read.is_reverse else hard + q

    placements, kinds, previous = {}, {}, None
    for q, r, n, gap, _ in _aligned_blocks(read):
        for k in range(n):
            if start <= original(q + k) < end:
                position = (read.reference_name, r + k, strand)
                placements[original(q + k)] = intern.setdefault(position, position)
        if previous is not None and gap is not None:
            kinds[tuple(sorted((original(previous), original(q))))] = _gap_kind(*gap)
        previous = q + n - 1
    return placements, kinds


def _segment_path(records, collector, intern):
    """Join one segment's linked records by original query offset."""
    bases, placements, sources, kinds = {}, defaultdict(set), defaultdict(set), {}
    for index, read in enumerate(records):
        view = collector.read_sequence_view(read)
        start, _ = view.sequenced_interval(0, len(view.sequence))
        sequence = (reverse_complement(view.sequence) if read.is_reverse else view.sequence).upper()
        for q, base in enumerate(sequence, start):
            if bases.setdefault(q, base) != base:
                return "conflicting_segment_sequence"
        mapping, gaps = _record_placements(read, view, intern)
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


def _build_segment(identity, records, collector, reasons, intern):
    """Both orientations of each linked grouping of one segment's records.

    ``intern`` shares one tuple per placement: deep reads cover the same bases.
    """
    observations = []
    for members in _segment_groupings(sorted(records, key=_record_id), reasons):
        if len(records) > 1 and len(members) == 1 and members[0].has_tag("SA"):
            reasons["SA_declared_piece_not_linked"] += 1
        built = _segment_path(members, collector, intern)
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
        positions = tuple(None if p is None else intern.setdefault(p, p) for p in positions)
        observations.append(RnaObservation(
            key + "-", sequence=sequence, positions=positions, reverse=True,
            breaks=tuple(sorted((n - 1 - j, n - 1 - i, kind) for i, j, kind in breaks)), **common))
    return observations


def _gap_kind(spliced, length):
    """A reference gap between aligned bases: N only for a plausible intron."""
    return "N" if spliced and length >= _MIN_INTRON else "D"


def _aligned_blocks(read):
    """(query offset, reference start, length, preceding gap, inserted bases) per aligned block.

    Query offsets index the stored SEQ. The preceding gap is ``(spliced,
    reference length)`` of the D/N operations since the previous block, or
    ``None``; inserted bases are the I operations since then.
    """
    q, r, gap, inserted = 0, read.reference_start, None, 0
    for operation, n in read.cigartuples:
        if operation in (0, 7, 8):
            yield q, r, n, gap, inserted
            q, r, gap, inserted = q + n, r + n, None, 0
        elif operation in (1, 4):
            q += n
            inserted += n if operation == 1 else 0
        elif operation in (2, 3):
            spliced, size = gap or (False, 0)
            gap, r = (spliced or operation == 3, size + n), r + n


def _record_gaps(read):
    """(left, right, kind, inserted sequence) of CIGAR N/D joins, in + orientation."""
    gaps, last = [], None
    for q, r, n, gap, inserted in _aligned_blocks(read):
        if last is not None and gap is not None:
            gaps.append(((read.reference_name, last, "+"), (read.reference_name, r, "+"), _gap_kind(*gap),
                         read.query_sequence[q - inserted:q].upper()))
        last = r + n - 1
    return gaps


def _offsets(query_runs, lo, key, runs):
    """(key, path offset of the observation start) for each shared placed base run."""
    found = set()
    for start, _, contig, low, high, strand in query_runs:
        for other, _, other_contig, other_low, other_high, other_strand in runs:
            if (other_contig, other_strand) != (contig, strand) or other_low >= high or low >= other_high:
                continue
            if strand == "+":
                shared = low if low > other_low else other_low
                path_q, observed_q = start + shared - low, other + shared - other_low
            else:
                shared = (high if high < other_high else other_high) - 1
                path_q, observed_q = start + high - 1 - shared, other + other_high - 1 - shared
            found.add((key, lo + path_q - observed_q))
    return found


class _Observations:
    """Segment observations, built only when a seed or a path end needs them.

    Building per-base placements for every record scales with total aligned
    bases, which deep long-read loci make prohibitive. Records are instead
    grouped by segment and binned by aligned block. A path end builds at most
    ``cap`` overlapping segments, those reaching furthest in its direction;
    reaching the cap is reported.
    """

    def __init__(self, records, collector, cap):
        self.collector, self.cap, self.limited = collector, cap, False
        self.groups, self.excluded, self.reasons = defaultdict(list), Counter(), Counter()
        self.extent, self.bins, self.gaps = defaultdict(dict), defaultdict(list), defaultdict(dict)
        self.gap_kinds = defaultdict(set)
        self.built, self.observations, self.joins, self.runs = {}, {}, defaultdict(list), {}
        self.intern, self.molecules, self.windows = {}, {}, {}
        for read in records:
            reason = _ineligible(read, collector)
            if reason:
                self.excluded[reason] += 1
                continue
            identity = segment_identity(read)
            self.groups[identity].append(read)
            contig = read.reference_name
            for start, end in read.get_blocks():
                for b in range(start // _RECORD_BIN, (end - 1) // _RECORD_BIN + 1):
                    self.bins[contig, b].append((start, end, identity))
            low, high = self.extent[identity].get(contig, (read.reference_start, read.reference_end))
            self.extent[identity][contig] = (min(low, read.reference_start), max(high, read.reference_end))

    def build(self, identity):
        if identity not in self.built:
            self.built[identity] = _build_segment(identity, self.groups[identity], self.collector, self.reasons,
                                                  self.intern)
            for o in self.built[identity]:
                self.observations[o.key] = o
                for i, j, kind in o.breaks:
                    self.joins[o.positions[i], o.positions[j]].append((o, i, j, kind))
        return self.built[identity]

    def molecule(self, identity):
        """(cell barcode, UMI) of a segment, when its primary record carries both."""
        if identity not in self.molecules:
            self.molecules[identity] = self._molecule(identity)
        return self.molecules[identity]

    def observation_runs(self, observation):
        if observation.key not in self.runs:
            self.runs[observation.key] = _runs(observation.positions)
        return self.runs[observation.key]

    def _molecule(self, identity):
        read = min(self.groups[identity], key=lambda r: (r.is_supplementary or r.is_secondary))
        cell = read.get_tag("CB") if read.has_tag("CB") else None
        umi = next((read.get_tag(t) for t in ("UB", "XM") if read.has_tag(t)), None)
        return (cell, umi) if cell and umi else None

    def classify(self, event, annotated, anchored):
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
            for left, right, kind, inserted in _record_gaps(read):
                self.gaps[left, right][identity] = inserted
                self.gap_kinds[left, right].add(kind)
                oriented = [(left, right), (_flip(right), _flip(left))]
                classes = [_join_class(event, annotated, *join, len(inserted), kind) for join in oriented]
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
                    lazy[next(r for use, r in classes if use == "lazy")][_canonical(left, right)].add(identity)
                elif kind != "D" and all(r is None for _, r in classes) and any(
                        anchored(p) for p in (left, right, _flip(left), _flip(right))):
                    lazy["regional_novel_junction"][_canonical(left, right)].add(identity)
            # SAM permits soft clips inside terminal hard clips.
            if (self.collector.use_soft_clipped_bases and any(op == 4 for op, _ in read.cigartuples)
                    and any(read.reference_name == contig and abs(edge - position) <= 1
                            for contig, position in breakpoints for edge in (read.reference_start, read.reference_end))):
                now.add(identity)
        return now, lazy, competing

    def join_segments(self, left, right):
        """Segments whose own path makes this join: built ones plus never-built single records.

        A built segment without the join (a failed build, or a trimmed read
        end) does not count.
        """
        unbuilt = set(self.gaps.get((left, right), ())) | set(self.gaps.get((_flip(right), _flip(left)), ()))
        return {o.identity for o, *_ in self.joins.get((left, right), ())} | {i for i in unbuilt if i not in self.built}

    def overlapping(self, positions, lo, hi):
        """(key, path offset) of observations sharing a placed base in positions[lo:hi].

        Forked branches re-query shared windows; results are memoized.
        """
        window = positions[lo:hi]
        if window not in self.windows:
            self.windows[window] = self._overlapping(window)
        return {(key, lo + offset) for key, offset in self.windows[window]}

    def _overlapping(self, window):
        runs = _runs(window)
        covered = Counter()  # Window bases each segment's aligned blocks cover.
        for _, _, contig, low, high, _ in runs:
            blocks = {block for b in range(low // _RECORD_BIN, (high - 1) // _RECORD_BIN + 1)
                      for block in self.bins.get((contig, b), ()) if block[0] < high and low < block[1]}
            for start, end, identity in blocks:
                covered[identity] += min(end, high) - max(start, low)
        candidates = set(covered)
        if len(candidates) > self.cap:
            # Extenders must overlap the whole window (``min_overlap``): prefer
            # segments covering most of it (e.g. both sides of a junction),
            # then those reaching furthest beyond it.
            self.limited = True
            _, _, contig, _, _, strand = runs[-1]

            def rank(identity):
                low, high = self.extent[identity].get(contig, (float("inf"), float("-inf")))
                return -covered[identity], -(high if strand == "+" else -low), identity
            candidates = nsmallest(self.cap, candidates, key=rank)
        found = set()
        for identity in candidates:
            for o in self.build(identity):
                found |= _offsets(runs, 0, o.key, self.observation_runs(o))
        return found


class _ObservationIndex:
    """A fixed set of observations (e.g. one seed's reads), binned by placed runs."""

    def __init__(self, observations):
        self.observations = {o.key: o for o in observations}
        self.runs = {o.key: _runs(o.positions) for o in observations}
        self.bins = defaultdict(set)
        for key, runs in self.runs.items():
            for _, _, contig, low, high, strand in runs:
                for b in range(low // _BIN, (high - 1) // _BIN + 1):
                    self.bins[contig, strand, b].add(key)

    def overlapping(self, positions, lo, hi):
        """(key, path offset of the observation start) sharing a placed base in positions[lo:hi]."""
        runs = _runs(positions[lo:hi])
        keys = {key for _, _, contig, low, high, strand in runs
                for b in range(low // _BIN, (high - 1) // _BIN + 1) for key in self.bins.get((contig, strand, b), ())}
        return {found for key in keys for found in _offsets(runs, lo, key, self.runs[key])}


def _agrees(sequence, positions, observation, offset, lo, hi):
    """Identical bases and no conflicting placement over path[lo:hi]."""
    a, b = lo - offset, hi - offset
    if observation.sequence[a:b] != sequence[lo:hi]:
        return False
    mine, theirs = positions[lo:hi], observation.positions[a:b]
    return mine == theirs or all(x is None or y is None or x == y for x, y in zip(mine, theirs))


def _supported(support, min_fragments, min_fraction, min_local_fraction):
    """Alternatives which meet the thresholds; empty when none does.

    Keys are (base, placement). An alternative placed within ``_LOCAL`` bases
    of the best one (or unplaced) is a base or small-indel variant of it, and
    also needs ``min_local_fraction`` of the best's support: systematic
    long-read errors often exceed ``min_fraction``, while a heterozygous
    allele approaches the best. Distant placements (other splice choices)
    need only ``min_fraction``.
    """
    best_key = max(support, key=lambda key: (support[key], repr(key)))
    best = support[best_key]

    def local(key):
        a, b = key[1], best_key[1]
        return a is None or b is None or (a[::2] == b[::2] and abs(a[1] - b[1]) <= _LOCAL)

    return {key for key, count in support.items()
            if count >= min_fragments and count >= min_fraction * best
            and (key == best_key or not local(key) or count >= min_local_fraction * best)}


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
    soft-clipped tails) where no new read can be anchored. Each finished
    path is returned with the keys of the observations which voted for it.
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
            done.append((sequence, positions, frozenset(used.values())))
            continue
        placed = {p for p in positions if p is not None}
        bases, added, active, forks, stopped, dropped = [], [], extenders, None, False, set()
        while True:
            votes = _votes(active, n + len(bases))
            if not votes or (bases and (len(votes) > 1 or 2 * sum(map(len, votes.values())) < len(extenders))):
                # Reads ending within this step leave fewer voters: decide a
                # disagreement, or continue once most extenders have ended,
                # only after re-collecting every read overlapping the new end
                # (so one long read cannot decide the rest of the path alone).
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
            voters = dict(used)
            for o in agreeing.values():
                voters.setdefault(o.identity, o.key)
            done.append((sequence, positions, frozenset(voters.values())))
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


def _distance(position, breakpoint, side):
    """Signed retained-partner bases between a base and a breakpoint.

    Negative when the base lies past the breakpoint, on the side the adjacency
    removes; ``None`` on another contig or strand.
    """
    contig, base, strand = position
    if (contig, strand) != (breakpoint.contig, breakpoint.strand):
        return None
    return breakpoint.position - base - 1 if (side == "donor") == (strand == "+") else base - breakpoint.position


def _forward_splice(left, right):
    """Whether the unrearranged reference could make this join by cis-splicing."""
    return left[::2] == right[::2] and (right[1] > left[1] if left[2] == "+" else right[1] < left[1])


class _Event:
    """The nominated oriented adjacency, compared with observed RNA joins."""

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
            _forward_splice(left, right) and left in self.donor_sites and right in self.acceptor_sites
            for left, right in (self._at(d, -d) for d in range(-max_shift, max_shift + 1)))

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

    def relation(self, left, right, unplaced):
        """Relation of one join to the adjacency, and its breakpoint assignment.

        Unplaced junction bases (homology or untemplated insertion) may belong
        to either partner, and an aligner may place up to ``max_shift``
        homologous bases past a breakpoint: every such join on the adjacency's
        diagonal is the breakpoint junction. The homology is not verified
        against a genome sequence. Every placement of a read-through-ambiguous
        adjacency is splice-ambiguous. A join just off the diagonal, with one
        side past a breakpoint, is a (noisy) event join, not an unrelated one.
        """
        d, a = _distance(left, self.donor, "donor"), _distance(right, self.acceptor, "acceptor")
        if d is None or a is None:
            return None, None
        if self.on_diagonal(d, a, unplaced):
            return ("splice_ambiguous_event_junction" if self.read_through else "breakpoint_junction"), [d, a]
        if min(d, a) < 0 and not self.near_breakpoint(d, a, unplaced):
            return None, None
        return ("splice_ambiguous_event_junction" if _forward_splice(left, right)
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
        placed = [q for q, p in enumerate(positions) if p is not None]
        found = []
        if placed[-1] < length - 1 and _distance(positions[placed[-1]], self.donor, "donor") == 0:
            found.append((placed[-1], placed[-1] + 1, [0, None]))
        if placed[0] > 0 and _distance(positions[placed[0]], self.acceptor, "acceptor") == 0:
            found.append((placed[0] - 1, placed[0], [None, 0]))
        return found


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


def _canonical(left, right):
    """One key for a join and its reverse complement."""
    return min((left, right), (_flip(right), _flip(left)))


def _join_class(event, annotated, left, right, unplaced, kind):
    """How one join, seen in this orientation, is used.

    Returns ``("annotated", relation)``, ``("event", relation)`` for joins at
    or across the adjacency, ``("lazy", relation)`` for ordinary-splice joins
    across the event (seeded later, by support), ``("wobble", None)`` for
    unannotated joins beside annotated ones, or ``(None, relation)``.
    """
    relation, assignment = event.relation(left, right, unplaced)
    if unplaced == 0 and (left, right) in annotated:
        return "annotated", relation
    if relation and (relation != "splice_ambiguous_event_junction" or event.on_diagonal(*assignment, unplaced)):
        return "event", relation
    if kind == "D":
        return None, relation  # Small deletions are not splices.
    if unplaced == 0 and event.near_annotated(left, right):
        return "wobble", None
    return ("lazy", relation) if relation else (None, None)


def _seeds(observations, event, annotated, anchored, competing, key=None, deferred=None):
    """Seeds at unannotated joins and breakpoint clips, strongest adjacency placement first.

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
            use, relation = _join_class(event, annotated, left, right, j - i - 1, kind)
            if use == "annotated":
                if relation:
                    competing[left, right].add(o.fragment)
                continue
            if use == "event" and key is None:
                cores.append((i, j, relation, kind))
                continue
            regional = (use is None and relation is None and (anchored(left) or anchored(right))
                        and not event.relation(_flip(right), _flip(left), j - i - 1)[0] and kind != "D")
            if use != "lazy" and not regional:
                continue
            relation = relation or "regional_novel_junction"
            if key is None:
                if deferred is not None:
                    deferred[relation][_canonical(left, right)].add(o.identity)
            elif _canonical(left, right) == key:
                cores.append((i, j, relation, kind))
        if key is None:
            cores += [(i, j, "breakpoint_clip_partner_unplaced", "clip")
                      for i, j, _ in event.clips(o.positions, len(o.sequence))]
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
        not item[1]["near_breakpoint"], 0 if item[1]["near_breakpoint"] else rank[item[1]["relation"]],
        -len(item[1]["fragments"]), item[0][0], repr(item[0][1])))


def _contains(path, core):
    sequence, positions = path
    start = sequence.find(core[0])
    while start >= 0:
        if all(x is None or x == y for x, y in zip(core[1], positions[start:start + len(core[0])])):
            return True
        start = sequence.find(core[0], start + 1)
    return False


def _path_result(sequence, positions, voters, store, event, annotated, adjacency, clip_support, frame_evidence):
    """Serialize one path with junction, sequence and frame evidence.

    Direct support for a junction is every segment whose own observed path
    makes that join, whatever its other bases (so noisy long reads count).
    All assignments of the nominated adjacency's junction bases support its
    breakpoint junction. Counts include single records not yet built, whose
    CIGAR makes the join; built observations are cited individually.
    Returns the result and the observations it cites.
    """
    path_runs, spans_seen, placed = defaultdict(list), {}, [0]
    for q0, _, contig, low, high, strand in _runs(positions):
        path_runs[contig, strand].append((q0, low, high))
    for position in positions:
        placed.append(placed[-1] + (position is not None))

    def span(observation):
        """First and last path offsets co-observed by one read, if it agrees with the path there.

        The read must observe every placed path base between them, and place
        no other base between its first and last shared ones (no skipped or
        extra exon). Intervening unplaced query spans must also agree in
        length and sequence; otherwise it contributes no interval.
        """
        if observation.key not in spans_seen:
            shared, shared_runs = 0, []
            runs = store.observation_runs(observation)
            for own, _, contig, other_low, other_high, strand in runs:
                for q0, low, high in path_runs.get((contig, strand), ()):
                    a, b = max(low, other_low), min(high, other_high)
                    if a < b:
                        shared += b - a
                        if strand == "+":
                            start, own_start = q0 + a - low, own + a - other_low
                        else:
                            start, own_start = q0 + high - b, own + other_high - b
                        shared_runs.append((start, start + b - a, own_start, own_start + b - a))
            result = ()
            if shared:
                shared_runs.sort()
                first, last = shared_runs[0][0], shared_runs[-1][1] - 1
                own_first = min(a for _, _, a, _ in shared_runs)
                own_last = max(b for _, _, _, b in shared_runs) - 1
                own_placed = sum(max(0, min(end, own_last + 1) - max(start, own_first)) for start, end, *_ in runs)
                if (shared == placed[last + 1] - placed[first] == own_placed
                        and all(start - end == own_start - own_end >= 0
                                and sequence[end:start] == observation.sequence[own_end:own_start]
                                for (_, end, _, own_end), (start, _, own_start, _)
                                in zip(shared_runs, shared_runs[1:]))):
                    result = (first, last)
            spans_seen[observation.key] = result
        return spans_seen[observation.key]
    rows = []
    for i, j in _breaks(positions):
        left, right = positions[i], positions[j]
        relation, assignment = event.relation(left, right, j - i - 1)
        is_annotated = j == i + 1 and (left, right) in annotated
        if relation is None and is_annotated:
            continue  # Ordinary annotated splicing.
        if relation == "breakpoint_junction":
            direct, segments = adjacency, {o.identity for o, *_ in adjacency}
        else:
            direct, segments = store.joins.get((left, right), []), store.join_segments(left, right)
        kinds = {kind for *_, kind in direct}
        if segments - {o.identity for o, *_ in direct}:
            kinds |= store.gap_kinds.get((left, right), set()) | store.gap_kinds.get((_flip(right), _flip(left)), set())
        # CIGAR deletions are reported only when they also cross the event.
        if relation is None and kinds <= {"D"}:
            continue
        rows.append((i, j, relation or "regional_novel_junction", assignment, is_annotated, direct, segments, kinds))
    for i, j, assignment in event.clips(positions, len(sequence)):
        direct = clip_support["donor" if assignment[0] == 0 else "acceptor"]
        rows.append((i, j, "breakpoint_clip_partner_unplaced", assignment, False, direct,
                     {o.identity for o, *_ in direct}, {"clip"}))
    junctions, evidence = [], {}
    for i, j, relation, assignment, is_annotated, direct, segments, kinds in sorted(rows, key=lambda row: row[:2]):
        evidence.update((o.key, o) for o, *_ in direct)
        unbuilt = segments - {o.identity for o, *_ in direct}
        junction_sequences = Counter(o.sequence[a + 1:b] for o, a, b, _ in direct)
        if unbuilt:
            left, right = positions[i], positions[j]
            junction_sequences.update(inserted for identity, inserted in store.gaps.get((left, right), {}).items()
                                      if identity in unbuilt)
            junction_sequences.update(reverse_complement(inserted) for identity, inserted in
                                      store.gaps.get((_flip(right), _flip(left)), {}).items() if identity in unbuilt)
        molecules = {store.molecule(identity) for identity in segments} - {None}
        # Bases co-observed with the join in single reads; annotated splices are context only.
        spans = [] if is_annotated else [q for q in (span(o) for o, *_ in direct) if q]
        clip = relation == "breakpoint_clip_partner_unplaced"
        junctions.append(dict(
            query_interval=[i, j], left=positions[i] and list(positions[i]),
            right=positions[j] and list(positions[j]), unplaced_bases=sequence[i + 1:j],
            forward_splice_geometry=bool(positions[i] and positions[j] and _forward_splice(positions[i], positions[j])),
            kinds=sorted(kinds), annotated=is_annotated, relation=relation, breakpoint_assignment=assignment,
            direct_segments=len(segments), direct_fragments=len({s[:2] for s in segments}),
            direct_molecules=len(molecules) if molecules else None,
            direct_observations=sorted({o.key for o, *_ in direct}),
            direct_junction_sequences=None if clip else [
                list(item) for item in sorted(junction_sequences.items(), key=lambda item: (-item[1], item[0]))],
            linked_interval=[min(min(q) for q in spans), max(max(q) for q in spans) + 1] if spans else None))
    departures = {q: j["relation"] for j in junctions if not j["annotated"]
                  for q in range(j["query_interval"][0] + 1,
                                 len(sequence) if j["right"] is None else j["query_interval"][1] + 1)}
    frame = frame_evidence(sequence, positions, departures)
    relations = {j["relation"] for j in junctions if not j["annotated"]}
    status = next((s for s in LINKAGE_STATUSES if s in relations), "no_novel_junction")
    voting = [store.observations[key] for key in sorted(voters)]
    evidence.update((o.key, o) for o in voting)
    segments = {o.identity for o in voting}
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
            voting_segments=len(segments), voting_fragments=len({s[:2] for s in segments}),
            missing_quality_segments=len({o.identity for o in voting if o.missing_qualities}),
            secondary_segments=len({o.identity for o in voting if o.secondary}),
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
                       min_alternative_fraction=SV_MIN_ALTERNATIVE_FRACTION,
                       min_local_variant_fraction=SV_MIN_LOCAL_VARIANT_FRACTION, max_records=SV_MAX_RECORDS,
                       max_queries=SV_MAX_QUERIES, max_paths=SV_MAX_PATHS,
                       max_extension_segments=SV_MAX_EXTENSION_SEGMENTS, assemble=SV_ASSEMBLE,
                       breakpoint_window=SV_BREAKPOINT_WINDOW, max_breakpoint_shift=SV_MAX_BREAKPOINT_SHIFT,
                       annotated_junction_tolerance=SV_ANNOTATED_JUNCTION_TOLERANCE,
                       peptide_lengths=FUSION_PEPTIDE_LENGTHS):
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
    min_local_variant_fraction : float
        A base or small-indel alternative (placed within 20 bases of the best,
        or a variant of a seeded junction's bases) also needs this fraction of
        the best's support, so systematic long-read errors do not fork paths.
    max_records, max_queries, max_paths : int
        Resource limits; reaching one is reported in ``limitations``.
    max_extension_segments : int
        Segments built to extend one path end; those reaching furthest are
        used. Reaching it is reported as ``extension_segment_limit``.
    assemble : bool
        Extend through overlapping observations which do not span the seed.
    breakpoint_window : int
        Bases searched to either side of each nominated breakpoint.
    max_breakpoint_shift : int
        Junction bases an aligner may place past a breakpoint (homology, not
        verified against a genome) while the join still counts as the
        breakpoint junction.
    annotated_junction_tolerance : int
        Unannotated joins within this many bases of an annotated junction are
        treated as alignment wobble and seed no ordinary-splicing or regional
        path (joins at the nominated adjacency are always used).
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
                  max_extension_segments, breakpoint_window):
        if type(value) is not int or value < 1:
            raise ValueError("SV evidence thresholds and search limits must be positive integers")
    for value in (max_breakpoint_shift, annotated_junction_tolerance):
        if type(value) is not int or value < 0:
            raise ValueError("Breakpoint shift and annotated-junction tolerance must be nonnegative integers")
    for value in (min_alternative_fraction, min_local_variant_fraction):
        if isinstance(value, bool) or not 0 <= value <= 1:
            raise ValueError("Alternative-support fractions must be between 0 and 1")
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
    event = _Event(donor, acceptor, max_breakpoint_shift, annotated, annotated_junction_tolerance)
    spans = defaultdict(list)
    for r in references:
        spans[r.contig, r.strand].append((r.exons[0][0], r.exons[-1][1]))

    def anchored(position):
        """Within an annotated gene span, on its strand."""
        return any(a <= position[1] < b for a, b in spans.get((position[0], position[2]), ()))

    records, acquisition = collect_sv_records(bam, regions, references, max_records, max_queries,
                                              (donor, acceptor), breakpoint_window)
    store = _Observations(records, collector, max_extension_segments)
    now, lazy, competing = store.classify(event, annotated, anchored)
    initial = [o for identity in sorted(now) for o in store.build(identity)]
    first = _seeds(initial, event, annotated, anchored, competing, deferred=lazy)
    queue = [(relation, key, ids) for relation in ("splice_ambiguous_event_junction", "regional_novel_junction")
             for key, ids in sorted(lazy[relation].items(),
                                    key=lambda item: (-len({i[:2] for i in item[1]}), repr(item[0])))
             if len({i[:2] for i in ids}) >= min_alternative_fragments]
    thresholds = (min_alternative_fragments, min_alternative_fraction, min_local_variant_fraction)
    paths, pruned, notes, seed_rows, strongest, seen = {}, [], set(acquisition["limitations"]), [], {}, set()

    def batches():
        """(seeds, gated): event seeds first; other joins are built only if reached."""
        yield first, False
        for k, (_, key, ids) in enumerate(queue):
            if len(paths) >= max_paths:
                notes.add("path_limit")
                seed_rows.extend(dict(relation=relation, left=list(left), right=list(right),
                                      fragments=len({i[:2] for i in ids}), unexplored=True, paths=[])
                                 for relation, (left, right), ids in queue[k:])
                return
            if len(ids) > max_extension_segments:
                notes.add("seed_segment_limit")
            chosen = sorted(ids)[:max_extension_segments]
            yield _seeds([o for identity in chosen for o in store.build(identity)], event, annotated, anchored,
                         competing, key=key), True

    for batch, gated in batches():
        for (core_sequence, core_positions), seed in batch:
            if (core_sequence, core_positions) in seen:
                continue
            seen.add((core_sequence, core_positions))
            row = dict(relation=seed["relation"], kinds=sorted(seed["kinds"]),
                       left=list(core_positions[0]) if core_positions[0] else None,
                       right=list(core_positions[-1]) if core_positions[-1] else None,
                       unplaced_bases=core_sequence[1:-1], fragments=len(seed["fragments"]),
                       segments=len({o.identity for o, _ in seed["observations"]}), paths=[])
            seed_rows.append(row)
            # Joins at or across the event may rest on one fragment; ordinary
            # splices and regional joins need support.
            if gated and row["fragments"] < min_alternative_fragments:
                row["below_min_fragments"] = True
                continue
            # Seeds come strongest first. A base variant of a seeded junction is
            # pruned like a weak extension branch (only against a supported
            # best). Another placement of the nominated adjacency (a different
            # junction-base assignment or noisy near-miss) is the same event:
            # it seeds only if it is itself supported.
            variant = "adjacency" if seed["near_breakpoint"] else core_positions
            if variant in strongest:
                best, best_positions = strongest[variant]
                strong = _supported({("", None): best, ("", "this"): row["fragments"]}, *thresholds)
                if ("", "this") not in strong and (("", None) in strong or core_positions != best_positions):
                    row["minor_variant_of_seeded_junction"] = True
                    continue
            strongest.setdefault(variant, (row["fragments"], core_positions))
            if any(_contains(path, (core_sequence, core_positions)) for path in paths):
                row["represented_by_earlier_path"] = True
                continue
            if len(paths) >= max_paths:
                notes.add("path_limit")
                row["unexplored"] = True
                continue
            pool = store
            if not assemble:
                spanning = [o for o, _ in seed["observations"]]
                mirrors = {o.key: m for o in spanning for m in store.built[o.identity]
                           if m.key[:-1] == o.key[:-1] and m.key != o.key}
                pool = _ObservationIndex(spanning + [mirrors[o.key] for o in spanning])
            for sequence, positions, voters in _extend(core_sequence, core_positions, pool, min_overlap,
                                                       thresholds, max_paths - len(paths), pruned, notes, False):
                for *mirror, more in _extend(*_mirror(sequence, positions), pool, min_overlap, thresholds,
                                             max(1, max_paths - len(paths)), pruned, notes, True):
                    path = _mirror(*mirror)
                    if path not in paths and len(paths) >= max_paths:
                        notes.add("path_limit")
                        break
                    rows, contributors = paths.setdefault(path, ([], set()))
                    rows.append(row)
                    contributors.update(voters | more)
    if store.limited:
        notes.add("extension_segment_limit")
    adjacency = [entry for (left, right), entries in store.joins.items() for entry in entries
                 if event.relation(left, right, entry[2] - entry[1] - 1)[0] == "breakpoint_junction"]
    clip_support = defaultdict(list)
    for o in store.observations.values():
        for i, j, assignment in event.clips(o.positions, len(o.sequence)):
            clip_support["donor" if assignment[0] == 0 else "acceptor"].append((o, i, j, "clip"))
    results, used = [], {}
    for (sequence, positions), (rows, voters) in paths.items():
        result, evidence = _path_result(
            sequence, positions, voters, store, event, annotated, adjacency, clip_support,
            lambda *path: _frame_evidence(*path, models, min_anchor_bases, reference_peptides, lengths))
        results.append(result)
        for row in rows:
            if result["path_id"] not in row["paths"]:
                row["paths"].append(result["path_id"])
        used.update(evidence)
    rank = {status: i for i, status in enumerate(LINKAGE_STATUSES)}
    results.sort(key=lambda r: (rank.get(r["event_linkage"]["status"], len(rank)),
                                -max([j["direct_fragments"] for j in r["junctions"]], default=0),
                                r["sequence"], r["path_id"]))
    record_ids = {rid for o in used.values() for rid in o.records}
    cited = {rid: r.to_string() for identity in {o.identity for o in used.values()}
             for r in store.groups[identity] if (rid := _record_id(r)) in record_ids}
    best = next((s for s in LINKAGE_STATUSES if s in {r["event_linkage"]["status"] for r in results}), None)
    return dict(
        schema="isovar.sv_rna_candidates.v2",
        status=("event_linked_candidates" if best in LINKAGE_STATUSES[:3]
                else "splice_ambiguous_candidates" if best == "splice_ambiguous_event_junction"
                else "regional_candidates_only" if results else "no_candidate_paths"),
        event_id=event_id, reference_name=reference_name, sample_id=sample_id, source=source,
        event_provenance=event_provenance, donor=asdict(donor), acceptor=asdict(acceptor),
        reference_models=[dict(transcript_id=r.transcript_id, annotation=r.annotation, contig=r.contig,
                               strand=r.strand, sequence_sha256=sha256(r.sequence.encode()).hexdigest())
                          for r in references],
        acquisition=acquisition, limitations=sorted(notes),
        parameters=dict(min_anchor_bases=min_anchor_bases, min_overlap=min_overlap,
                        min_alternative_fragments=min_alternative_fragments,
                        min_alternative_fraction=min_alternative_fraction,
                        min_local_variant_fraction=min_local_variant_fraction, max_records=max_records,
                        max_queries=max_queries, max_paths=max_paths,
                        max_extension_segments=max_extension_segments, assemble=assemble,
                        breakpoint_window=breakpoint_window, max_breakpoint_shift=max_breakpoint_shift,
                        annotated_junction_tolerance=annotated_junction_tolerance,
                        peptide_lengths=lengths, genetic_code=1,
                        use_soft_clipped_bases=collector.use_soft_clipped_bases,
                        min_mapping_quality=collector.min_mapping_quality),
        seeds=seed_rows,
        competing_annotated_junctions=[dict(left=list(left), right=list(right), fragments=len(fragments))
                                       for (left, right), fragments in sorted(competing.items())],
        pruned_branches=pruned, paths=results,
        observations={key: dict(identity=list(o.identity), records=list(o.records), reverse_complement=o.reverse,
                                original_query_interval=list(o.query_interval),
                                missing_qualities=o.missing_qualities, secondary=o.secondary)
                      for key, o in sorted(used.items())},
        observation_counts=dict(eligible_segments=len(store.groups), built_segments=len(store.built)),
        excluded_records=dict(store.excluded), segment_path_notes=dict(store.reasons),
        original_records=cited)


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
