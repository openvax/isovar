"""Frozen candidate support with bounded transient observation geometry."""

from collections import Counter, defaultdict

from .sv_rna import _PathSpans, _breaks
from .sv_rna_orfs import _witnesses


# Internal working-set bound, independent of discovery and scientific settings.
# Templates are indivisible; a single oversized template is consumed on its own.
_SUPPORT_BATCH_BASES = 250_000


def _observation_row(o):
    return dict(identity=list(o.identity), records=list(o.records),
                reverse_complement=o.reverse_complement,
                original_query_interval=list(o.query_interval),
                missing_qualities=o.missing_qualities, secondary=o.secondary)


def _merge_interval(previous, interval):
    if not interval:
        return previous
    return (min(previous[0], interval[0]), max(previous[1], interval[1])) if previous else interval


class _SupportStream:
    """Transpose batches into output-sized junction and interval evidence.

    No observation, sequence, quality vector or coordinate tuple from a batch
    survives ``consume``. The retained rows are original observation metadata,
    direct-junction sequence counts, and exact witnesses required by the output.
    """

    def __init__(self, paths, positions, event, annotated):
        self.paths, self.positions, self.event = paths, positions, event
        self.queries, self.path_queries = defaultdict(list), []
        self.pairs = defaultdict(set)
        self.cache_keys, self.rows, self.counts, self.kinds = {}, {}, defaultdict(Counter), defaultdict(set)
        self.direct = defaultdict(dict)
        self.spans = [_PathSpans(p['sequence'], pos) for p, pos in zip(paths, positions)]
        self.linked = [dict() for _ in paths]
        self.witnesses = [defaultdict(list) for _ in paths]
        for path_index, (path, pos) in enumerate(zip(paths, positions)):
            queries = []
            for i, j in _breaks(pos):
                left, right = pos[i], pos[j]
                relation, _ = event.relation(left, right, j - i - 1)
                if relation is None and j == i + 1 and (left, right) in annotated:
                    continue
                cache_key = (left, right, relation)
                key = 'event' if relation == 'breakpoint_junction' else (left, right)
                self.cache_keys[cache_key] = key
                if key != 'event':
                    self.pairs[left, right].add(key)
                queries.append((i, j, key))
            for i, j, assignment in event.clips(pos, len(path['sequence'])):
                key = ('clip', 'donor' if assignment[0] == 0 else 'acceptor')
                self.cache_keys[key] = key
                queries.append((i, j, key))
            self.path_queries.append([
                dict(query_interval=[i + 1, j], left=pos[i] and list(pos[i]), right=pos[j] and list(pos[j]), key=key)
                for i, j, key in sorted(queries, key=lambda row: row[:2])])
            for key in {key for _, _, key in queries}:
                self.queries[key].append(path_index)

    def consume(self, observations):
        """Record a batch's contributions, releasing all per-base objects."""
        direct = defaultdict(list)
        for o in observations:
            for i, j, kind in o.breaks:
                left, right = o.positions[i], o.positions[j]
                keys = self.pairs.get((left, right), ())
                for key in keys:
                    direct[key].append((o, i, j, kind))
                if 'event' in self.queries and self.event.relation(left, right, j - i - 1)[0] == 'breakpoint_junction':
                    direct['event'].append((o, i, j, kind))
            for i, j, assignment in self.event.clips(o.positions, len(o.sequence)):
                key = ('clip', 'donor' if assignment[0] == 0 else 'acceptor')
                if key in self.queries:
                    direct[key].append((o, i, j, 'clip'))
        by_key, affected, groups = {}, set(), {}
        for key, entries in direct.items():
            evidence = {o.key: o for o, *_ in entries}
            by_key[key] = evidence
            for observation_key, o in evidence.items():
                row = self.rows.setdefault(observation_key, _observation_row(o))
                self.direct[key][observation_key] = row
            self.kinds[key].update(kind for *_, kind in entries)
            self.counts[key].update(o.sequence[i + 1:j] for o, i, j, _ in entries)
            equivalent = {o.comparison_key: o for o in evidence.values()}
            for path_index in self.queries[key]:
                affected.add(path_index)
                span = self.spans[path_index]
                interval = self.linked[path_index].get(key, ())
                for o in equivalent.values():
                    interval = _merge_interval(interval, span(o))
                self.linked[path_index][key] = interval
        for path_index in affected:
            path, pos = self.paths[path_index], self.positions[path_index]
            junctions = [dict(q, direct_observations=sorted(by_key[q['key']]))
                         for q in self.path_queries[path_index] if q['key'] in by_key]
            evidence = {key: o for q in junctions for key, o in by_key[q['key']].items()}
            intervals = {(0, len(pos))} | {tuple(c['query_interval']) for c in path['exploratory_orfs']['candidates']}
            for start, end in intervals:
                self.witnesses[path_index][start, end].extend(
                    _witnesses(path['sequence'], pos, start, end, junctions, evidence, group_cache=groups))
            self.spans[path_index].cache.clear()

    def summaries(self, cell_umi, lineage):
        """Resolve original identities once, after all batches are accounted for."""
        summaries, cache = {}, {}
        for cache_key, key in self.cache_keys.items():
            if len(cache_key) == 3 and cache_key[2] is None and self.kinds[key] <= {'D'}:
                cache[cache_key] = None
                continue
            if key not in summaries:
                evidence = self.direct[key]
                identities = {tuple(row['identity']) for row in evidence.values()}
                summaries[key] = dict(
                    evidence=evidence, spanning=[], kinds=sorted(self.kinds[key]),
                    direct_support=cell_umi.support(identities), direct_read_lineage=lineage.support(identities),
                    direct_observations=sorted(evidence),
                    direct_junction_sequences=[list(item) for item in sorted(
                        self.counts[key].items(), key=lambda item: (-item[1], item[0]))])
            cache[cache_key] = summaries[key]
        return cache, [{id(summaries[key]): [interval] for key, interval in linked.items()
                        if key in summaries and interval} for linked in self.linked]

    def interval_witnesses(self, path_index, start, end, junctions):
        """Keep only links to the final reported/crossed junctions, in original order."""
        selected = {tuple(j['query_interval']) for j in junctions}
        rows = []
        for witness in self.witnesses[path_index][start, end]:
            links = [link for link in witness['junction_links'] if tuple(link['junction_query_interval']) in selected]
            if links:
                rows.append(dict(witness, junction_links=links))
        return sorted(rows, key=lambda row: (row['observation'], row['observation_interval'][0]))
