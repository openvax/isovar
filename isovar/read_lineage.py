"""Producer-aware signal ancestry, separate from alignment and molecule IDs."""

from collections import Counter
import shlex

from .read_metadata import unique_header_entries


def _is_dorado_basecaller(program):
    # Dorado also has aligner/trim tools, which cannot establish tag origin.
    if program.get("PN", "").lower() != "dorado":
        return False
    try:
        command = shlex.split(program.get("CL", ""))
    except ValueError:
        return False
    return (len(command) >= 2 and command[0].rsplit("/", 1)[-1].lower() == "dorado"
            and command[1] in ("basecaller", "duplex"))


class ReadLineage:
    """Interpret Dorado ancestry lazily for one input's eligible segment groups.

    Parameters
    ----------
    groups : mapping
        ``(RG, QNAME, segment bits)`` to original pysam records, including all
        retained alternative placements. No sequences are combined here.
    header : dict
        SAM header (``AlignmentHeader.to_dict()``).

    Notes
    -----
    Signal groups are not biological molecules. Unknown producers, ambiguous
    metadata and duplex relationships remain unresolved. Read groups never
    merge, even if their library labels match. Callers must keep inputs separate.
    """

    def __init__(self, groups, header):
        self.groups = groups
        self.read_groups = unique_header_entries(header.get("RG", []))
        self.programs = unique_header_entries(header.get("PG", []))
        self.program_is_dorado = {key: self._dorado_chain(key) for key in self.programs}
        # Multiple unrelated chains (e.g. a merged BAM) cannot identify a record's
        # producer without an explicit PG pointer. All roots must be Dorado.
        parents = {p.get("PP") for p in self.programs.values()}
        leaves = set(self.programs) - parents
        self.header_is_dorado = (bool(leaves) and len(self.programs) == len(header.get("PG", []))
                                 and all(self.program_is_dorado.values()))
        self.cache = {}

    def _dorado_chain(self, key):
        seen, found = set(), False
        while key is not None:
            if key in seen or key not in self.programs:
                return False
            seen.add(key)
            program = self.programs[key]
            found |= _is_dorado_basecaller(program)
            key = program.get("PP")
        return found

    def _producer(self, read, group):
        if group.get("PL", "").upper() != "ONT":
            return False
        program = read.get_tag("PG") if read.has_tag("PG") else group.get("PG")
        if program is None:
            return self.header_is_dorado
        return isinstance(program, str) and self.program_is_dorado.get(program, False)

    def _metadata(self, identity):
        group = self.read_groups.get(identity[0], {})
        row = dict(identity=list(identity), producer=None, status="unknown_producer",
                   parent_read_id=None, duplex_status=None, split_signal_start=None, signal_group=None)
        records = self.groups[identity]
        if not all(self._producer(read, group) for read in records):
            return row
        row["producer"] = "dorado"
        values = [tuple(read.get_tag(tag) if read.has_tag(tag) else None for tag in ("pi", "dx", "sp"))
                  for read in records]
        parent, duplex, start = values[0]
        if any(any(type(a) is not type(b) or a != b for a, b in zip(value, values[0]))
               for value in values[1:]):
            row["status"] = "conflicting_tags"
        elif (parent is not None and (not isinstance(parent, str) or not parent or parent == identity[1])
              or duplex is not None and (type(duplex) is not int or duplex not in (-1, 0, 1))
              or start is not None and (type(start) is not int or start < 0 or parent is None)):
            row["status"] = "invalid_tags"
        else:
            row.update(parent_read_id=parent, duplex_status=duplex, split_signal_start=start)
            row["status"] = ("unsupported_segment" if identity[2] or any(r.is_paired for r in records)
                             else "duplex_parent" if duplex == -1
                             else "duplex_consensus" if duplex == 1
                             else "split_read" if parent is not None else "simplex_read")
        return row

    def segment(self, identity):
        """Return JSON-ready ancestry evidence, checking visible parent chains."""
        pending, seen, current = [], set(), identity
        while current not in self.cache:
            if current in seen:
                break
            seen.add(current)
            row = self._metadata(current)
            pending.append((current, row))
            if row["status"] not in ("split_read", "simplex_read"):
                break
            parent = row["parent_read_id"]
            parent_identity = (current[0], parent, 0)
            if parent is None or parent_identity not in self.groups:
                row["signal_group"] = [current[0], parent or current[1]]
                break
            current = parent_identity
        signal_group = self.cache[current]["signal_group"] if current in self.cache else pending[-1][1]["signal_group"]
        for key, row in reversed(pending):
            if row["status"] in ("split_read", "simplex_read"):
                row["signal_group"] = signal_group
                if signal_group is None:
                    row["status"] = "unresolved_parent"
            self.cache[key] = row
        return self.cache[identity]

    def support(self, identities):
        """Count resolved signal groups only among the supplied witnesses.

        Every supplied segment is listed for audit, including unresolved ones.
        A known subset count must not be treated as the total support count.
        """
        rows = [self.segment(identity) for identity in sorted(set(identities))]
        return summarize_lineage_rows(rows)

    def evidence(self):
        """Evidence for requested segments and any consulted visible parents."""
        return [self.cache[key] for key in sorted(self.cache)]


def summarize_lineage_rows(rows):
    """Summarize unique serialized segment rows using the reconstruction policy."""
    rows = list(rows)
    groups = {tuple(row["signal_group"]) for row in rows if row["signal_group"] is not None}
    unresolved = sum(row["signal_group"] is None for row in rows)
    return dict(segment_ids=[row["identity"] for row in rows],
                resolved_signal_groups=len(groups), unresolved_segments=unresolved,
                all_segments_resolved=not unresolved,
                status_counts=dict(sorted(Counter(row["status"] for row in rows).items())))
