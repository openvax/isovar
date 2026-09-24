"""Input-scoped cell/UMI labels, separate from independent molecule counts."""

from collections import Counter, defaultdict
import shlex

from .read_metadata import unique_header_entries


def _isoseq_step(program):
    """Recognize named Iso-Seq steps, including its PN-less emitted PG format."""
    try:
        command = shlex.split(program.get("CL", ""))
    except ValueError:
        return None
    if not command:
        return None
    executable = command[0].rsplit("/", 1)[-1]
    if program.get("PN") in ("isoseq", "isoseq3") and executable in ("isoseq", "isoseq3"):
        return command[1] if len(command) > 1 else None
    if not program.get("PN") and program.get("ID") == "isoseq " + executable:
        return executable
    return None


class CellUmiEvidence:
    """Resolve reported CB/UMI labels for one input's eligible alignments.

    Parameters
    ----------
    groups : mapping
        ``(RG, QNAME, segment bits)`` to original records, including retained
        alternative placements. No records or sequences are changed.
    header : dict
        SAM header. Matching SM/LB identify a shared declared library across RGs.
    sample_id, source : str
        Explicit input sample and source namespaces; not inferred from paths.

    Notes
    -----
    CB/UB are used as supplied identifiers, without raw-tag fallback, barcode
    correction or suffix removal. XM is accepted only with attributable
    Iso-Seq correction provenance. Matching labels do not establish independent
    molecules, expression, cell prevalence or phase. Missing LB remains unknown.
    """

    policy = "isovar.cell_umi_labels.v1"

    def __init__(self, groups, header, sample_id, source):
        self.groups, self.sample_id, self.source = groups, sample_id, source
        self.read_groups = unique_header_entries(header.get("RG", []))
        self.ambiguous_groups = {r.get("ID") for r in header.get("RG", [])} - self.read_groups.keys()
        self.programs = unique_header_entries(header.get("PG", []))
        self.programs_unique = len(self.programs) == len(header.get("PG", []))
        self.program_cache, self.cache = {}, {}
        self.header_xm_producer = False
        self.templates = defaultdict(list)
        for identity in groups:
            self.templates[identity[:2]].append(identity)

    def _xm_program(self, key):
        if key in self.program_cache:
            return self.program_cache[key]
        current, seen, state = key, set(), None
        while current is not None:
            if current in seen or current not in self.programs:
                state = None
                break
            seen.add(current)
            program = self.programs[current]
            if "bismark" in program.get("PN", program["ID"]).lower():
                state = None
                break
            step = _isoseq_step(program)
            # The latest tag/correct step determines whether XM is raw.
            if state in (None, "isoseq_correction_unknown"):
                if step == "correct":
                    state = "isoseq_corrected"
                elif step == "tag":
                    state = "isoseq_raw"
                elif step in ("dedup", "groupdedup"):
                    state = "isoseq_correction_unknown"
            current = program.get("PP")
        self.program_cache[key] = state
        return state

    def _xm_producer(self, read, group):
        program = read.get_tag("PG") if read.has_tag("PG") else group.get("PG")
        if program is not None:
            return self._xm_program(program) if isinstance(program, str) else None
        # Without a record/RG pointer, every header lineage must agree. A
        # disconnected aligner/merged input cannot identify this tag's producer.
        if self.header_xm_producer is False:
            parents = {p.get("PP") for p in self.programs.values()}
            leaves = self.programs.keys() - parents
            states = {self._xm_program(key) for key in leaves}
            covered, pending = set(), list(leaves)
            while pending:
                key = pending.pop()
                if key in self.programs and key not in covered:
                    covered.add(key)
                    pending.append(self.programs[key].get("PP"))
            self.header_xm_producer = (next(iter(states)) if self.programs_unique and len(states) == 1
                                      and covered == self.programs.keys() else None)
        return self.header_xm_producer

    def _metadata(self, identity):
        group = self.read_groups.get(identity[0], {})
        library, sample = group.get("LB") or None, group.get("SM") or None
        scope_known = bool(library) and identity[0] not in self.ambiguous_groups
        scope = dict(source=self.source, sample_id=self.sample_id, header_sample=sample, library=library,
                     read_group=None if scope_known else identity[0],
                     basis="sample_library" if scope_known else "read_group" if identity[0] else "input")
        row = dict(identity=list(identity), scope=scope, library_scope_known=scope_known,
                   status="missing_label", cell_barcode=None, umi=None, umi_tags=[],
                   xm_semantics=[], records_without_complete_label=0, label=None)
        cells, umis, pairs, tags, xm_states = set(), set(), set(), set(), set()
        invalid, raw_cell, raw_umi = False, False, False
        for read in self.groups[identity]:
            values = {tag: read.get_tag(tag) if read.has_tag(tag) else None for tag in ("CB", "UB", "XM")}
            raw_cell |= read.has_tag("CR") or read.has_tag("XC")
            raw_umi |= any(read.has_tag(tag) for tag in ("UR", "RX", "OX"))
            state = self._xm_producer(read, group) if values["XM"] is not None else None
            if values["XM"] is not None:
                xm_states.add(state or "unknown_producer")
            selected = ["UB"] + (["XM"] if state == "isoseq_corrected" else [])
            relevant = [values[tag] for tag in ("CB", *selected) if values[tag] is not None]
            if any(not isinstance(v, str) or not v or v.isspace() for v in relevant):
                invalid = True
                row["records_without_complete_label"] += 1
                continue
            cell = values["CB"]
            if cell is not None:
                cells.add(cell)
            record_umis = {values[tag] for tag in selected if values[tag] is not None}
            tags.update(tag for tag in selected if values[tag] is not None)
            umis.update(record_umis)
            if cell is not None and len(record_umis) == 1:
                pairs.add((cell, next(iter(record_umis))))
            else:
                row["records_without_complete_label"] += 1
        row.update(umi_tags=sorted(tags), xm_semantics=sorted(xm_states),
                   cell_barcode=next(iter(cells)) if len(cells) == 1 else None,
                   umi=next(iter(umis)) if len(umis) == 1 else None)
        if identity[0] in self.ambiguous_groups:
            row["status"] = "ambiguous_read_group"
        elif invalid:
            row["status"] = "invalid_tags"
        elif len(cells) > 1 or len(umis) > 1:
            row["status"] = "conflicting_tags"
        elif "XM" in tags and xm_states != {"isoseq_corrected"}:
            row["status"] = "unresolved_xm"
        elif len(pairs) == 1:
            row["status"] = "resolved_label"
            row["label"] = [*scope.values(), *next(iter(pairs))]
        elif cells and umis:
            row["status"] = "partial_label"
        elif not cells:
            row["status"] = "raw_cell_barcode_only" if raw_cell else "missing_cell_barcode"
        elif xm_states:
            row["status"] = "unresolved_xm"
        else:
            row["status"] = "raw_umi_only" if raw_umi else "missing_umi"
        return row

    def segment(self, identity):
        """Resolve a segment and check conflicting labels on its visible mates."""
        if identity not in self.cache:
            rows = [self._metadata(key) for key in self.templates[identity[:2]]]
            # One template cannot acquire two cell/UMI labels merely because
            # its mates or alternative placements disagree.
            conflict = (len(rows) > 1 and (any(r["status"] == "conflicting_tags" for r in rows)
                        or any(len({r[field] for r in rows if r[field] is not None}) > 1
                               for field in ("cell_barcode", "umi"))))
            for row in rows:
                if conflict:
                    row.update(status="conflicting_template_labels", label=None)
                self.cache[tuple(row["identity"])] = row
        return self.cache[identity]

    def support(self, identities):
        """Count reported labels only among the supplied segment witnesses.

        ``observed_labels`` is a subset count. ``complete_label_count`` is null
        unless every witness has a label and known library scope. Neither is
        an independent molecule count; no clustering or collision model is used.
        """
        rows = [self.segment(key) for key in sorted(set(identities))]
        return summarize_cell_umi_rows(rows)

    def evidence(self):
        """JSON-ready policy and evidence for consulted segments and mates."""
        return dict(policy=self.policy, segments=[self.cache[key] for key in sorted(self.cache)])


def summarize_cell_umi_rows(rows):
    """Summarize unique serialized segment rows using the reconstruction policy."""
    rows = list(rows)
    labels = {tuple(row["label"]) for row in rows if row["label"] is not None}
    unresolved = sum(row["label"] is None for row in rows)
    unknown_scope = sum(not row["library_scope_known"] for row in rows)
    return dict(unit="cell_umi_label", segment_ids=[row["identity"] for row in rows],
                observed_labels=len(labels), unresolved_segments=unresolved,
                unknown_library_segments=unknown_scope, all_segments_labeled=not unresolved,
                complete_label_count=len(labels) if rows and not unresolved and not unknown_scope else None,
                independent_molecules=None,
                status_counts=dict(sorted(Counter(row["status"] for row in rows).items())))
