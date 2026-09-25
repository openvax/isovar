"""Input-scoped cell/UMI labels, separate from independent molecule counts."""

from collections import Counter, defaultdict
import shlex

from .read_metadata import ProgramHistory, unique_header_entries


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
        self.history = ProgramHistory(header)
        self.program_cache, self.cache = {}, {}
        self.header_xm_producer = False
        self.templates = defaultdict(list)
        for identity in groups:
            self.templates[identity[:2]].append(identity)

    def _xm_program(self, key):
        if key not in self.program_cache:
            self.program_cache[key] = self._xm_state(self.history.chain(key))
        return self.program_cache[key]

    @staticmethod
    def _xm_state(chain):
        """Iso-Seq XM semantics from a program chain, newest first."""
        state = None
        for program in chain or ():
            if "bismark" in program.get("PN", program["ID"]).lower():
                return None
            step = _isoseq_step(program)
            # The latest tag/correct step determines whether XM is raw.
            if state in (None, "isoseq_correction_unknown"):
                if step == "correct":
                    state = "isoseq_corrected"
                elif step == "tag":
                    state = "isoseq_raw"
                elif step in ("dedup", "groupdedup"):
                    state = "isoseq_correction_unknown"
        return state

    def _xm_producer(self, read, group):
        program = self.history.pointer(read, group)
        if program is not None:
            return self._xm_program(program) if isinstance(program, str) else None
        # Without a record/RG pointer, every header lineage must agree. A
        # disconnected aligner/merged input cannot identify this tag's producer.
        if self.header_xm_producer is False:
            states = {self._xm_program(key) for key in self.history.leaves}
            self.header_xm_producer = (
                next(iter(states)) if self.history.unique and len(states) == 1
                and self.history.leaves_cover_all() else None)
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
        """The RNA support record for these reads, with their cell/UMI counts.

        Reads without an eligible record here count as unlabelled
        (``metadata_unavailable``). See `isovar.rna_evidence.rna_support`.
        """
        from .rna_evidence import rna_support

        return rna_support(identities, self.row)

    def row(self, identity):
        """The label row for one read, or a metadata_unavailable row without a record."""
        return self.segment(identity) if identity in self.groups else missing_label_row(identity)

    def shared_cells(self, first, second):
        """How many trusted cells have reads among both sets of read identities."""
        return len(trusted_cells(map(self.row, first)) & trusted_cells(map(self.row, second)))

    def evidence(self):
        """JSON-ready policy and per-read label rows for consulted reads and mates."""
        return dict(policy=self.policy, reads=[self.cache[key] for key in sorted(self.cache)])


# Statuses whose cell barcode cannot be trusted, even when one is present.
UNTRUSTED_CELL_STATUSES = frozenset(("conflicting_tags", "conflicting_template_labels", "invalid_tags",
                                     "ambiguous_read_group", "metadata_unavailable"))


def missing_label_row(identity):
    """A read whose record was not available to resolve its label."""
    return dict(identity=list(identity), scope=None, label=None, cell_barcode=None,
                library_scope_known=False, status="metadata_unavailable")


def trusted_cells(rows):
    """Cells, keyed by label scope and barcode, from rows with a trusted barcode."""
    return {(*row["scope"].values(), row["cell_barcode"]) for row in rows
            if row["cell_barcode"] is not None and row["status"] not in UNTRUSTED_CELL_STATUSES}


def cell_umi_counts(rows):
    """
    Cell/UMI fields of an RNA support record from label rows, one per read.

    ``umis`` counts distinct cell barcode and UMI pairs, and ``cells`` distinct
    cell barcodes, each within its declared library. A barcode without a UMI
    still identifies its cell. ``umis_complete`` and ``cells_complete`` say
    whether there are reads and every one contributed, with a known library,
    so that the count is exact; without reads they are False. Neither count
    is a molecule or prevalence estimate.
    """
    rows = list(rows)
    umis = {tuple(row["label"]) for row in rows if row["label"] is not None}
    cells = trusted_cells(rows)
    untrusted = sum(row["cell_barcode"] is None or row["status"] in UNTRUSTED_CELL_STATUSES for row in rows)
    unlabeled = sum(row["label"] is None for row in rows)
    unknown_library = sum(not row["library_scope_known"] for row in rows)
    return dict(umis=len(umis), cells=len(cells),
                umis_complete=bool(rows) and not unlabeled and not unknown_library,
                cells_complete=bool(rows) and not untrusted and not unknown_library,
                unlabeled_reads=unlabeled, unknown_library_reads=unknown_library,
                label_statuses=dict(sorted(Counter(row["status"] for row in rows).items())))
