"""Cell barcodes and UMIs behind each small variant's allele reads.

Single-cell libraries tag reads with a cell barcode (``CB``) and a UMI (``UB``,
or an attributable Iso-Seq ``XM``). `CellUmiAlleles` counts, for each variant,
the distinct cells and cell/UMI labels behind the reference, alternate and
other alleles, and behind each protein hypothesis. Labels are resolved with the
same policy as SV reconstruction (``isovar.cell_umi_labels.v1``; see
docs/cell-umi-evidence.md). Sample-level reconstruction is unchanged: this only
counts labels on reads Isovar already used.
"""

from .cell_umi import CellUmiEvidence, summarize_cell_umi_rows
from .logging import get_logger
from .read_collector import ReadCollector
from .read_identity import segment_identity, source_read_ids
from .variant_helpers import base0_interval_for_variant

logger = get_logger(__name__)

ALLELES = ("ref", "alt", "other")

# Tags the label policy reads; keeping only these avoids holding whole records.
_LABEL_TAGS = ("CB", "UB", "XM", "CR", "XC", "UR", "RX", "OX", "PG")

# Statuses whose cell barcode cannot be trusted, even when one is present.
_UNRESOLVED_CELL = frozenset(("conflicting_tags", "conflicting_template_labels", "invalid_tags",
                              "ambiguous_read_group"))


class _Tags(object):
    """The label-relevant tags of one record, with pysam's has_tag/get_tag."""

    __slots__ = ("tags",)

    def __init__(self, record):
        self.tags = {tag: record.get_tag(tag) for tag in _LABEL_TAGS if record.has_tag(tag)}

    def has_tag(self, tag):
        return tag in self.tags

    def get_tag(self, tag):
        return self.tags[tag]


def _missing_record_row(identity):
    """A read whose record is not among the eligible records at the locus."""
    return dict(identity=list(identity), scope=None, label=None, cell_barcode=None,
                library_scope_known=False, status="metadata_unavailable")


class CellUmiAlleles(object):
    """
    Resolve cell/UMI labels for reads at variant loci in one alignment file.

    Parameters
    ----------
    alignment_file : pysam.AlignmentFile
    sample_id, source : str
        The explicit evidence scope; labels never merge across it.
    read_collector : ReadCollector, optional
        The collector the reads came from, for record eligibility. A read
        whose record this collector rejects is reported as
        ``metadata_unavailable``.
    """

    def __init__(self, alignment_file, *, sample_id, source, read_collector=None):
        self.alignment_file = alignment_file
        self.header = alignment_file.header.to_dict()
        self.scope = [sample_id, source]
        self.read_collector = read_collector or ReadCollector()

    def _labels(self, variant, templates):
        """Label resolution over eligible records of the needed templates at the variant."""
        chromosome = self.read_collector._infer_chromosome_name(
            variant.contig, set(self.alignment_file.references))
        groups = {}
        if chromosome is not None:
            start, end = base0_interval_for_variant(variant)
            # The window ReadCollector fetches, so reads ending at an insertion are included.
            for record in self.alignment_file.fetch(chromosome, max(0, start - 1), end + 1):
                identity = segment_identity(record)
                if identity[:2] in templates and self.read_collector.alignment_filter_reason(record) is None:
                    groups.setdefault(identity, []).append(_Tags(record))
        return CellUmiEvidence(groups, self.header, *self.scope)

    @staticmethod
    def _summary(labels, identities):
        rows = [labels.segment(identity) if identity in labels.groups else _missing_record_row(identity)
                for identity in identities]
        summary = summarize_cell_umi_rows(rows)
        summary.pop("segment_ids")
        trusted = [row for row in rows if row["cell_barcode"] is not None and row["status"] not in _UNRESOLVED_CELL]
        cells = {(*row["scope"].values(), row["cell_barcode"]) for row in trusted}
        complete = bool(rows) and len(trusted) == len(rows) and all(row["library_scope_known"] for row in rows)
        summary.update(segments=len(rows), observed_cells=len(cells),
                       complete_cell_count=len(cells) if complete else None)
        return summary, cells

    def evidence(self, variant, read_evidence, protein_sequences=()):
        """
        Cell/UMI label support for one variant's alleles and protein hypotheses.

        Parameters
        ----------
        variant : varcode.Variant
        read_evidence : ReadEvidence
        protein_sequences : sequence of ProteinSequence
            Proteins to summarize, in order.

        Returns
        -------
        dict
            ``policy``; ``alleles`` with a summary for ``ref``, ``alt`` and
            ``other``; ``cells_with_ref_and_alt``, cells with a trusted barcode
            on reads of both alleles; and one summary per protein in
            ``protein_hypotheses``. Each summary gives ``segments``,
            ``observed_labels`` and ``observed_cells``, ``complete_label_count``
            and ``complete_cell_count`` (null unless every read is resolved in
            a known library), unresolved and unknown-library segments and status
            counts. Cells come from the barcode alone, with or without a UMI.
            Counts are not molecules or cell prevalence, and must not be added
            across alleles or proteins that share cells.
        """
        by_allele = {allele: sorted({key for read in getattr(read_evidence, allele + "_reads")
                                     for key in source_read_ids(read)}) for allele in ALLELES}
        by_protein = [sorted({key for read in protein.supporting_reads for key in source_read_ids(read)})
                      for protein in protein_sequences]
        templates = {identity[:2] for identities in [*by_allele.values(), *by_protein] for identity in identities}
        labels = self._labels(variant, templates)
        missing = sum(identity not in labels.groups for identities in by_allele.values() for identity in identities)
        if missing:
            logger.warning(
                "%d read(s) at %s had no eligible record for cell/UMI labels; pass the ReadCollector "
                "the reads were collected with", missing, variant.short_description)
        alleles, cells = {}, {}
        for allele in ALLELES:
            alleles[allele], cells[allele] = self._summary(labels, by_allele[allele])
        return dict(
            policy=CellUmiEvidence.policy, alleles=alleles,
            cells_with_ref_and_alt=len(cells["ref"] & cells["alt"]),
            protein_hypotheses=[self._summary(labels, identities)[0] for identities in by_protein])


def warn_if_unlabelled(evidence):
    """Log a warning when no read carried a cell barcode: likely not single-cell data."""
    summaries = [summary for item in evidence for summary in item["alleles"].values()]
    if any(s["segments"] for s in summaries) and not any(s["observed_cells"] for s in summaries):
        logger.warning("No read carried a usable CB cell barcode; cell counts are all zero")


def cell_umi_allele_evidence(isovar_results, alignment_file, *, sample_id, source, read_collector=None):
    """
    Cell/UMI label evidence for each result's alleles and protein hypotheses.

    Parameters
    ----------
    isovar_results : iterable of IsovarResult
    alignment_file : pysam.AlignmentFile
        The alignments the results were collected from.
    sample_id, source : str
        Evidence scope, as in `export_protein_hypotheses`.
    read_collector : ReadCollector, optional
        The collector used for the results.

    Returns
    -------
    list of dict
        One `CellUmiAlleles.evidence` per result, in order.
    """
    labels = CellUmiAlleles(alignment_file, sample_id=sample_id, source=source, read_collector=read_collector)
    evidence = [labels.evidence(r.variant, r.read_evidence, r.sorted_protein_sequences) for r in isovar_results]
    warn_if_unlabelled(evidence)
    return evidence
