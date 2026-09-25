"""Cell barcodes and UMIs behind each small variant's allele reads.

Single-cell libraries tag reads with a cell barcode (``CB``) and a UMI (``UB``,
or an attributable Iso-Seq ``XM``). `CellUmiAlleles` reports, for each variant,
how many distinct labels and cells carry the reference, alternate and other
alleles, and the reads behind each protein hypothesis. Labels are resolved with
the same policy as SV reconstruction (``isovar.cell_umi_labels.v1``; see
docs/cell-umi-evidence.md). Sample-level reconstruction is unchanged: this only
counts labels among reads Isovar already used.
"""

from .cell_umi import CellUmiEvidence, summarize_cell_umi_rows
from .read_collector import ReadCollector
from .read_identity import source_read_ids, segment_identity
from .rna_evidence import segment_support
from .variant_helpers import base0_interval_for_variant

ALLELES = ("ref", "alt", "other")


class CellUmiAlleles(object):
    """
    Resolve cell/UMI labels for reads at variant loci in one alignment file.

    Parameters
    ----------
    alignment_file : pysam.AlignmentFile
    sample_id, source : str
        The explicit evidence scope; labels never merge across it.
    read_collector : ReadCollector, optional
        Record eligibility, as used to collect the reads.
    """

    def __init__(self, alignment_file, *, sample_id, source, read_collector=None):
        self.alignment_file = alignment_file
        self.header = alignment_file.header.to_dict()
        self.scope = [sample_id, source]
        self.read_collector = read_collector or ReadCollector()

    def _labels(self, variant):
        """Label resolution over the eligible records overlapping the variant."""
        chromosome = self.read_collector._infer_chromosome_name(
            variant.contig, set(self.alignment_file.references))
        groups = {}
        if chromosome is not None:
            start, end = base0_interval_for_variant(variant)
            for record in self.alignment_file.fetch(chromosome, start, max(end, start + 1)):
                if self.read_collector.alignment_filter_reason(record) is None:
                    groups.setdefault(segment_identity(record), []).append(record)
        return CellUmiEvidence(groups, self.header, *self.scope)

    def _summary(self, labels, reads):
        identities = sorted({key for read in reads for key in source_read_ids(read)})
        rows = [labels.segment(identity) for identity in identities if identity in labels.groups]
        summary = summarize_cell_umi_rows(rows)
        summary.pop("segment_ids")
        cells = {tuple(row["label"][:-1]) for row in rows if row["label"] is not None}
        summary.update(
            segments=len(identities),
            # Reads from records no longer at the locus cannot be labelled here.
            segments_without_records=len(identities) - len(rows),
            observed_cells=len(cells),
            complete_cell_count=len(cells) if summary["complete_label_count"] is not None
            and len(rows) == len(identities) else None,
            evidence_set_id=segment_support(self.scope, identities)["evidence_set_id"] if identities else None)
        if len(rows) != len(identities):
            summary["complete_label_count"] = None
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
            ``other``; ``cells_with_ref_and_alt``, cells with a resolved label
            on reads of both alleles; and one summary per protein in
            ``protein_hypotheses``. Each summary gives ``observed_labels`` and
            ``observed_cells`` (lower bounds), ``complete_label_count`` and
            ``complete_cell_count`` (null unless every segment has a resolved
            label with known library), unresolved and unknown-library segments
            and status counts. Labels are not molecules or cell prevalence.
        """
        labels = self._labels(variant)
        alleles, cells = {}, {}
        for allele in ALLELES:
            alleles[allele], cells[allele] = self._summary(labels, getattr(read_evidence, allele + "_reads"))
        return dict(
            policy=CellUmiEvidence.policy, alleles=alleles,
            cells_with_ref_and_alt=len(cells["ref"] & cells["alt"]),
            protein_hypotheses=[self._summary(labels, protein.supporting_reads)[0]
                                for protein in protein_sequences])


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
    return [labels.evidence(r.variant, r.read_evidence, r.sorted_protein_sequences) for r in isovar_results]
