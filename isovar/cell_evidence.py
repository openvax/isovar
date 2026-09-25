"""Cell barcodes and UMIs behind each small variant's allele reads.

Single-cell libraries tag reads with a cell barcode (``CB``) and a UMI (``UB``,
or an attributable Iso-Seq ``XM``). `CellUmiAlleles` gives, for each variant,
the RNA support record (reads, fragments, UMIs and cells; see
`isovar.rna_evidence.rna_support`) of the reference, alternate and other
alleles and of each protein hypothesis. Labels are resolved with the same
policy as SV reconstruction (``isovar.cell_umi_labels.v1``; see
docs/cell-umi-evidence.md). Sample-level reconstruction is unchanged: this only
counts labels on reads Isovar already used.
"""

from .cell_umi import CellUmiEvidence
from .logging import get_logger
from .read_collector import ReadCollector
from .read_identity import segment_identity, source_read_ids
from .variant_helpers import base0_interval_for_variant

logger = get_logger(__name__)

ALLELES = ("ref", "alt", "other")

# Tags the label policy reads; keeping only these avoids holding whole records.
_LABEL_TAGS = ("CB", "UB", "XM", "CR", "XC", "UR", "RX", "OX", "PG")


class _Tags(object):
    """The label-relevant tags of one record, with pysam's has_tag/get_tag."""

    __slots__ = ("tags",)

    def __init__(self, record):
        self.tags = {tag: record.get_tag(tag) for tag in _LABEL_TAGS if record.has_tag(tag)}

    def has_tag(self, tag):
        return tag in self.tags

    def get_tag(self, tag):
        return self.tags[tag]


def _identities(reads):
    return sorted({key for read in reads for key in source_read_ids(read)})


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
        whose record this collector rejects is unlabelled
        (``metadata_unavailable``).
    """

    def __init__(self, alignment_file, *, sample_id, source, read_collector=None):
        self.alignment_file = alignment_file
        self.header = alignment_file.header.to_dict()
        self.scope = [sample_id, source]
        self.read_collector = read_collector or ReadCollector()

    def labels(self, variant, identities):
        """
        The label resolver (`isovar.cell_umi.CellUmiEvidence`) for these reads
        at a variant. Its ``support(identities)`` gives RNA support records.
        """
        start, end = base0_interval_for_variant(variant)
        return self.labels_at(variant.contig, start, end, identities)

    def labels_at(self, contig, start, end, identities):
        """
        The label resolver for these reads over the 0-based interval [start, end).

        Logs a warning when some of the reads have no eligible record there,
        which usually means they were collected with a different ReadCollector.
        """
        identities = {tuple(identity) for identity in identities}
        templates = {identity[:2] for identity in identities}
        chromosome = self.read_collector._infer_chromosome_name(contig, set(self.alignment_file.references))
        groups = {}
        if chromosome is not None:
            # The window ReadCollector fetches, so reads ending at an insertion are included.
            for record in self.alignment_file.fetch(chromosome, max(0, start - 1), end + 1):
                identity = segment_identity(record)
                if identity[:2] in templates and self.read_collector.alignment_filter_reason(record) is None:
                    groups.setdefault(identity, []).append(_Tags(record))
        missing = len(identities - groups.keys())
        if missing:
            logger.warning(
                "%d read(s) at %s:%d-%d had no eligible record for cell/UMI labels; pass the ReadCollector "
                "the reads were collected with", missing, contig, start + 1, end)
        return CellUmiEvidence(groups, self.header, *self.scope)

    def evidence(self, variant, read_evidence, protein_sequences=()):
        """
        RNA support records, with cell/UMI counts, for one variant's alleles and proteins.

        Parameters
        ----------
        variant : varcode.Variant
        read_evidence : ReadEvidence
        protein_sequences : sequence of ProteinSequence
            Proteins to summarize, in order.

        Returns
        -------
        dict
            ``policy``; ``alleles``, the record of each of ``ref``, ``alt`` and
            ``other``; ``cells_with_ref_and_alt``, cells with reads of both
            alleles; and ``protein_hypotheses``, one record per protein. Counts
            are not molecules or cell prevalence, and must not be added across
            alleles or proteins that share cells.
        """
        by_allele = {allele: _identities(getattr(read_evidence, allele + "_reads")) for allele in ALLELES}
        by_protein = [_identities(protein.supporting_reads) for protein in protein_sequences]
        labels = self.labels(variant, [i for identities in [*by_allele.values(), *by_protein] for i in identities])
        return dict(
            policy=CellUmiEvidence.policy,
            alleles={allele: labels.support(by_allele[allele]) for allele in ALLELES},
            cells_with_ref_and_alt=labels.shared_cells(by_allele["ref"], by_allele["alt"]),
            protein_hypotheses=[labels.support(identities) for identities in by_protein])


def warn_if_unlabelled(supports):
    """
    Log a warning when reads with a record carried no usable cell barcode:
    likely not single-cell data. ``supports`` are RNA support records.
    """
    supports = list(supports)
    with_record = sum(s["reads"] - s["label_statuses"].get("metadata_unavailable", 0) for s in supports)
    if with_record and not any(s["cells"] for s in supports):
        logger.warning("No read carried a usable CB cell barcode; cell counts are all zero")


def cell_umi_allele_evidence(isovar_results, alignment_file, *, sample_id, source, read_collector=None):
    """
    RNA support records with cell/UMI counts for each result's alleles and proteins.

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
    warn_if_unlabelled(support for e in evidence for support in e["alleles"].values())
    return evidence
