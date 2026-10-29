"""Versioned gene nomination and shared antigen reference preparation."""

from collections.abc import Iterable
from dataclasses import dataclass
from importlib.metadata import version
import json
from pathlib import Path
from typing import Any
import warnings

import pysam
from pyensembl import Genome

from .antigen_reference import (
    GDC_REFERENCE, IsoformTarget, _annotation_cache, gdc_genome,
    isoform_targets_from_genome, reference_file_identity, validate_alignment_reference,
)

CTA_GENE_SET = "oncoref-cta"


@dataclass(frozen=True)
class AntigenReference:
    """Prepared transcript competitors and their selection/reference provenance.

    Attributes
    ----------
    targets : tuple of IsoformTarget
        Detached transcript models usable after the annotation is closed.
    identity : dict
        Alignment/reference validation and local FASTA/GTF file identities.
    target_selection : dict
        Selected genes, provider decisions, sources and annotation mapping.
        Attach this ledger to the reconstruction result to retain provenance.
    """

    targets: tuple[IsoformTarget, ...]
    identity: dict[str, Any]
    target_selection: dict[str, Any]


def antigen_gene_selection(*, genes: Iterable[str] | None = None,
                           gene_set: str | None = None) -> dict[str, Any]:
    """Resolve exactly one explicit gene list or a versioned CTA catalogue.

    Parameters
    ----------
    genes : iterable of str, optional
        Explicit stable gene IDs or symbols. Symbols are resolved only against
        the supplied annotation during reference preparation.
    gene_set : str, optional
        ``oncoref-cta`` selects all human gene loci in OncoRef's CTA category.
        This research nomination does not apply Vaxrank admission exclusions.

    Returns
    -------
    dict
        JSON-compatible gene selection, provider identity, decisions and source
        inventory. The installed OncoRef snapshot is used without downloads.
    """
    if (genes is None) == (gene_set is None):
        raise ValueError("Provide exactly one of genes or gene_set")
    if genes is not None:
        if isinstance(genes, str):
            raise TypeError("genes must be an iterable of gene IDs/names, not a string")
        selected = tuple(genes)
        if not selected or any(not isinstance(g, str) or not g.strip() for g in selected):
            raise ValueError("genes must contain nonempty gene IDs/names")
        return dict(kind="explicit", genes=sorted(set(selected)))
    if gene_set != CTA_GENE_SET:
        raise ValueError("Unknown antigen gene set: %r" % gene_set)
    from oncoref import cta_selection

    coverage = cta_selection.coverage(species="human")
    if coverage["selection_status"] != "available":
        raise ValueError("OncoRef's human CTA selection is unavailable")
    rows = json.loads(cta_selection.results(species="human", category="CTA").to_json(orient="records"))
    selected = sorted({gene for row in rows for gene in row["gene_ids"]})
    if not selected:
        raise ValueError("OncoRef's human CTA selection contains no gene loci")
    return dict(kind="catalogue", gene_set=gene_set, genes=selected,
                provider=dict(name="oncoref", version=version("oncoref"),
                              api="cta_selection.results(species='human', category='CTA')"),
                coverage=coverage, decisions=rows,
                sources=json.loads(cta_selection.sources(species="human").to_json(orient="records")))


def _annotation_genes(
        genome: Genome, selection: dict[str, Any]) -> tuple[list[str], list[dict[str, Any]]]:
    """Resolve IDs exactly; only explicitly supplied symbols permit name lookup."""
    resolved, decisions = set(), []
    for query in selection["genes"]:
        try:
            matches = [genome.gene_by_id(query)]
            method = "gene_id"
        except ValueError:
            matches = genome.genes_by_name(query) if selection["kind"] == "explicit" else []
            method = "gene_name" if matches else "absent_from_annotation"
        ids = sorted({gene.id for gene in matches})
        resolved.update(ids)
        decisions.append(dict(query=query, gene_ids=ids, status=method))
    missing = [row["query"] for row in decisions if not row["gene_ids"]]
    if missing:
        warnings.warn("Antigen genes absent from the supplied annotation: " + ", ".join(missing),
                      UserWarning, stacklevel=2)
    if not resolved:
        raise ValueError("No selected antigen genes occur in the supplied annotation")
    return sorted(resolved), decisions


def prepare_antigen_reference(
        bam: pysam.AlignmentFile, *, genome_fasta: str | Path, gtf: str | Path,
        annotation_version: str | int, genes: Iterable[str] | None = None,
        gene_set: str | None = None, annotation_name: str = "GENCODE",
        assembly: str = "GRCh38", alignment_reference: str = GDC_REFERENCE,
        reference_cache: str | Path | None = None) -> AntigenReference:
    """Load all target transcripts using one explicit annotation/reference pair.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Open sample alignment; its sequence dictionary is validated.
    genome_fasta, gtf : str or Path
        Local alignment FASTA and transcript/CDS annotation; hashes are retained.
    annotation_version : str or int
        Annotation release, independent of the alignment reference profile.
    genes, gene_set : optional
        Exactly one selection, as described by :func:`antigen_gene_selection`.
    annotation_name, assembly, alignment_reference : str
        Annotation identity and either the GDC profile or ``custom``. The human
        OncoRef catalogue accepts GRCh37/GRCh38 stable IDs without liftover.
    reference_cache : str or Path, optional
        Content-addressed annotation cache root, shared across sample sources.

    Returns
    -------
    AntigenReference
        Detached models, reference identity and full selection/mapping ledger.
        Missing loci and unavailable transcripts remain explicit. Reference
        models are competitors; reconstruction does not fill RNA gaps with them.
    """
    selection = antigen_gene_selection(genes=genes, gene_set=gene_set)
    if gene_set is not None and assembly not in {"GRCh37", "GRCh38"}:
        raise ValueError("The OncoRef human CTA catalogue requires GRCh37 or GRCh38")
    if alignment_reference == GDC_REFERENCE:
        if assembly != "GRCh38":
            raise ValueError("The GDC profile uses GRCh38; use custom for other assemblies")
        genome = gdc_genome(genome_fasta=genome_fasta, gtf=gtf,
                            annotation_version=annotation_version, annotation_name=annotation_name,
                            cache_directory=reference_cache)
    elif alignment_reference == "custom":
        genome = Genome(assembly, annotation_name, annotation_version,
                        gtf_path_or_url=str(Path(gtf).expanduser().resolve()),
                        genome_fasta_path_or_url=str(Path(genome_fasta).expanduser().resolve()),
                        cache_directory_path=_annotation_cache(gtf, directory=reference_cache))
    else:
        raise ValueError("Unknown alignment reference profile: %r" % alignment_reference)
    with genome:
        genome.index()
        selected, mapping = _annotation_genes(genome, selection)
        targets, unavailable = isoform_targets_from_genome(genome, selected, contig_names=bam.references)
        if not targets:
            raise ValueError("No selected antigen transcripts have usable reference models")
        identity = validate_alignment_reference(
            bam, genome, profile=None if alignment_reference == "custom" else GDC_REFERENCE,
            contigs=[target.reference.contig for target in targets])
        identity["files"] = reference_file_identity(genome_fasta, gtf)
    selection.update(resolved_gene_ids=selected, gene_mapping=mapping, unavailable_transcripts=unavailable)
    return AntigenReference(tuple(targets), identity, selection)
