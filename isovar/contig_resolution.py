"""Resolve input names against immutable alignment reference vocabularies.

Varcode owns spelling rules. This module excludes known incompatible reference
sequences and owns alignment-context lifetime and unresolved-name diagnostics.
Name resolution does not establish sequence identity or perform liftover.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import asdict, replace
import os

from pysam import AlignmentFile
from varcode import ContigResolution, ContigResolver, Variant, normalize_contig_name


def variant_reference_name(variant: Variant) -> str | None:
    """Return assembly spelling context without opening annotation data."""
    genome = getattr(variant, "genome", None)
    return getattr(genome, "assembly", None) or getattr(variant, "reference_name", None)


def variant_source_contig(variant: Variant) -> str | int:
    """Retain the supplied name so exact source matches precede aliases."""
    return getattr(variant, "original_contig", variant.contig)


def _grch37_mitochondrial_group(name: str | int, reference_name: str | None) -> str | None:
    if (reference_name or "").strip().lower() in {"grch37", "hg19", "b37", "ncbi37"}:
        return {"M": "chrM", "CHRM": "chrM", "MT": "MT", "CHRMT": "MT"}.get(str(name).upper())
    return None


def reference_contig_key(name: str | int, *, reference_name: str | None = None) -> str:
    """Normalize allele identity without merging GRCh37's two mitochondria.

    Parameters
    ----------
    name : str or int
        Source contig spelling.
    reference_name : str, optional
        Assembly context. GRCh37/hg19/b37/NCBI37 keeps M/chrM separate from
        MT/chrMT; other names follow Varcode's spelling normalization.

    Returns
    -------
    str
        Comparison key, not a resolved target name or proof of sequence identity.
    """
    normalized = normalize_contig_name(name, reference_name=reference_name)
    return _grch37_mitochondrial_group(name, reference_name) or normalized


def resolve_reference_contig(resolver: ContigResolver, name: str | int) -> ContigResolution:
    """Resolve coordinate names without aliasing distinct GRCh37 mitochondria.

    Parameters
    ----------
    resolver : varcode.ContigResolver
        Reusable vocabulary index with the input's assembly context.
    name : str or int
        Source contig spelling. An exact target match always takes precedence.

    Returns
    -------
    varcode.ContigResolution
        Exact, alias, ambiguous, missing or reference-unavailable result. For
        GRCh37/hg19/b37/NCBI37, aliases stay within M/chrM (NC_001807) or
        MT/chrMT (rCRS, NC_012920), case-insensitively. Other builds retain
        Varcode's spelling rules. No sequence validation or liftover occurs.

    Notes
    -----
    UCSC documents the distinct hg19 sequences at
    https://genome.ucsc.edu/FAQ/FAQgenes.html. A GRCh37 chrM
    already using rCRS coordinates needs an explicitly reconciled source name;
    a name or length alone does not establish sequence identity.
    """
    result = resolver.resolve(name)
    group = _grch37_mitochondrial_group(name, resolver.reference_name)
    if group is None:
        return result
    result = replace(result, normalized_name=group)
    if result.status not in ("alias", "ambiguous"):
        return result
    candidates = tuple(candidate for candidate in result.candidates
                       if _grch37_mitochondrial_group(candidate, resolver.reference_name) == group)
    status = "alias" if len(candidates) == 1 else "ambiguous" if candidates else "missing"
    return replace(result, candidates=candidates, status=status,
                   resolved_name=candidates[0] if status == "alias" else None)


class UnresolvedContigError(ValueError):
    """An alignment contig was missing, ambiguous or unavailable.

    Attributes
    ----------
    resolution : varcode.ContigResolution
        Supplied name, reference context, status and ambiguous candidates.
    source : str or None
        Alignment filename when available. This is not an RNA assessment.
    """

    def __init__(self, resolution: ContigResolution, source: str | None = None) -> None:
        self.resolution = resolution
        self.source = source
        suffix = " (alignment %r)" % source if source is not None else ""
        super().__init__("unresolved_chromosome: %s%s" % (resolution.diagnostic, suffix))

    def to_dict(self) -> dict[str, object]:
        """Return JSON-ready diagnostics for an unassessed input record."""
        return dict(asdict(self.resolution), source=self.source)


def require_resolved_contig(resolution: ContigResolution, alignment_file: AlignmentFile | None = None) -> str:
    """Return an actual target name or raise without implying zero reads."""
    if resolution.resolved_name is None:
        filename = getattr(alignment_file, "filename", None)
        source = None if filename is None else os.fsdecode(filename)
        raise UnresolvedContigError(resolution, source)
    return resolution.resolved_name


class AlignmentContigCache:
    """Reuse Varcode indexes for up to eight active header/assembly contexts.

    Pysam headers are immutable and have stable identity, including when two
    open files have different vocabularies. The cache retains headers, not file
    handles. Reopening a file supplies a new header and therefore a new index.
    Header-less reader adapters use their object identity and must expose an
    immutable ``references`` vocabulary for their lifetime. The bounded cache
    prevents a long-lived collector from retaining every historical header.
    """

    def __init__(self) -> None:
        self._indexes: OrderedDict[tuple[int, str | None], tuple[object, ContigResolver]] = OrderedDict()

    def resolver(self, alignment_file: AlignmentFile, reference_name: str | None) -> ContigResolver:
        """Index the vocabulary once per active input/reference context."""
        if getattr(alignment_file, "closed", False):
            raise ValueError("Cannot resolve contigs from a closed alignment file")
        header = getattr(alignment_file, "header", None)
        identity = alignment_file if header is None else header
        key = (id(identity), reference_name)
        if key in self._indexes:
            self._indexes.move_to_end(key)
            return self._indexes[key][1]
        resolver = ContigResolver(getattr(alignment_file, "references", None), reference_name=reference_name)
        self._indexes[key] = (identity, resolver)
        if len(self._indexes) > 8:
            self._indexes.popitem(last=False)
        return resolver
