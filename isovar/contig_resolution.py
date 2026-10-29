"""Resolve input names against immutable alignment reference vocabularies.

Varcode owns spelling rules. This module owns alignment-context lifetime and
the distinction between an unresolved name and a successfully queried locus.
Name resolution does not establish sequence identity or perform liftover.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import asdict
import os

from pysam import AlignmentFile
from varcode import ContigResolution, ContigResolver, Variant


def variant_reference_name(variant: Variant) -> str | None:
    """Return assembly spelling context without opening annotation data."""
    genome = getattr(variant, "genome", None)
    return getattr(genome, "assembly", None) or getattr(variant, "reference_name", None)


def variant_source_contig(variant: Variant) -> str | int:
    """Retain the supplied name so exact source matches precede aliases."""
    return getattr(variant, "original_contig", variant.contig)


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
