"""Supplied somatic and matched germline variants, for attributing assembled edits.

An assembled cDNA can contain differences from the reference transcript besides
the nominated variant. Isovar attributes such an edit to a supplied variant only
when the edit is that variant, on the same transcript. It never guesses an
origin from allele fraction.
"""

from bisect import bisect_left, bisect_right


def _variant_key(variant):
    return (variant.contig, variant.start, variant.ref, variant.alt)


class KnownVariants(object):
    """
    Supplied somatic variants and an optional matched germline collection.

    Parameters
    ----------
    somatic : iterable of varcode.Variant
        The nominated somatic variants, usually every variant given to
        `run_isovar`. Edits matching one other than a result's own variant
        are co-somatic.
    germline : iterable of varcode.Variant
        Variants from a matched normal sample.

    Raises
    ------
    ValueError
        If a variant is in both collections, or they use different reference
        assemblies.
    """

    def __init__(self, somatic=(), germline=()):
        self.somatic = tuple(sorted(set(somatic), key=_variant_key))
        self.germline = tuple(sorted(set(germline), key=_variant_key))
        both = set(self.somatic) & set(self.germline)
        if both:
            raise ValueError("Variants cannot be both somatic and germline: %s" % ", ".join(
                sorted(v.short_description for v in both)))
        references = {v.reference_name for v in self.somatic + self.germline}
        if len(references) > 1:
            raise ValueError("Somatic and germline variants use different references: %s" % ", ".join(
                sorted(str(r) for r in references)))
        self._origins = dict.fromkeys(self.somatic, "somatic")
        self._origins.update(dict.fromkeys(self.germline, "germline"))
        by_contig = {}
        for variant in self.somatic + self.germline:
            by_contig.setdefault(variant.contig, []).append(variant)
        self._index = {}
        for contig, variants in by_contig.items():
            variants.sort(key=lambda v: v.start)
            self._index[contig] = ([v.start for v in variants], variants,
                                   max(max(len(v.ref), 1) for v in variants))
        self._hash = hash((self.somatic, self.germline))

    def origin(self, variant):
        """``"somatic"``, ``"germline"`` or None for a variant not supplied."""
        return self._origins.get(variant)

    def overlapping(self, contig, start, end):
        """Supplied variants touching 1-based positions ``start`` to ``end``."""
        if contig not in self._index:
            return []
        starts, variants, span = self._index[contig]
        return variants[bisect_left(starts, start - span):bisect_right(starts, end)]

    def __len__(self):
        return len(self._origins)

    def __eq__(self, other):
        return isinstance(other, KnownVariants) and (self.somatic, self.germline) == (
            other.somatic, other.germline)

    def __hash__(self):
        return self._hash

    def __repr__(self):
        return "KnownVariants(somatic=%d, germline=%d)" % (len(self.somatic), len(self.germline))
