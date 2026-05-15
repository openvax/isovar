# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Adapter that exposes Isovar's per-variant RNA-read co-occurrence as a
narrow read-phasing interface.

Implements the ``has_evidence(variant)`` + ``partners_in_cis(variant)``
shape that varcode's ``ReadPhasingSource`` Protocol adopts. An instance
of this class plugs directly into ``varcode.ReadPhaseResolver``. The
companion ``MutantTranscriptSource`` concern (assembled-cDNA-derived
mutant transcripts) is intentionally not covered here -- it's a
separate axis, served by :class:`isovar.IsovarMutantTranscript`.

Answers only the phasing question: "are these two variants observed on
the same RNA reads?" -- derived from
:attr:`isovar.IsovarResult.phased_variants_in_supporting_reads`.

Contract
--------
``IsovarReadPhasing`` reads exactly three attributes from each element
of its input iterable:

* ``variant`` -- a hashable identity used as the index key.
* ``num_alt_fragments`` -- int; non-zero means
  :meth:`IsovarReadPhasing.has_evidence` returns ``True``.
* ``phased_variants_in_supporting_reads`` -- iterable of variants
  co-observed on supporting RNA reads.

Any object exposing those three attributes is a valid input; the
adapter is intentionally duck-typed so tests, mocks, and alternative
RNA-phasing producers (e.g. long-read tools) can target the same shape.

Co-occurrence symmetry is inherited from the underlying field: Isovar's
``compute_phasing_counts`` builds counts symmetrically and thresholds
both directions with the same minimum, so
``v2 in phasing.partners_in_cis(v1)`` implies
``v1 in phasing.partners_in_cis(v2)`` whenever both variants are in the
input set.

``partners_in_cis`` and ``has_evidence`` answer independent questions and
are not cross-checked: ``partners_in_cis`` returns
``phased_variants_in_supporting_reads`` as-is from the input results, so
a caller constructing inputs by hand could observe ``has_evidence(v)``
return ``False`` while ``partners_in_cis(v)`` returns a non-empty tuple.
Outputs from ``run_isovar``/``annotate_phased_variants`` never exhibit
this -- the upstream pipeline only populates partners for variants that
have shared supporting reads, which implies positive alt-fragment counts
on both sides.
"""

from .isovar_result_provider import IsovarResultProvider


class IsovarReadPhasing(IsovarResultProvider):
    """
    Index a collection of IsovarResult-shaped objects by variant and
    answer RNA-read phasing queries from
    ``phased_variants_in_supporting_reads``.

    Examples
    --------
    >>> from types import SimpleNamespace
    >>> from varcode import Variant
    >>> v1 = Variant("1", 10, "A", "C", normalize_contig_names=False)
    >>> v2 = Variant("1", 11, "G", "T", normalize_contig_names=False)
    >>> results = [
    ...     SimpleNamespace(
    ...         variant=v1,
    ...         num_alt_fragments=3,
    ...         phased_variants_in_supporting_reads={v2}),
    ...     SimpleNamespace(
    ...         variant=v2,
    ...         num_alt_fragments=2,
    ...         phased_variants_in_supporting_reads={v1}),
    ... ]
    >>> phasing = IsovarReadPhasing(results)
    >>> phasing.has_evidence(v1)
    True
    >>> phasing.partners_in_cis(v1) == (v2,)
    True
    """

    def __init__(self, isovar_results, *args, **kwargs):
        """
        Parameters
        ----------
        isovar_results : Iterable[IsovarResult]
            Results from a finished Isovar run (e.g. the output of
            ``run_isovar``). The iterable is consumed once.
        """
        isovar_results = tuple(isovar_results)
        super().__init__(isovar_results, *args, **kwargs)
        self._by_variant = {
            result.variant: result
            for result in isovar_results
        }

    def __repr__(self):
        return "IsovarReadPhasing(%d variants)" % len(self._by_variant)

    def has_evidence(self, variant):
        """
        True iff Isovar saw at least one alt-supporting RNA fragment
        for ``variant``.

        Returns ``False`` both when the variant was not in the input
        set and when it was present but had no alt fragments.
        """
        result = self._by_variant.get(variant)
        return result is not None and result.num_alt_fragments > 0

    def partners_in_cis(self, variant):
        """
        Variants co-observed on supporting RNA reads with ``variant``.

        Returns an empty tuple when ``variant`` is unknown or has no
        co-observed partners. Partners are returned in a deterministic
        order (by chromosome, start, ref, alt) so callers can rely on
        ordering for snapshot tests.
        """
        result = self._by_variant.get(variant)
        if result is None:
            return ()
        return tuple(sorted(
            result.phased_variants_in_supporting_reads,
            key=_variant_sort_key))


def _variant_sort_key(variant):
    return (variant.contig, variant.start, variant.ref, variant.alt)
