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

Distinct from :class:`isovar.VarcodeAdapter`, which implements varcode's
full ``IsovarAssemblyProvider`` protocol (assembled contigs +
reference-keyed mutant transcripts). This adapter answers only the
phasing question: "are these two variants observed on the same RNA
reads?" — derived from
:attr:`isovar.IsovarResult.phased_variants_in_supporting_reads`.
"""


class IsovarReadPhasing(object):
    """
    Indexes a collection of IsovarResult objects by variant and answers
    read-phasing queries via Isovar's per-variant
    ``phased_variants_in_supporting_reads``.

    Example
    -------
    >>> from isovar import run_isovar, IsovarReadPhasing
    >>> results = run_isovar(
    ...     variants="tumor.vcf",
    ...     alignment_file="tumor.rna.bam")
    >>> phasing = IsovarReadPhasing(results)
    >>> phasing.has_evidence(some_variant)
    True
    >>> phasing.partners_in_cis(some_variant)
    (Variant(...), ...)
    """

    def __init__(self, isovar_results):
        """
        Parameters
        ----------
        isovar_results : Iterable[IsovarResult]
            Results from a finished Isovar run (e.g. the output of
            ``run_isovar``). The iterable is consumed once.
        """
        self._by_variant = {
            result.variant: result
            for result in isovar_results
        }

    def has_evidence(self, variant):
        """
        True iff Isovar saw at least one alt-supporting RNA fragment
        for ``variant``.

        Returns False both when the variant was not in the input set
        and when it was present but had no alt fragments.
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
        partners = result.phased_variants_in_supporting_reads
        return tuple(sorted(partners, key=_variant_sort_key))


def _variant_sort_key(variant):
    return (variant.contig, variant.start, variant.ref, variant.alt)
