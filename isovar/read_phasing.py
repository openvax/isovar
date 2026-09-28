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

Answers only the phasing question. ``partners_in_cis`` lists variants
observed on the same RNA fragments (from
:attr:`isovar.IsovarResult.phased_variants_in_supporting_reads`), which
establishes cis. Not being a partner is not evidence of trans: Isovar
only compares variants in its run, so a variant outside it was never
examined. :meth:`IsovarReadPhasing.in_cis` therefore answers cis *and*
trans from fragments that cover both loci, including matched germline
variants whose reads ``run_isovar`` collected, and varcode's
``MolecularPhaseResolver`` uses it in preference to the partner lists.

Contract
--------
``IsovarReadPhasing`` reads three attributes from each element of its
input iterable:

* ``variant`` -- a hashable identity used as the index key.
* ``num_alt_fragments`` -- int; non-zero means
  :meth:`IsovarReadPhasing.has_evidence` returns ``True``.
* ``phased_variants_in_supporting_reads`` -- iterable of variants
  co-observed on supporting RNA reads.

:meth:`IsovarReadPhasing.in_cis` also reads ``alt_reads`` and ``ref_reads``
(or ``*_read_names``), and optionally ``germline_read_evidence`` and
``germline_variants_in_top_protein_sequence``.

Any object exposing those attributes is a valid input; the
adapter is intentionally duck-typed so tests, mocks, and alternative
RNA-phasing producers (e.g. long-read tools) can target the same shape.

Co-occurrence symmetry is inherited from the underlying field: Isovar's
phasing counts shared fragments symmetrically and thresholds both
directions with the same minimum, so
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

from .default_parameters import (
    MAX_P_VALUE_FOR_PHASING,
    MIN_SHARED_FRAGMENTS_FOR_PHASING,
    PHASING_ERROR_RATE,
)
from .isovar_result_provider import IsovarResultProvider
from .phasing import (
    _allele_table,
    _four_gamete,
    _validate_phasing_rates,
    _variant_sort_key,
)


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

    def __init__(
            self,
            isovar_results,
            *args,
            min_shared_fragments_for_phasing=MIN_SHARED_FRAGMENTS_FOR_PHASING,
            phasing_error_rate=PHASING_ERROR_RATE,
            max_p_value_for_phasing=MAX_P_VALUE_FOR_PHASING,
            **kwargs):
        """
        Parameters
        ----------
        isovar_results : Iterable[IsovarResult]
            Results from a finished Isovar run (e.g. the output of
            ``run_isovar``). The iterable is consumed once.

        min_shared_fragments_for_phasing : int
            Fragments needed for :meth:`in_cis` to call cis or trans; the
            same default as :func:`isovar.run_isovar`.

        phasing_error_rate : float
            Chance that a read shows the wrong allele at a variant's locus,
            from sequencing error, RNA editing or mismapping. Raise it for
            error-prone reads, such as ONT's.

        max_p_value_for_phasing : float
            One-sided binomial p-value at or below which :meth:`in_cis`
            counts a combination of two variants' alleles as more than such
            errors; about the chance, per pair, of a false cis or trans when
            ``phasing_error_rate`` is the true error rate.
        """
        _validate_phasing_rates(phasing_error_rate, max_p_value_for_phasing)
        isovar_results = tuple(isovar_results)
        super().__init__(isovar_results, *args, **kwargs)
        self._by_variant = {
            result.variant: result
            for result in isovar_results
        }
        self.min_shared_fragments_for_phasing = min_shared_fragments_for_phasing
        self.phasing_error_rate = phasing_error_rate
        self.max_p_value_for_phasing = max_p_value_for_phasing
        self._calls = {}

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

    def in_cis(self, v1, v2, transcript=None):
        """
        Whether ``v1`` and ``v2`` are on the same RNA molecules.

        For two variants in the run, the fragments with compatible
        placements that cover both loci carry one of four combinations: both
        alt alleles, either variant's alt allele alone, or neither. A
        combination is present when it's on at least
        ``min_shared_fragments_for_phasing`` fragments, more than reads
        showing the wrong allele at either locus leak into it from its two
        neighbouring combinations: a one-sided binomial test at
        ``phasing_error_rate`` and ``max_p_value_for_phasing``. The pair is
        cis when both alt alleles are present together and not each alone,
        and trans (never on one molecule) when each is present alone and not
        together. Otherwise, including when all three are present, the
        answer is ``None``.

        This is the four-gamete test (Hudson and Kaplan 1985) applied to
        read-phased variant pairs, as Nik-Zainal et al. 2012 (Cell 149:994)
        read them: reads with the earlier variant alone and with both mean
        nesting; reads with either variant alone and never both mean
        separate molecules (another copy, or a sibling subclone). A variant
        that arose in a subclone of the other's cells, on the other's copy,
        has all its alt molecules cis with the other's alt allele; the
        earlier variant's alt allele alone is not trans (#393). Errors put a
        variant's alt allele on a small fraction of the other's alt
        fragments, which is not cis (#410). All three combinations present
        fit no single lineage.

        For a variant in the run and a matched germline variant
        (``run_isovar(germline_variants=...)``) that its alt reads cover,
        only fragments carrying the somatic alt allele count: with the
        germline alt allele they are cis, with its reference allele trans.
        A fragment with the somatic reference allele says nothing, since the
        germline alt allele also comes from cells without the somatic
        variant and from both copies of a homozygous variant. The germline
        reads are in ``IsovarResult.germline_read_evidence``. When the reads
        do not decide, a germline edit in the variant's top assembled
        protein is cis; when they say trans but the protein has the edit,
        the answer is ``None``.

        Returns ``True`` or ``False`` as above, and ``None`` otherwise,
        including when no fragment covers both loci. For a germline variant,
        a count must reach ``min_shared_fragments_for_phasing`` and exceed
        the other.

        ``transcript`` is accepted for varcode's phase-resolver interface;
        the evidence does not depend on the isoform.
        """
        key = tuple(sorted((v1, v2), key=_variant_sort_key))
        if key not in self._calls:
            self._calls[key] = self._in_cis(v1, v2)
        return self._calls[key]

    def _in_cis(self, v1, v2):
        r1 = self._by_variant.get(v1)
        r2 = self._by_variant.get(v2)
        if r1 is None and r2 is None:
            return None
        if r1 is not None and r2 is not None:
            return _four_gamete(
                *_allele_table(v1, r1, v2, r2),
                min_shared_fragments_for_phasing=self.min_shared_fragments_for_phasing,
                phasing_error_rate=self.phasing_error_rate,
                max_p_value_for_phasing=self.max_p_value_for_phasing)
        result, germline = (r1, v2) if r2 is None else (r2, v1)
        evidence = (getattr(result, "germline_read_evidence", None) or {}).get(germline)
        call = None
        if evidence is not None:
            # Only the somatic alt allele's fragments: with the germline alt allele, or its reference.
            cis, trans = _allele_table(result.variant, result, germline, evidence)[:2]
            call = self._decide(cis, trans)
        in_protein = germline in getattr(result, "germline_variants_in_top_protein_sequence", ())
        if in_protein:
            return None if call is False else True
        return call

    def _decide(self, cis, trans):
        if cis >= self.min_shared_fragments_for_phasing and cis > trans:
            return True
        if trans >= self.min_shared_fragments_for_phasing and trans > cis:
            return False
        return None
