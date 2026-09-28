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

from collections import defaultdict
from itertools import combinations
from math import exp, lgamma, log, log1p

from .default_parameters import (
    MAX_P_VALUE_FOR_PHASING,
    MIN_SHARED_FRAGMENTS_FOR_PHASING,
    PHASING_ERROR_RATE,
)
from .phase_group import PhaseGroup
from .read_identity import alignment_constraints, fragment_ids
from .chimeric_alignment import compatible_phasing_alignments
from .transcript_edit_helpers import transcript_assembly_edit_sort_key


def _binomial_tail(count, trials, rate):
    """P(X >= count) for X ~ Binomial(trials, rate), with 0 < rate < 1."""
    if count <= 0:
        return 1.0
    total, scale = 0.0, lgamma(trials + 1)
    for successes in range(count, trials + 1):
        term = exp(scale - lgamma(successes + 1) - lgamma(trials - successes + 1)
                   + successes * log(rate) + (trials - successes) * log1p(-rate))
        total += term
        # Past the mean the terms only shrink.
        if successes > trials * rate and term <= 1e-17 * total:
            break
    return min(total, 1.0)


def _allele_reads(result, allele):
    """``result``'s reads for ``"ref"`` or ``"alt"``; names for legacy inputs."""
    reads = getattr(result, allele + "_reads", None)
    if reads is None:
        reads = getattr(result, allele + "_read_names", ())
    return reads


def _validate_phasing_rates(phasing_error_rate, max_p_value_for_phasing):
    for name, value in (("phasing_error_rate", phasing_error_rate),
                        ("max_p_value_for_phasing", max_p_value_for_phasing)):
        if not 0 < value < 1:
            raise ValueError("%s must be between 0 and 1, not %r" % (name, value))


def _four_gamete(
        both, first, second, neither,
        min_shared_fragments_for_phasing=MIN_SHARED_FRAGMENTS_FOR_PHASING,
        phasing_error_rate=PHASING_ERROR_RATE,
        max_p_value_for_phasing=MAX_P_VALUE_FOR_PHASING):
    """Cis, trans or ``None`` for two variants in the run, from their fragment table."""
    def present(count, neighbour, other_neighbour):
        # Enough fragments, and more than errors at either locus leak in
        # from the two neighbouring combinations' fragments together.
        return (count >= min_shared_fragments_for_phasing and _binomial_tail(
            count, count + neighbour + other_neighbour, phasing_error_rate)
            <= max_p_value_for_phasing)
    together = present(both, first, second)
    first_alone, second_alone = present(first, neither, both), present(second, neither, both)
    if together and not (first_alone and second_alone):
        return True
    if first_alone and second_alone and not together:
        return False
    return None


def _variant_sort_key(variant):
    return (
        variant.contig,
        variant.start,
        variant.ref,
        variant.alt,
    )


def create_variant_to_top_protein_sequence_dict(isovar_results):
    """
    Create dictionary from variant to its top translated protein sequence.

    Variants without an assembled protein sequence are mapped to None.
    """
    return {
        isovar_result.variant: (
            isovar_result.top_protein_sequence
            if isovar_result.has_mutant_protein_sequence_from_rna
            else None
        )
        for isovar_result in isovar_results
    }


def _phasing_support(variant_to_reads):
    """Build pairwise edges from compatible observations, once per fragment.

    Alternative observations are hypotheses, not simultaneous constraints.
    A pair needs at least one jointly compatible choice. Reuse the assembly
    rule for linear placements and complementary mates. Different placements
    of one segment require reciprocal, compatible supplementary-path evidence.

    Metadata-free reads/IDs retain legacy name-only semantics. Public names
    are kept separately and never substituted for scoped evidence IDs.
    """
    observations = defaultdict(lambda: defaultdict(list))
    names_by_id = {}
    for variant, reads in variant_to_reads.items():
        for read in reads:
            constraints = alignment_constraints((read,))
            for fragment in fragment_ids((read,)):
                observations[fragment][variant].append((constraints, read))
                names_by_id[fragment] = getattr(read, "name", read)
    support = defaultdict(dict)
    for fragment, by_variant in observations.items():
        placements = defaultdict(set)
        for reads in by_variant.values():
            for constraints, _ in reads:
                for segment, placement in constraints.items():
                    placements[segment].add(placement)
        for (variant, reads), (other, other_reads) in combinations(by_variant.items(), 2):
            if any(compatible_phasing_alignments(constraints, read, other_read, placements)
                   for constraints, read in reads for _, other_read in other_reads):
                shared = support[variant].setdefault(other, set())
                shared.add(fragment)
                support[other][variant] = shared
    return support, names_by_id


def _allele_table(v1, reads1, v2, reads2):
    """Unambiguous fragment counts: both alt, first only, second only, neither."""
    support, _ = _phasing_support({
        (allele, variant): _allele_reads(reads, allele)
        for variant, reads in ((v1, reads1), (v2, reads2)) for allele in ("alt", "ref")
    })
    cells = [support[(a, v1)].get((b, v2), set())
             for a, b in (("alt", "alt"), ("alt", "ref"), ("ref", "alt"), ("ref", "ref"))]
    # Unmerged, disagreeing mates can support multiple combinations. A
    # physical fragment has only one combination; keep these observations
    # in the results but exclude them from phase counts. Consider all four
    # cells even when the caller only uses the somatic-alt row.
    seen, ambiguous = set(), set()
    for cell in cells:
        ambiguous.update(seen & cell)
        seen.update(cell)
    return tuple(len(cell - ambiguous) for cell in cells)


def _phase_annotations(
        variant_to_reads,
        min_shared_fragments_for_phasing,
        variant_to_top_protein_sequence_dict=None,
        in_cis=None):
    """Derive neighbors and groups from the same validated fragment edges.

    ``variant_to_reads`` maps each variant to its reads, or to plain read
    names for callers without alignment metadata. Returns the phased
    neighbors of every variant and the PhaseGroup of each grouped variant.
    When supplied, ``in_cis`` must confirm each edge against the complete
    reference/alternate allele evidence.
    """
    support, names_by_id = _phasing_support(variant_to_reads)
    phased_neighbors = {
        variant: {other for other, ids in support[variant].items()
                  if len(ids) >= min_shared_fragments_for_phasing
                  and (in_cis is None or in_cis(variant, other))}
        for variant in variant_to_reads
    }

    visited = set()
    variant_to_phase_group = {}
    for variant in sorted(variant_to_reads, key=_variant_sort_key):
        if variant in visited:
            continue

        component = set()
        pending = [variant]
        while pending:
            current_variant = pending.pop()
            if current_variant in component:
                continue
            component.add(current_variant)
            visited.add(current_variant)
            pending.extend(phased_neighbors.get(current_variant, set()) - component)

        if len(component) <= 1:
            continue

        supporting_read_names = {
            names_by_id[fragment]
            for grouped_variant in component
            for other in phased_neighbors[grouped_variant]
            for fragment in support[grouped_variant][other]
        }

        cdna_sequences = set()
        mutant_protein_sequences = set()
        transcript_ids = set()
        transcript_names = set()
        known_somatic_transcript_edits = set()
        known_germline_transcript_edits = set()
        unexplained_transcript_edits = set()
        for grouped_variant in component:
            protein_sequence = (variant_to_top_protein_sequence_dict or {}).get(grouped_variant)
            if protein_sequence is None:
                continue
            if hasattr(protein_sequence, "_transcript_assembly_edits_by_category"):
                categorized_edits = protein_sequence._transcript_assembly_edits_by_category()
            else:
                categorized_edits = {
                    "known_somatic": getattr(protein_sequence, "known_somatic_transcript_edits", ()),
                    "known_germline": getattr(protein_sequence, "known_germline_transcript_edits", ()),
                    "unexplained": getattr(protein_sequence, "unexplained_transcript_edits", ()),
                }
            cdna_sequences.update(protein_sequence.cdna_sequences)
            mutant_protein_sequences.add(protein_sequence.amino_acids)
            transcript_ids.update(protein_sequence.transcript_ids)
            transcript_names.update(protein_sequence.transcript_names)
            known_somatic_transcript_edits.update(categorized_edits["known_somatic"])
            known_germline_transcript_edits.update(categorized_edits["known_germline"])
            unexplained_transcript_edits.update(categorized_edits["unexplained"])

        phase_group = PhaseGroup(
            somatic_variants=tuple(sorted(component, key=_variant_sort_key)),
            # Germline variants whose edits the group's assemblies contain.
            germline_variants=tuple(sorted(
                {edit.source_variant for edit in known_germline_transcript_edits
                 if edit.source_variant is not None},
                key=_variant_sort_key)),
            supporting_read_names=supporting_read_names,
            cdna_sequences=tuple(sorted(cdna_sequences)),
            mutant_protein_sequences=tuple(sorted(mutant_protein_sequences)),
            transcript_ids=tuple(sorted(transcript_ids)),
            transcript_names=tuple(sorted(transcript_names)),
            known_somatic_transcript_edits=tuple(sorted(
                known_somatic_transcript_edits,
                key=transcript_assembly_edit_sort_key,
            )),
            known_germline_transcript_edits=tuple(sorted(
                known_germline_transcript_edits,
                key=transcript_assembly_edit_sort_key,
            )),
            unexplained_transcript_edits=tuple(sorted(
                unexplained_transcript_edits,
                key=transcript_assembly_edit_sort_key,
            )),
        )
        for grouped_variant in component:
            variant_to_phase_group[grouped_variant] = phase_group
    return phased_neighbors, variant_to_phase_group


def _variant_reads(isovar_results, protein=False):
    """Retain placement evidence; accept names-only legacy result objects."""
    by_variant = {}
    for result in isovar_results:
        if protein:
            if result.has_mutant_protein_sequence_from_rna:
                sequence = result.top_protein_sequence
                reads = getattr(sequence, "supporting_reads", None)
                if reads is None:
                    reads = sequence.read_names_supporting_protein_sequence
            else:
                reads = ()
        else:
            reads = getattr(result, "alt_reads", None)
            if reads is None:
                reads = result.alt_read_names
        by_variant[result.variant] = reads
    return by_variant


def annotate_phased_variants(
        unphased_isovar_results,
        min_shared_fragments_for_phasing=MIN_SHARED_FRAGMENTS_FOR_PHASING,
        phasing_error_rate=PHASING_ERROR_RATE,
        max_p_value_for_phasing=MAX_P_VALUE_FOR_PHASING):
    """
    Annotate IsovarResult objects with phasing information. Pairs must pass
    the same error-aware test as IsovarReadPhasing.in_cis, using all reference
    and alternate reads. Protein-sequence edges additionally require enough
    shared protein-supporting fragments with compatible placements, scoped
    by SAM read group. Public read-name sets remain display
    names, not evidence IDs. Legacy caller-created objects without alignment
    metadata retain name-only phasing;
    they are not equated with scoped collected fragments.

    Parameters
    ----------
    unphased_isovar_results : list of IsovarResult

    min_shared_fragments_for_phasing : int
        Minimum shared fragments for an edge.

    phasing_error_rate : float
        Probability of an incorrect allele at either locus.

    max_p_value_for_phasing : float
        Maximum binomial p-value for an allele combination to exceed errors.

    Returns
    -------
    list of IsovarResult
    """

    _validate_phasing_rates(phasing_error_rate, max_p_value_for_phasing)
    unphased_isovar_results = tuple(unphased_isovar_results)
    # Assess every pair against the full allele evidence, including reference
    # reads. Selecting only protein-supporting reads would discard the error
    # background and could turn two sequencing errors into a cis edge.
    by_variant = {result.variant: result for result in unphased_isovar_results}
    calls = {}

    def in_cis(v1, v2):
        key = frozenset((v1, v2))
        if key not in calls:
            calls[key] = _four_gamete(
                *_allele_table(v1, by_variant[v1], v2, by_variant[v2]),
                min_shared_fragments_for_phasing=min_shared_fragments_for_phasing,
                phasing_error_rate=phasing_error_rate,
                max_p_value_for_phasing=max_p_value_for_phasing)
        return calls[key] is True

    updates = {result.variant: {} for result in unphased_isovar_results}
    for source in ("supporting_reads", "protein_sequence"):
        protein = source == "protein_sequence"
        neighbors, groups = _phase_annotations(
            _variant_reads(unphased_isovar_results, protein=protein),
            min_shared_fragments_for_phasing,
            create_variant_to_top_protein_sequence_dict(unphased_isovar_results) if protein else None,
            in_cis=in_cis)
        for variant, fields in updates.items():
            fields["phased_variants_in_" + source] = neighbors[variant]
            fields["phase_group_from_" + source] = groups.get(variant)
    return [result.clone_with_updates(**updates[result.variant])
            for result in unphased_isovar_results]
