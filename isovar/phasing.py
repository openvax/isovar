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

from collections import defaultdict, Counter

from .default_parameters import MIN_SHARED_FRAGMENTS_FOR_PHASING
from .phase_group import PhaseGroup


def _variant_sort_key(variant):
    return (
        variant.contig,
        variant.start,
        variant.ref,
        variant.alt,
    )


def create_read_names_to_variants_dict(variant_to_read_names_dict):
    """
    Invert a variant -> read-name mapping into read-name -> variants.
    """
    read_names_to_variants = defaultdict(set)

    for variant, read_names in variant_to_read_names_dict.items():
        for read_name in read_names:
            read_names_to_variants[read_name].add(variant)
    return read_names_to_variants


def create_variant_to_alt_read_names_dict(isovar_results):
    """
    Create dictionary from variant to names of alt reads supporting
    that variant in an IsovarResult.

    Parameters
    ----------
    isovar_results : list of IsovarResult

    Returns
    -------
    Dictionary from varcode.Variant to set(str) of read names

    """
    return {
        isovar_result.variant: set(isovar_result.alt_read_names)
        for isovar_result in isovar_results
    }


def create_variant_to_protein_sequence_read_names_dict(isovar_results):
    """
    Create dictionary from variant to names of alt reads used to create
    mutant protein sequence in an IsovarResult.

    Parameters
    ----------
    isovar_results : list of IsovarResult

    Returns
    -------
    Dictionary from varcode.Variant to set(str) of read names
    """
    variant_to_read_names = {}
    for isovar_result in isovar_results:
        if isovar_result.has_mutant_protein_sequence_from_rna:
            protein_sequence = isovar_result.top_protein_sequence
            read_names = set(
                protein_sequence.read_names_supporting_protein_sequence)
        else:
            read_names = set()
        variant_to_read_names[isovar_result.variant] = read_names
    return variant_to_read_names


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


def compute_phasing_counts(variant_to_read_names_dict):
    """

    Parameters
    ----------
    variants_to_read_names : dict
        Dictionary mapping varcode.Variant to set of read names

    Returns
    -------
    Dictionary from variant to Counter(Variant)
    """
    read_names_to_variants = create_read_names_to_variants_dict(
        variant_to_read_names_dict
    )

    # now count up how many reads are shared between pairs of variants
    phasing_counts = defaultdict(Counter)
    for variant, read_names in variant_to_read_names_dict.items():
        for read_name in read_names:
            for other_variant in read_names_to_variants[read_name]:
                if variant != other_variant:
                    phasing_counts[variant][other_variant] += 1
    return phasing_counts

def threshold_phased_variant_counts(counts_dict, min_count):
    """
    Choose set of phased variants by keeping any variants with associated
    counts greater than or equal the given threshold.

    Parameters
    ----------
    counts_dict : variant -> int dict

    min_count : int

    Returns
    -------
    set of varcode.Variant
    """
    return {
        variant
        for (variant, count)
        in counts_dict.items()
        if count >= min_count
    }


def create_phase_groups(
        variant_to_read_names_dict,
        min_shared_fragments_for_phasing,
        variant_to_top_protein_sequence_dict=None):
    """
    Group variants into connected components of the phasing graph.

    If top translated protein sequences are provided then each resulting
    PhaseGroup is also annotated with directly observed cDNA, protein, and
    transcript metadata from those assemblies.

    Returns
    -------
    dict
        Mapping from variant to PhaseGroup. Variants without phased partners are
        omitted.
    """
    phasing_counts = compute_phasing_counts(variant_to_read_names_dict)
    phased_neighbors = {
        variant: threshold_phased_variant_counts(
            phasing_counts[variant],
            min_count=min_shared_fragments_for_phasing)
        for variant in variant_to_read_names_dict
    }

    read_names_to_variants = create_read_names_to_variants_dict(
        variant_to_read_names_dict
    )

    visited = set()
    variant_to_phase_group = {}
    for variant in sorted(variant_to_read_names_dict, key=_variant_sort_key):
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
            read_name
            for read_name, read_variants in read_names_to_variants.items()
            if len(component.intersection(read_variants)) >= 2
        }

        if variant_to_top_protein_sequence_dict is None:
            cdna_sequences = ()
            mutant_protein_sequences = ()
            transcript_ids = ()
            transcript_names = ()
        else:
            cdna_sequences = set()
            mutant_protein_sequences = set()
            transcript_ids = set()
            transcript_names = set()
            for grouped_variant in component:
                protein_sequence = variant_to_top_protein_sequence_dict.get(
                    grouped_variant
                )
                if protein_sequence is None:
                    continue
                cdna_sequences.update(protein_sequence.cdna_sequences)
                mutant_protein_sequences.add(protein_sequence.amino_acids)
                transcript_ids.update(protein_sequence.transcript_ids)
                transcript_names.update(protein_sequence.transcript_names)

        phase_group = PhaseGroup(
            somatic_variants=tuple(sorted(component, key=_variant_sort_key)),
            germline_variants=(),
            supporting_read_names=supporting_read_names,
            cdna_sequences=tuple(sorted(cdna_sequences)),
            mutant_protein_sequences=tuple(sorted(mutant_protein_sequences)),
            transcript_ids=tuple(sorted(transcript_ids)),
            transcript_names=tuple(sorted(transcript_names)),
        )
        for grouped_variant in component:
            variant_to_phase_group[grouped_variant] = phase_group
    return variant_to_phase_group


def annotate_phased_variants(
        unphased_isovar_results,
        min_shared_fragments_for_phasing=MIN_SHARED_FRAGMENTS_FOR_PHASING):
    """
    Annotate IsovarResult objects with phasing information. Phasing
    is determined by looking at RNA fragments used for assembled protein
    sequences and counting the number of shared fragments.

    Parameters
    ----------
    unphased_isovar_results : list of IsovarResult

    min_shared_fragments_for_phasing : int

    Returns
    -------
    list of IsovarResult
    """

    phasing_counts_from_supporting_reads = compute_phasing_counts(
        create_variant_to_alt_read_names_dict(unphased_isovar_results)
    )

    phase_groups_from_supporting_reads = create_phase_groups(
        create_variant_to_alt_read_names_dict(unphased_isovar_results),
        min_shared_fragments_for_phasing=min_shared_fragments_for_phasing,
    )

    phasing_counts_from_protein_sequences = compute_phasing_counts(
        create_variant_to_protein_sequence_read_names_dict(
            unphased_isovar_results)
    )

    phase_groups_from_protein_sequences = create_phase_groups(
        create_variant_to_protein_sequence_read_names_dict(
            unphased_isovar_results),
        min_shared_fragments_for_phasing=min_shared_fragments_for_phasing,
        variant_to_top_protein_sequence_dict=create_variant_to_top_protein_sequence_dict(
            unphased_isovar_results),
    )

    results_with_phasing = []
    for isovar_result in unphased_isovar_results:
        variant = isovar_result.variant
        phase_group_from_supporting_reads = phase_groups_from_supporting_reads.get(
            variant
        )
        phased_variants_in_supporting_reads = threshold_phased_variant_counts(
            phasing_counts_from_supporting_reads[variant],
            min_count=min_shared_fragments_for_phasing,
        )

        phase_group_from_protein_sequence = phase_groups_from_protein_sequences.get(
            variant
        )
        phased_variants_in_protein_sequence = threshold_phased_variant_counts(
            phasing_counts_from_protein_sequences[variant],
            min_count=min_shared_fragments_for_phasing,
        )
        results_with_phasing.append(
            isovar_result.clone_with_updates(
                phase_group_from_supporting_reads=phase_group_from_supporting_reads,
                phase_group_from_protein_sequence=phase_group_from_protein_sequence,
                phased_variants_in_supporting_reads=phased_variants_in_supporting_reads,
                phased_variants_in_protein_sequence=phased_variants_in_protein_sequence))
    return results_with_phasing
