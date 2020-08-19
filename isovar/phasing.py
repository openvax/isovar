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

from __future__ import print_function, division, absolute_import

from collections import defaultdict, Counter

from .default_parameters import MIN_SHARED_FRAGMENTS_FOR_PHASING


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
    read_names_to_variants = defaultdict(set)

    for variant, read_names in variant_to_read_names_dict.items():
        for read_name in read_names:
            read_names_to_variants[read_name].add(variant)

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

    # create dictionary counting how often each variant co-occurs with others
    # in any reads supporting those variants
    phasing_counts_from_supporting_reads = compute_phasing_counts(
        create_variant_to_alt_read_names_dict(unphased_isovar_results))

    # create dictionary counting how often each variant co-occurs with others
    # in any reads used to construct their protein sequences
    phasing_counts_from_protein_sequences = compute_phasing_counts(
        create_variant_to_protein_sequence_read_names_dict(
            unphased_isovar_results))

    results_with_phasing = []
    for isovar_result in unphased_isovar_results:
        variant = isovar_result.variant
        phased_variants_in_supporting_reads = threshold_phased_variant_counts(
            phasing_counts_from_supporting_reads[variant],
            min_count=min_shared_fragments_for_phasing)
        phased_variants_in_protein_sequence = \
            threshold_phased_variant_counts(
                phasing_counts_from_protein_sequences[variant],
                min_count=min_shared_fragments_for_phasing)
        results_with_phasing.append(
            isovar_result.clone_with_updates(
                phased_variants_in_supporting_reads=phased_variants_in_supporting_reads,
                phased_variants_in_protein_sequence=phased_variants_in_protein_sequence))
    return results_with_phasing
