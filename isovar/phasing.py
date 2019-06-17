# Copyright (c) 2019. Mount Sinai School of Medicine
#
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


def find_phased_variants(
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
    read_names_to_variants = defaultdict(set)
    variants_to_read_names = {}

    # first associate every read name with the variants which use that read
    # to assemble to assemble a protein sequence
    for isovar_result in unphased_isovar_results:
        variant = isovar_result.variant
        read_names = isovar_result.alt_read_names
        variants_to_read_names[variant] = read_names
        for read_name in read_names:
            read_names_to_variants[read_name].add(variant)

    # now count up how many reads are shared between pairs of variants
    phasing_counts = defaultdict(Counter)
    for variant, read_names in variants_to_read_names.items():
        for read_name in read_names:
            for other_variant in read_names_to_variants[read_name]:
                if variant != other_variant:
                    phasing_counts[variant][other_variant] += 1

    results_with_phasing = []
    for isovar_result in unphased_isovar_results:
        variant = isovar_result.variant
        other_variant_counts = phasing_counts[variant]
        phased_variants = {
            other_variant
            for (other_variant, count)
            in other_variant_counts.items()
            if count >= min_shared_fragments_for_phasing
        }
        results_with_phasing.append(
            isovar_result.clone_with_updates(phased_variants=phased_variants)
        )
    return results_with_phasing
