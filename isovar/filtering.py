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

"""
CoverageStats object is used to collect allow relevant
counts and fractions used to filter variants by RNA evidence.
"""

from __future__ import print_function, division, absolute_import

import operator

from collections import OrderedDict

from .default_parameters import (
    MIN_NUM_RNA_ALT_READS,
    MIN_NUM_RNA_ALT_FRAGMENTS,
    MIN_FRACTION_RNA_ALT_READS,
    MIN_FRACTION_RNA_ALT_FRAGMENTS,
    MAX_NUM_RNA_REF_READS,
    MAX_NUM_RNA_REF_FRAGMENTS,
    MAX_FRACTION_RNA_REF_READS,
    MAX_FRACTION_RNA_REF_FRAGMENTS,
    MAX_NUM_RNA_OTHER_READS,
    MAX_NUM_RNA_OTHER_FRAGMENTS,
    MAX_FRACTION_RNA_OTHER_READS,
    MAX_FRACTION_RNA_OTHER_FRAGMENTS,
    MIN_RATIO_RNA_ALT_TO_OTHER_FRAGMENTS,
)

DEFAULT_FILTER_THRESHOLDS =  OrderedDict([
    # alt allele
    ("min_num_alt_reads", MIN_NUM_RNA_ALT_READS),
    ("min_num_alt_fragments", MIN_NUM_RNA_ALT_FRAGMENTS),
    ("min_fraction_alt_reads", MIN_FRACTION_RNA_ALT_READS),
    ("min_fraction_alt_fragments", MIN_FRACTION_RNA_ALT_FRAGMENTS),

    # ref allele coverage and VAF
    ("max_num_ref_reads", MAX_NUM_RNA_REF_READS),
    ("max_num_ref_fragments", MAX_NUM_RNA_REF_FRAGMENTS),
    ("max_fraction_ref_reads", MAX_FRACTION_RNA_REF_READS),
    ("max_fraction_ref_fragments", MAX_FRACTION_RNA_REF_FRAGMENTS),

    # other alleles
    ("max_num_other_reads", MAX_NUM_RNA_OTHER_READS),
    ("max_num_other_fragments", MAX_NUM_RNA_OTHER_FRAGMENTS),
    ("max_fraction_other_reads", MAX_FRACTION_RNA_OTHER_READS),
    ("max_fraction_other_fragments", MAX_FRACTION_RNA_OTHER_FRAGMENTS),

    # misc. filters
    ("min_ratio_alt_to_other_fragments", MIN_RATIO_RNA_ALT_TO_OTHER_FRAGMENTS)
])




def apply_filters(
        result_generator,
        filter_thresholds=DEFAULT_FILTER_THRESHOLDS):
    """
    Given a generator of IsovarResult objects and thresholds for various 
    filters, this function generates dictionaries of filter values for each
    result along with a boolean value indicating whether any of the filters
    failed.
    
    Parameters
    ----------
    result_generator : generator of IsovarResult
    
    filter_thresholds : dict or OrderedDict
        Dictionary of filter names such as "min_num_alt_fragments" with 
        associated threshold values. The filter names must follow the
        pattern:
            {min|max}_{num|fraction}_{ref|alt|other}_{reads|fragments}
        Filters whose first component is "min" are filtered using the >=
        operator whereas those which start with "max" are filtered using
        the <= operator. 
        There is currently one special case which doesn't fit this 
        pattern, "min_ratio_alt_to_other_fragments", which is a threshold
        on the ratio of the number of alt fragments to non-ref/non-alt fragments.
    Generator of (IsovarResult, OrderedDict, bool) tuples
    """
    for isovar_result in result_generator:
        filter_dict = create_dictionary_of_filter_results(
            isovar_result=isovar_result,
            filter_thresholds=filter_thresholds)
        all_passed = all(filter_dict.values())
        yield (isovar_result, filter_dict, all_passed)


def split_by_filters(
        result_generator,
        filter_thresholds=DEFAULT_FILTER_THRESHOLDS):
    """
    Split a series of IsovarResult objects into two lists, the first
    is the results which passed all the filters, whereas the second
    list is failing results. Each IsovarResult is paired with a dictionary
    mapping each filter name to whether it passed or failed.

    Parameters
    ----------
    result_generator : generator of IsovarResult

    filter_thresholds : dict
        Names such as "min_num_alt_reads" mapped to cutoffs. Everything
        after "min_" or "max_" is expected to be a property of the IsovarResult
        object.
    """
    passing_list = []
    failing_list = []
    for (isovar_result, filter_result_dict, all_passed) in apply_filters(
            result_generator, filter_thresholds):
        pair = (isovar_result, filter_result_dict)
        if all_passed:
            passing_list.append(pair)
        else:
            failing_list.append(pair)
    return passing_list, failing_list
