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

import numpy as np


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

def create_dictionary_of_filter_results(
        isovar_result,
        filter_thresholds):
    """
    Creates a dictionary whose keys are named of different
    filter conditions and values are booleans, where True
    indicates whether this set of coverage stats passes
    the filter and False indicates that it failed.

    Parameters
    ----------
    isovar_result : IsovarResult

    filter_thresholds : dict or OrderedDict
        Every argument is supposed to be something like "max_alt_reads"
        where the first three characters are "min" or "max" and the
        rest of the name is either a field of IsovarResult or
        a numeric field like "num_alt_reads". The name of each filter
        maps to a cutoff value. Filters starting with "max"
        require that the corresponding field on CoverageStats
        is <= cutoff, whereas filters starting with
        "min" require >= cutoff.

    Returns
    -------
    Dictionary of filter names mapped to boolean value indicating
    whether this locus passed the filter.
    """
    filter_values_dict = OrderedDict()
    for name, threshold in filter_thresholds.items():
        parts = name.split("_")
        min_or_max = parts[0]
        field_name = "_".join(parts[1:])
        if min_or_max == "min":
            comparison_fn = operator.ge
        elif min_or_max == "max":
            comparison_fn = operator.le
        else:
            raise ValueError(
                "Invalid filter '%s', must start with 'min' or 'max'" % name)
        if hasattr(isovar_result, field_name):
            field_value = getattr(isovar_result, field_name)
        else:
            raise ValueError(
                "Invalid filter '%s' IsovarResult does not have property '%s'" % (
                    name,
                    field_name))
        filter_values_dict[name] = comparison_fn(field_value, threshold)
    return filter_values_dict

def apply_filters(
        result_generator,
        filter_thresholds=DEFAULT_FILTER_THRESHOLDS)
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


def apply_filters(self, **kwargs):
    """
    Apply dictionary of cutoff values to fields of this CoverageStats
    object and return enough information to quickly determine
    whether all filters passed and which, if any, failed.

    Parameters
    ----------
    **kwargs : dict
        Dictionary of filter names (e.g. "max_alt_reads") and
        their corresponding cutoffs. The filter names can apply
        to any filter of CoverageStats which either starts with
        "n_" (like "n_alt_reads") or ends with "fraction"
        (like "alt_fragment_fraction").

    Returns
    -------
    A tuple with two elements:
        - boolean indicating whether all filters passed
        - list of names of which filters failed
    """
    filter_name_to_value = self.create_dictionary_of_filter_results(**kwargs)
    all_passed = True
    failing_filters = []
    for (filter_name, filter_value) in filter_name_to_value.items():
        if not filter_value:
            all_passed = False
            failing_filters.append(filter_name)
    return all_passed, failing_filters

