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
Functions used to annotate IsovarResult objects with filters.
"""

from __future__ import print_function, division, absolute_import

from collections import OrderedDict
import operator


def evaluate_threshold_filters(isovar_result, filter_thresholds):
    """
    Helper method used by apply_filters

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

    Returns OrderedDict
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


def evaluate_boolean_filters(isovar_result, filter_flags):
    """
    Helper function used by apply_filters.

    Parameters
    ----------
    isovar_result : IsovarResult

    filter_flags : list of str
        Every element should be a boolean property of IsovarResult
        or "not_" and the name of a property to be negated.

    Returns OrderedDict
    """
    filter_values = OrderedDict()
    for boolean_filter_name in filter_flags:
        if boolean_filter_name.startswith("not_"):
            boolean_field_name = boolean_filter_name[4:]
            negate = True
        else:
            boolean_field_name = boolean_filter_name
            negate = False
        if hasattr(isovar_result, boolean_field_name):
            field_value = getattr(isovar_result, boolean_field_name)
        else:
            raise ValueError(
                "IsovarResult does not have field name '%s'" % boolean_field_name)
        if field_value is None:
            field_value = False
        elif field_value not in {True, False}:
            raise ValueError("Expected filter '%s' to be boolean but got %s" % (
                boolean_filter_name,
                field_value))
        filter_values[boolean_filter_name] = (
            not field_value if negate else field_value
        )
    return filter_values


def evaluate_filters(
        isovar_result,
        filter_thresholds,
        filter_flags=[]):
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

    filter_flags : list of str
        Every element should be a boolean property of IsovarResult
        or "not_" and the name of a property to be negated.

    Returns
    -------
    Dictionary of filter names mapped to boolean value indicating
    whether this locus passed the filter.
    """
    filter_values_dict = evaluate_boolean_filters(isovar_result, filter_flags)
    filter_values_dict.update(
        evaluate_threshold_filters(isovar_result, filter_thresholds))
    return filter_values_dict


def apply_filters(
        isovar_result,
        filter_thresholds={},
        filter_flags=[]):
    """
    Given an IsovarResult object, evaluates given filters
    for each object, and returns a copy of the IsovarResult with new fiter
    values.

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

    filter_flags : list of str
        Every element should be a boolean property of IsovarResult
        or "not_" and the name of a property to be negated.

    Returns IsovarResult
    """
    filter_values = OrderedDict(isovar_result.filter_values.items())
    new_filter_values = evaluate_filters(
        isovar_result,
        filter_thresholds=filter_thresholds,
        filter_flags=filter_flags)
    filter_values.update(new_filter_values)
    return isovar_result.clone_with_updates(filter_values=filter_values)
