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

from collections import OrderedDict
import numpy as np

from .common import safediv
from .value_object import ValueObject


class CoverageStats(ValueObject):
    """
    Read and fragment counts for ref, alt, and other alleles at this locus,
    which can be used for filtering this variant/locus based on criteria such
    as "min_alt_fragments", "max_other_to_alt_fragment_ratio", &c.
    """
    __slots__  = [
        "n_ref_reads",
        "n_ref_fragments",
        "n_alt_reads",
        "n_alt_fragments",
        "n_other_reads",
        "n_other_fragments"
    ]

    @classmethod
    def from_reads(cls, ref_reads, alt_reads, other_reads):
        """
        Create a CoverageStats object from sets of AlleleRead objects
        corresponding to the ref, alt, and other alleles.

        Parameters
        ----------
        ref_reads : set or list of AlleleRead
        alt_reads : set or list of AlleleRead
        other_reads : set or list of AlleleRead

        Returns CoverageStats

        """
        n_ref_reads = len(ref_reads)
        n_alt_reads = len(alt_reads)
        n_other_reads = len(other_reads)
        n_ref_fragments = len({r.name for r in ref_reads})
        n_alt_fragments = len({r.name for r in alt_reads})
        n_other_fragments = len({r.name for r in other_reads})
        return cls(
            n_ref_reads=n_ref_reads,
            n_ref_fragments=n_ref_fragments,
            n_alt_reads=n_alt_reads,
            n_alt_fragments=n_alt_fragments,
            n_other_reads=n_other_reads,
            n_other_fragments=n_other_fragments)

    @property
    def n_reads(self):
        """
        Total read coverage at this site.
        """
        return self.n_ref_reads + self.n_alt_reads + self.n_other_reads


    @property
    def n_fragments(self):
        """
        Number of fragments supporting reference allele.
        """
        return self.n_ref_fragments + self.n_alt_fragments + self.n_other_fragments


    @property
    def ref_read_fraction(self):
        """
        Allelic fraction of the reference allele among all reads at this site.
        """
        return self.n_ref_reads / self.n_reads


    @property
    def alt_fragment_fraction(self):
        """
        Allelic fraction of the reference allele among all fragments at this site.
        """
        return self.n_ref_fragments / self.n_fragments


    @property
    def alt_read_fraction(self):
        """
        Allelic fraction of the variant allele among all reads at this site.
        """
        return self.n_alt_reads / self.n_reads


    @property
    def alt_fragment_fraction(self):
        """
        Allelic fraction of the variant allele among all fragments at this site.
        """
        return self.n_alt_fragments / self.n_fragments


    @property
    def other_read_fraction(self):
        """
        Allelic fraction of the "other" (non-ref, non-alt) alleles among all
        reads at this site.
        """
        return self.n_other_reads / self.n_reads


    @property
    def other_fragment_fraction(self):
        """
        Allelic fraction of the "other" (non-ref, non-alt) alleles among all
        reads at this site.
        """
        return self.n_other_fragments / self.n_fragments


    @property
    def other_to_ref_read_ratio(self):
        return safediv(self.n_other_reads, self.n_ref_reads)


    @property
    def other_to_alt_read_ratio(self):
        return safediv(self.n_other_reads, self.n_alt_reads)


    @property
    def other_to_alt_fragment_ratio(self):
        return safediv(self.n_other_reads, self.n_ref_reads)


    @property
    def other_to_ref_fragment_ratio(self):
        return safediv(self.n_other_reads, self.n_ref_reads)


    # List of which properties get included in the result of stats()
    # Only using these since others can be derived from these quantities.
    _extra_dict_fields = [
        "n_reads",
        "n_fragments",
        "other_to_alt_read_ratio",
        "other_to_alt_fragment_ratio",
    ]


    def to_dict(self):
        """
        Returns OrderedDict containing all fields of thsi CoverageStats object.
        """
        result = OrderedDict()
        for field in list(self.__slots__) + self._extra_dict_fields:
            result[field] = getattr(self, field)
        return result


    def create_dictionary_of_filter_results(
            self,
            **kwargs):
        """
        Creates a dictionary whose keys are named of different
        filter conditions and values are booleans, where True
        indicates whether this set of coverage stats passes
        the filter and False indicates that it failed.

        Parameters
        ----------
        **kwargs : dict
            Every argument is supposed to be something like "max_alt_reads"
            where the first three characters are "min" or "max" and the
            rest of the name is either a field of CoverageStats or
            a count like "n_alt_reads". The name of each filter
            maps to a cutoff value. Filters starting with "max"
            require that the corresponding field on CoverageStats
            is <= cutoff, whereas filters starting with
            "min" require >= cutoff.

        Returns
        -------
        Dictionary of filter names mapped to boolean value indicating
        whether this locus passed the filter.
        """
        result = {}
        for name, cutoff in kwargs.items():
            parts = name.split("_")
            min_or_max = parts[0]
            field_name = "_".join(parts[1:])
            if min_or_max not in {"min", "max"}:
                raise ValueError(
                    "Invalid filter name: '%s', must start with 'min' or 'max'" % name)

            if "reads" in field_name or "fragments" in field_name:
                field_name = "n_" + field_name
            elif "fraction" not in field_name:
                raise ValueError("Invalid filter name: '%s'" % name)
            self_value = getattr(self, field_name)
            if min_or_max == "min":
                filter_value = (self_value >= cutoff)
            else:
                filter_value = (self_value <= cutoff)
            result[name] = filter_value
        return result

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

# these are the arguments used as defaults for
# CoverageStats.create_dictionary_of_filter_values
# but since we also want to expose the same arguments
# on the commandline, they're pulled out into their
# own dictionary


COVERAGE_FILTER_DEFAULT_VALUES = OrderedDict([
    ("min_reads", 0),
    ("min_fragments", 0),
    ("min_alt_reads", 0),
    ("min_alt_fragments", 0),
    ("max_ref_reads", np.inf),
    ("max_ref_fragments", np.inf),
    ("max_other_reads", np.inf),
    ("max_other_fragments", np.inf),
    ("min_alt_read_fraction", 0.0),
    ("min_alt_fragment_fraction", 0.0),
    ("max_ref_read_fraction", 1.0),
    ("max_ref_fragment_fraction", 1.0),
    ("max_other_read_fraction", 1.0),
    ("max_other_fragment_fraction", 1.0)])