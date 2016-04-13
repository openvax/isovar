# Copyright (c) 2016. Mount Sinai School of Medicine
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
from collections import OrderedDict

import pandas as pd

class DataFrameBuilder(object):
    """
    Helper class for constructing a DataFrame which always has fields
    of a variant (chr/pos/ref/alt) as well as some subset of the fields
    from a namedtuple.
    """
    def __init__(
            self,
            element_class,
            exclude_field_names=set([]),
            transform_fields={},
            convert_sequences_to_string=True):
        """
        Parameters
        ----------
        element_class : type
            Expected to a have a class-member named '_fields' which is a list
            of field names.

        exclude_field_names : set
            Field names from element_class which should be used as columns for
            the DataFrame we're building

        transform_fields : dict
            Dictionary of names mapping to functions. These functions will be
            applied to each element of a column before it's added to the
            DataFrame.

        convert_sequences_to_string : bool
            If a value being added to a column is a list or tuple, should
            it be converted to a semi-colon separated string?
        """
        self.element_class = element_class
        self.convert_sequences_to_string = convert_sequences_to_string
        # remove specified field names without changing the order of the others
        self.field_names = [
            x
            for x in element_class._fields
            if x not in exclude_field_names
        ]
        columns_list = [
            # fields related to variant
            ("chr", []),
            ("pos", []),
            ("ref", []),
            ("alt", []),
        ]

        for name in self.field_names:
            columns_list.append((name, []))
        self.columns_dict = OrderedDict(columns_list)
        self.transform_fields = transform_fields
        for name in transform_fields:
            if name not in self.columns_dict:
                raise ValueError("No field named '%s' in %s" % (
                    name, element_class))

    def add(self, variant, element):
        self.columns_dict["chr"].append(variant.contig)
        self.columns_dict["pos"].append(variant.original_start)
        self.columns_dict["ref"].append(variant.original_ref)
        self.columns_dict["alt"].append(variant.original_alt)

        for name in self.field_names:
            value = getattr(element, name)

            if name in self.transform_fields:
                fn = self.transform_fields[name]
                value = fn(value)

            if isinstance(value, (list, tuple)):
                if self.convert_sequences_to_string:
                    value = ";".join(value)
            self.columns_dict[name].append(value)

    def to_dataframe(self):
        return pd.DataFrame(self.columns_dict)
