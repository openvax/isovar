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
from six import integer_types, text_type, binary_type

from varcode import Variant
import pandas as pd

VALID_ELEMENT_TYPES = integer_types + (text_type, binary_type, float, bool)

class DataFrameBuilder(object):
    """
    Helper class for constructing a DataFrame which always has fields
    of a variant (chr/pos/ref/alt) as well as some subset of the fields
    from a namedtuple.
    """
    def __init__(
            self,
            element_class,
            exclude=set([]),
            converters={},
            rename_dict={},
            variant_columns=True,
            convert_collections_to_size=True):
        """
        Parameters
        ----------
        element_class : type
            Expected to a have a class-member named '_fields' which is a list
            of field names.

        exclude : set
            Field names from element_class which should be used as columns for
            the DataFrame we're building

        converters : dict
            Dictionary of names mapping to functions. These functions will be
            applied to each element of a column before it's added to the
            DataFrame.

        rename_dict : dict
            Dictionary mapping element_class field names to desired column names
            in the produced DataFrame.

        variant_columns : bool
            If True, then add four columns for fields of a Variant: chr/pos/ref/alt

        convert_collections_to_size : bool
            If a value is a built-in collection (list, tuple, or set) then
            transform it to the size of that collection. If this option is False
            then collection values cause a runtime error.
        """
        self.element_class = element_class
        self.rename_dict = rename_dict
        self.converters = converters
        self.variant_columns = variant_columns
        self.convert_collections_to_size = convert_collections_to_size

        # remove specified field names without changing the order of the others
        self.original_field_names = [
            x
            for x in element_class._fields
            if x not in exclude
        ]

        for name in converters:
            if name not in self.original_field_names:
                raise ValueError("No field named '%s', valid names: %s" % (
                    name,
                    self.original_field_names))

        self.renamed_field_names = [
            self.rename_dict.get(x, x)
            for x in self.original_field_names
        ]
        if self.variant_columns:
            columns_list = [
                # fields related to variant
                ("chr", []),
                ("pos", []),
                ("ref", []),
                ("alt", []),
            ]
        else:
            columns_list = []

        for name in self.renamed_field_names:
            columns_list.append((name, []))

        self.columns_dict = OrderedDict(columns_list)

    def add(self, variant, element):
        if self.variant_columns:
            assert isinstance(variant, Variant)
            self.columns_dict["chr"].append(variant.contig)
            self.columns_dict["pos"].append(variant.original_start)
            self.columns_dict["ref"].append(variant.original_ref)
            self.columns_dict["alt"].append(variant.original_alt)
        else:
            assert variant is None

        assert isinstance(element, self.element_class)

        for name in self.original_field_names:
            value = getattr(element, name)

            if name in self.converters:
                fn = self.converters[name]
                value = fn(value)

            if isinstance(value, (tuple, list, set)) and self.convert_collections_to_size:
                value = len(value)
            elif not isinstance(value, VALID_ELEMENT_TYPES):
                raise ValueError(
                    "Please provider converter for field '%s' : %s to make a scalar or string" % (
                        name,
                        type(value)))

            if name in self.rename_dict:
                name = self.rename_dict[name]
            self.columns_dict[name].append(value)

    def add_many(self, variant, elements):
        for element in elements:
            self.add(variant, element)

    def to_dataframe(self):
        return pd.DataFrame(self.columns_dict)

def dataframe_from_generator(element_class, variant_and_elements_generator, **kwargs):
    builder = DataFrameBuilder(element_class, **kwargs)
    for variant, elements in variant_and_elements_generator:
        builder.add_many(variant, elements)
    return builder.to_dataframe()
