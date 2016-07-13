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

# values of these types are automatically converted to their size or length
# unless some other conversion function is provided
COLLECTION_TYPES = (tuple, list, set, frozenset)

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
            extra_column_fns={},
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

        extra_column_fns : dict
            Dictionary mapping column names to functions which take a variant
            and element (such as an AlleleRead instance) and return a single
            value for each row.

        variant_columns : bool
            If True, then add four columns for fields of a Variant: chr/pos/ref/alt
            along with a "gene" column indicating which gene name(s) the variant
            overlaps.

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

        self.extra_column_fns = extra_column_fns
        for column_name in self.extra_column_fns:
            columns_list.append((column_name, []))

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

            if isinstance(value, COLLECTION_TYPES) and self.convert_collections_to_size:
                value = len(value)
            elif not isinstance(value, VALID_ELEMENT_TYPES):
                raise ValueError(
                    "Please provider converter for field '%s' : %s to make a scalar or string" % (
                        name,
                        type(value)))

            if name in self.rename_dict:
                name = self.rename_dict[name]
            self.columns_dict[name].append(value)

        for column_name, fn in self.extra_column_fns.items():
            self.columns_dict[column_name].append(fn(variant, element))

    def add_many(self, variant, elements):
        for element in elements:
            self.add(variant, element)

    def _check_column_lengths(self):
        """
        Make sure columns are of the same length or else DataFrame construction
        will fail.
        """
        column_lengths_dict = {
            name: len(xs)
            for (name, xs)
            in self.columns_dict.items()
        }
        unique_column_lengths = set(column_lengths_dict.values())
        if len(unique_column_lengths) != 1:
            raise ValueError(
                "Mismatch between lengths of columns: %s" % (column_lengths_dict,))

    def to_dataframe(self):
        self._check_column_lengths()
        return pd.DataFrame(self.columns_dict)

def dataframe_from_generator(element_class, variant_and_elements_generator, **kwargs):
    builder = DataFrameBuilder(element_class, **kwargs)
    for variant, elements in variant_and_elements_generator:
        builder.add_many(variant, elements)
    return builder.to_dataframe()
