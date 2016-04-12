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
            exclude_field_names=set([])):
        self.element_class = element_class

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

    def add(self, variant, element):
        self.columns_dict["chr"].append(variant.contig)
        self.columns_dict["pos"].append(variant.original_start)
        self.columns_dict["ref"].append(variant.original_ref)
        self.columns_dict["alt"].append(variant.original_alt)

        for name in self.field_names:
            self.columns_dict[name].append(getattr(element, name))

    def to_dataframe(self):
        return pd.DataFrame(self.columns_dict)
