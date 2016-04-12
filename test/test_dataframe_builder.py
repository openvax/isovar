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
from collections import namedtuple, OrderedDict
from nose.tools import eq_
from isovar.dataframe_builder import DataFrameBuilder
from varcode import Variant
import pandas as pd

def test_dataframe_builder():
    TestClass = namedtuple("TestClass", "a b c")
    df_builder = DataFrameBuilder(TestClass)
    test_obj = TestClass(a=1, b="s", c=3.0)
    test_variant = Variant("X", 10, "CC", "C")
    df_builder.add(test_variant, test_obj)
    df = df_builder.to_dataframe()
    expected = pd.DataFrame(OrderedDict([
        ("chr", ["X"]),
        ("pos", [10]),
        ("ref", ["CC"]),
        ("alt", ["C"]),
        ("a", [test_obj.a]),
        ("b", [test_obj.b]),
        ("c", [test_obj.c]),
    ]))
    eq_(len(df.columns), len(expected.columns))
    assert all(x == y for (x, y) in zip(df.columns, expected.columns)), \
        (df.columns, expected.columns)
    # need to aggregate twice because the first time just rolls up each column,
    # giving us a boolean Series
    assert (df == expected).all().all(), "Expected %s == %s" % (df, expected)
