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
from os.path import join, dirname

from nose.tools import eq_
from pysam import AlignmentFile
import varcode


def assert_equal_fields(result, expected):
    """
    Assert that fields of two namedtuple objects have same field values.
    """
    eq_(result.__class__, expected.__class__)
    for field in result.__class__._fields:
        result_value = getattr(result, field)
        expected_value = getattr(expected, field)
        assert result_value == expected_value, \
            "Wrong value for '%s', expected %s but got %s" % (
                field,
                expected_value,
                result_value)


def data_path(path):
    """
    Return the absolute path to a file in the test directory.
    The name specified should be relative to test.
    """
    return join(dirname(__file__), path)


def load_bam(bam_path):
    return AlignmentFile(data_path(bam_path))


def load_vcf(vcf_path, genome=None):
    return varcode.load_vcf(data_path(vcf_path), genome=genome)
