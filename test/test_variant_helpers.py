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

from isovar.variant_helpers import (
    trim_variant_fields,
    trim_variant,
    base0_interval_for_variant,
    interbase_range_affected_by_variant_on_transcript
)
from nose.tools import eq_
from varcode import Variant
from pyensembl import ensembl_grch38


def test_trim_variant_substitution():
    loc, ref, alt = trim_variant(Variant("chr1", 10, "C", "G"))
    eq_(loc, 10)
    eq_(ref, "C")
    eq_(alt, "G")

def test_trim_variant_fields_substitution():
    loc, ref, alt = trim_variant_fields(10, "C", "G")
    eq_(loc, 10)
    eq_(ref, "C")
    eq_(alt, "G")

def test_trim_variant_insertion():
    loc, ref, alt = trim_variant(Variant("chr1", 10, "C", "CG"))
    eq_(loc, 10)
    eq_(ref, "")
    eq_(alt, "G")


def test_trim_variant_fields_insertion():
    loc, ref, alt = trim_variant_fields(10, "C", "CG")
    eq_(loc, 10)
    eq_(ref, "")
    eq_(alt, "G")


def test_trim_variant_deletion():
    loc, ref, alt = trim_variant(Variant("chr1", 10, "CG", "C"))
    eq_(loc, 11)
    eq_(ref, "G")
    eq_(alt, "")

def test_trim_variant_fields_deletion():
    loc, ref, alt = trim_variant_fields(10, "CG", "C")
    eq_(loc, 11)
    eq_(ref, "G")
    eq_(alt, "")

def test_base0_interval_for_variant_substitution():
    (start, end) = base0_interval_for_variant(Variant("chr1", 10, "C", "G"))
    eq_(start, 9)
    eq_(end, 10)

def test_base0_interval_for_variant_insertion():
    (start, end) = base0_interval_for_variant(Variant("chr1", 10, "C", "CG"))
    eq_(start, 10)
    eq_(end, 10)

def test_base0_interval_for_variant_deletion():
    (start, end) = base0_interval_for_variant(Variant("chr1", 10, "CG", "C"))
    eq_(start, 10)
    eq_(end, 11)

def test_interbase_range_for_brca2_utr_substitution():
    # rs769125639 is a simple T>A substitution in the 6th nucleotide of
    # BRCA2-001's 5' UTR
    brca2_variant_rs769125639 = Variant(
        "13", 32315479, "T", "A", ensembl_grch38)
    brca2_001 = ensembl_grch38.transcripts_by_name("BRCA2-001")[0]
    interbase_range = interbase_range_affected_by_variant_on_transcript(
        variant=brca2_variant_rs769125639,
        transcript=brca2_001)
    print(interbase_range)
    eq_(interbase_range, (5, 6))

def test_interbase_range_for_brca2_utr_insertion():
    # T>TC insertion after the 6th nucleotide of BRCA2-001's 5' UTR
    brca2_insertion = Variant("13", 32315479, "T", "TC", ensembl_grch38)
    brca2_001 = ensembl_grch38.transcripts_by_name("BRCA2-001")[0]
    interbase_range = interbase_range_affected_by_variant_on_transcript(
        variant=brca2_insertion,
        transcript=brca2_001)
    print(interbase_range)
    eq_(interbase_range, (6, 6))

def test_interbase_range_for_brca2_utr_deletion():
    # Deletion of the 6th nucleotide of BRCA2-001's 5' UTR
    brca2_deletion = Variant("13", 32315479, "T", "", ensembl_grch38)
    brca2_001 = ensembl_grch38.transcripts_by_name("BRCA2-001")[0]
    interbase_range = interbase_range_affected_by_variant_on_transcript(
        variant=brca2_deletion,
        transcript=brca2_001)
    print(interbase_range)
    eq_(interbase_range, (5, 6))
