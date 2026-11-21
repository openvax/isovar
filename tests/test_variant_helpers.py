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

from isovar.variant_helpers import (
    trim_variant,
    base0_interval_for_variant,
    interbase_range_affected_by_variant_on_transcript
)

from varcode import Variant

from .genomes_for_testing import grch38
from .common import eq_


def test_trim_variant_substitution():
    loc, ref, alt = trim_variant(Variant("chr1", 10, "C", "G", grch38))
    eq_(loc, 10)
    eq_(ref, "C")
    eq_(alt, "G")


def test_trim_variant_insertion():
    loc, ref, alt = trim_variant(Variant("chr1", 10, "C", "CG"))
    eq_(loc, 10)
    eq_(ref, "")
    eq_(alt, "G")


def test_trim_variant_deletion():
    loc, ref, alt = trim_variant(Variant("chr1", 10, "CG", "C"))
    eq_(loc, 11)
    eq_(ref, "G")
    eq_(alt, "")


def test_trim_variant_relies_on_varcode_prefix_suffix_and_insertion_normalization():
    # Varcode trims shared flanking prefixes and suffixes during Variant
    # construction and anchors insertions on the preceding base.
    prefixed_and_suffixed_sub = Variant("chr1", 10, "ACT", "AGT")
    eq_(trim_variant(prefixed_and_suffixed_sub), (11, "C", "G"))

    prefixed_and_suffixed_del = Variant("chr1", 10, "ACGT", "AT")
    eq_(trim_variant(prefixed_and_suffixed_del), (11, "CG", ""))

    prefixed_and_suffixed_ins = Variant("chr1", 10, "AT", "ACGT")
    eq_(trim_variant(prefixed_and_suffixed_ins), (10, "", "CG"))
    assert prefixed_and_suffixed_ins.is_insertion

    # Duck-typed VCF-style variant records normalize via Varcode's shared
    # flanking-string helper with identical coordinates.
    from types import SimpleNamespace
    eq_(trim_variant(SimpleNamespace(start=10, ref="CG", alt="C")), (11, "G", ""))
    eq_(trim_variant(SimpleNamespace(start=10, ref="C", alt="CG")), (10, "", "G"))
    eq_(trim_variant(SimpleNamespace(start=10, ref="ACGT", alt="AT")), (11, "CG", ""))


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
    brca2_variant_rs769125639 = Variant("13", 32315479, "T", "A", grch38)
    brca2_001 = grch38.transcripts_by_name("BRCA2-001")[0]
    interbase_range = interbase_range_affected_by_variant_on_transcript(
        variant=brca2_variant_rs769125639,
        transcript=brca2_001)
    print(interbase_range)
    eq_(interbase_range, (5, 6))


def test_interbase_range_for_brca2_utr_insertion():
    # T>TC insertion after the 6th nucleotide of BRCA2-001's 5' UTR
    brca2_insertion = Variant("13", 32315479, "T", "TC", grch38)
    brca2_001 = grch38.transcripts_by_name("BRCA2-001")[0]
    interbase_range = interbase_range_affected_by_variant_on_transcript(
        variant=brca2_insertion,
        transcript=brca2_001)
    print(interbase_range)
    eq_(interbase_range, (6, 6))


def test_interbase_range_for_brca2_utr_deletion():
    # Deletion of the 6th nucleotide of BRCA2-001's 5' UTR
    brca2_deletion = Variant("13", 32315479, "T", "", grch38)
    brca2_001 = grch38.transcripts_by_name("BRCA2-001")[0]
    interbase_range = interbase_range_affected_by_variant_on_transcript(
        variant=brca2_deletion,
        transcript=brca2_001)
    print(interbase_range)
    eq_(interbase_range, (5, 6))
