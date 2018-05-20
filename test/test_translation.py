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


from isovar.translation import translate_variants
from isovar.variant_reads import reads_supporting_variants

from nose.tools import eq_

from testing_helpers import load_bam, load_vcf


def test_translate_variant_collection():
    variants = load_vcf("data/b16.f10/b16.vcf")
    samfile = load_bam("data/b16.f10/b16.combined.sorted.bam")

    result = list(translate_variants(reads_supporting_variants(variants, samfile)))
    eq_(
        len(result),
        4,
        "Expected %d translated variants but got %d: %s" % (
            len(variants),
            len(result),
            result))
