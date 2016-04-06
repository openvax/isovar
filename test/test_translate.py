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


from isovar.translate import translate_variants

from testing_helpers import load_bam, load_vcf


VCF = "data/CELSR1/vcfs/vcf_1.vcf"

BAM = "data/CELSR1/bams/bam_1.bam"

GENOME = "hg19"


def test_translate_variant_collection():
    variants = load_vcf(VCF, genome=GENOME)
    samfile = load_bam(BAM)
    result = translate_variants(variants, samfile)
    print(result)
    assert len(result) > 0, result

if __name__ == "__main__":
    test_translate_variant_collection()
