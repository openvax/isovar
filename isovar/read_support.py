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

import varlens
import varlens.support

def load_variants(vcf_path, rna_bam_path, genome=None):
    variants_df = varlens.variants_util.load_as_dataframe(
        vcf_path,
        genome=genome)
    variant_to_locus_dict = {}
    for variant in variants_df["variant"]:
        locus = varlens.read_evidence.pileup_collection.to_locus(variant)
        variant_to_locus_dict[variant] = locus

    variant_loci = list(variant_to_locus_dict.values())
