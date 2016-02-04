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

import argparse

import varlens
import varlens.support

parser = argparse.ArgumentParser()

parser.add_argument(
    "--vcf",
    default="/Users/iskander/code/varlens/test/data/CELSR1/vcfs/vcf_1.vcf")

parser.add_argument(
    "--bam",
    default="/Users/iskander/code/varlens/test/data/CELSR1/bams/bam_1.bam")

parser.add_argument(
    "--genome",
    default="hg19")

if __name__ == "__main__":
    args = parser.parse_args()
    variants_df = varlens.variants_util.load_as_dataframe(
        args.vcf,
        genome=args.genome)
    print(variants_df)

    read_source = varlens.reads_util.load_bam(args.bam)

    variant_to_locus_dict = {}
    for variant in variants_df["variant"]:
        locus = varlens.read_evidence.pileup_collection.to_locus(variant)
        variant_to_locus_dict[variant] = locus

    variant_loci = list(variant_to_locus_dict.values())
    for read in read_source.reads(variant_loci):
        print(read.__dict__)
        print(read)
