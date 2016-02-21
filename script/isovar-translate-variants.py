#!/usr/bin/env python

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

import argparse

import varcode
from pysam import AlignmentFile

from isovar import translate_variant_collection

parser = argparse.ArgumentParser()

parser.add_argument(
    "--vcf",
    default="/Users/iskander/code/varlens/test/data/CELSR1/vcfs/vcf_1.vcf")

parser.add_argument(
    "--bam",
    default="/Users/iskander/code/varlens/test/data/CELSR1/bams/bam_1.bam")

parser.add_argument(
    "--genome",
    default=None)

parser.add_argument(
    "--min-read-count",
    type=int,
    default=3)

parser.add_argument(
    "--sequence-length",
    default=105,
    type=int)

parser.add_argument(
    "--max-sequences-per-variant",
    type=int,
    default=5)

if __name__ == "__main__":
    args = parser.parse_args()
    variants = varcode.load_vcf(args.vcf, genome=args.genome)

    samfile = AlignmentFile(args.bam)
    print(variants)
    print(translate_variant_collection(variants, samfile))
