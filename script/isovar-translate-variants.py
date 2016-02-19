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
import skbio
from pysam import AlignmentFile

from isovar import gather_variant_reads, sequence_counts

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
    "--min-count", type=int, default=3)

parser.add_argument(
    "--context-size",
    default=45,
    type=int)


if __name__ == "__main__":
    args = parser.parse_args()
    variants = varcode.load_vcf(args.vcf, genome=args.genome)

    samfile = AlignmentFile(args.bam)

    for variant in variants:
        variant_reads = gather_variant_reads(
            samfile=samfile,
            chromosome="chr" + variant.contig,
            base1_location=variant.start,
            ref=variant.ref,
            alt=variant.alt)
        if len(variant_reads) < args.min_count:
            continue
        sequence_count_info = sequence_counts(
            variant_reads, context_size=args.context_size)
        for ((prefix, suffix), count) in sorted(
                sequence_count_info.full_read_counts.items(),
                key=lambda x: -x[1]):
            if count < args.min_count:
                break

            variant = sequence_count_info.variant_nucleotides
            print("\t%s_%s_%s: %d" % (
                prefix,
                variant,
                suffix,
                count))

            # translate in three reading frames:
            seq = "%s%s%s" % (prefix, variant, suffix)
            for offset in range(3):
                dna = skbio.DNA(seq[offset:])
                print("\t\tframe=%d: %s" % (offset, dna.translate()))
