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
import numpy as np
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

    for variant in variants:
        print(variant)
        variant_reads = gather_variant_reads(
            samfile=samfile,
            chromosome="chr" + variant.contig,
            base1_location=variant.start,
            ref=variant.ref,
            alt=variant.alt)
        if len(variant_reads) < args.min_read_count:
            continue

        # the number of context nucleotides on either side of the variant
        # is half the desired length (minus the number of variant nucleotides)
        context_size = int(
            np.ceil((args.sequence_length - len(variant.alt)) / 2.0))
        sequence_count_info = sequence_counts(
            variant_reads,
            context_size=context_size)
        for i, ((prefix, suffix), count) in enumerate(sorted(
                sequence_count_info.full_read_counts.items(),
                key=lambda x: -x[1])):
            if i >= args.max_sequences_per_variant:
                break

            if count < args.min_read_count:
                break

            variant_seq = sequence_count_info.variant_nucleotides

            print("\t%s_%s_%s: %d" % (
                prefix,
                variant_seq,
                suffix,
                count))

            # translate in three reading frames:
            seq = "%s%s%s" % (prefix, variant_seq, suffix)
            for offset in range(3):
                dna = skbio.DNA(seq[offset:])
                print("\t\tframe=%d: %s" % (offset, dna.translate()))
