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

"""
Translates cDNA sequences associated with variants into proteins using
all 3 possible open reading frames.
"""

from __future__ import print_function, division, absolute_import

import argparse

import varcode
import skbio
from pysam import AlignmentFile

from isovar.variant_reads import gather_variant_reads, sequence_counts
from isovar.default_parameters import (
    MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
    CDNA_CONTEXT_SIZE
)

parser = argparse.ArgumentParser()

parser.add_argument(
    "--vcf",
    default="../test/data/CELSR1/vcfs/vcf_1.vcf")

parser.add_argument(
    "--bam",
    default="../test/data/CELSR1/bams/bam_1.bam")

parser.add_argument(
    "--genome",
    default=None)

parser.add_argument(
    "--min-reads",
    type=int,
    default=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE)

parser.add_argument(
    "--context-size",
    default=CDNA_CONTEXT_SIZE,
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
        if len(variant_reads) == 0:
            continue
        result = sequence_counts(variant_reads, context_size=args.context_size)
        for ((prefix, suffix), count) in sorted(
                result.full_read_counts.items(),
                key=lambda x: -x[1]):
            if count < args.min_reads:
                continue

            variant = result.alt
            print("\t%s|%s|%s: %d" % (
                prefix,
                variant,
                suffix,
                count))

            # translate in three reading frames:
            seq = "%s%s%s" % (prefix, variant, suffix)
            for offset in range(3):
                dna = skbio.DNA(seq[offset:])
                print("\t\tframe=%d: %s" % (offset, dna.translate()))
