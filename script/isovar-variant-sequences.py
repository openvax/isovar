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

from isovar.args import add_somatic_vcf_args
from isovar.variant_sequence import variant_sequences_dataframe
from isovar.default_parameters import (
    VARIANT_CDNA_SEQUENCE_LENGTH,
    MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE
)

parser = argparse.ArgumentParser()
parser = add_somatic_vcf_args(parser)

parser.add_argument(
    "--bam",
    required=True,
    help="Path to BAM or SAM file containing RNAseq reads")

parser.add_argument(
    "--min-reads",
    type=int,
    default=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
    help="Minimum number of reads supporting a variant sequence")

parser.add_argument(
    "--sequence-length",
    default=VARIANT_CDNA_SEQUENCE_LENGTH,
    type=int)

parser.add_argument(
    "--output",
    default="isovar-variant-sequences-results.csv",
    help="Name of CSV file which contains predicted sequences")

if __name__ == "__main__":
    args = parser.parse_args()

    print(args)

    variants = varcode.load_vcf(
        args.vcf,
        genome=args.genome)

    samfile = AlignmentFile(args.bam)

    df = variant_sequences_dataframe(
        variants=variants,
        samfile=samfile,
        sequence_length=args.sequence_length,
        min_reads=args.min_reads)

    print(df)
    df.to_csv(args.output)
