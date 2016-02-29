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

from isovar import variant_protein_fragments_dataframe

parser = argparse.ArgumentParser()

parser.add_argument(
    "--vcf",
    default="/Users/iskander/code/varlens/test/data/CELSR1/vcfs/vcf_1.vcf")

parser.add_argument(
    "--bam",
    default="/Users/iskander/code/varlens/test/data/CELSR1/bams/bam_9.bam")

parser.add_argument(
    "--genome",
    default=None)

parser.add_argument(
    "--min-reads",
    type=int,
    default=3)

parser.add_argument(
    "--protein-fragment-length",
    default=30,
    type=int)

parser.add_argument(
    "--max-sequences-per-variant",
    type=int,
    default=5)

parser.add_argument(
    "--max-reference-transcript-mismatches",
    type=int,
    default=2)

parser.add_argument(
    "--min-transcript-prefix-length",
    type=int,
    default=15,
    help=(
        "Number of nucleotides before the variant we try to match against "
        "a reference transcript. Values greater than zero exclude variants "
        "near the start codon of transcripts without 5' UTRs."))

parser.add_argument(
    "--output",
    default="isovar-results.csv",
    help="Name of CSV file which contains predicted sequences")

if __name__ == "__main__":
    args = parser.parse_args()
    variants = varcode.load_vcf(
        args.vcf,
        genome=args.genome)

    samfile = AlignmentFile(args.bam)
    df = variant_protein_fragments_dataframe(
        variants=variants,
        samfile=samfile,
        protein_fragment_length=args.protein_fragment_length,
        min_reads_supporting_rna_sequence=args.min_reads,
        min_transcript_prefix_length=args.min_transcript_prefix_length,
        max_transcript_mismatches=args.max_reference_transcript_mismatches,
        max_sequences_per_variant=args.max_sequences_per_variant)
    print(df)
    df.to_csv(args.output)
