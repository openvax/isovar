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

from collections import OrderedDict
import argparse

import varcode
from pysam import AlignmentFile
import pandas as pd

from isovar import translate_variant_collection

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
    variants = varcode.load_vcf(args.vcf, genome=args.genome)

    samfile = AlignmentFile(args.bam)
    variant_to_proteins = translate_variant_collection(
        variants,
        samfile,
        protein_fragment_length=args.protein_fragment_length,
        min_reads_supporting_rna_sequence=args.min_reads,
        min_transcript_prefix_length=args.min_transcript_prefix_length,
        max_transcript_mismatches=args.max_reference_transcript_mismatches,
        max_sequences_per_variant=args.max_sequences_per_variant)

    # construct a dictionary incrementally which we'll turn into a
    # DataFrame
    column_dict = OrderedDict([
        ("chr", []),
        ("base1_start_pos", []),
        ("base1_end_pos", []),
        ("ref", []),
        ("alt", []),
        ("variant_protein_sequence", []),
        ("variant_protein_sequence_length", []),
        ("reference_transcript_ids", []),
        ("reference_transcript_names", []),
        ("reference_protein_sequences", []),
        ("cdna_sequences", []),
        ("cdna_sequence_length", []),
        ("number_supporting_reads", [])
    ])
    for (variant, protein_fragments) in variant_to_proteins.items():
        print(variant)
        for protein_fragment in protein_fragments:
            column_dict["chr"].append(variant.contig)
            column_dict["base1_start_pos"].append(variant.start)
            column_dict["base1_end_pos"].append(variant.end)
            column_dict["ref"].append(variant.ref)
            column_dict["alt"].append(variant.alt)
            column_dict["variant_protein_sequence"].append(
                protein_fragment.variant_protein_sequence)
            column_dict["variant_protein_sequence_length"].append(
                len(protein_fragment.variant_protein_sequence))
            column_dict["reference_transcript_ids"].append(
                ";".join(protein_fragment.reference_transcript_ids))
            column_dict["reference_transcript_names"].append(
                ";".join(protein_fragment.reference_transcript_names))
            column_dict["reference_protein_sequences"].append(
                ";".join(protein_fragment.reference_protein_sequences))
            column_dict["cdna_sequences"].append(
                ";".join([
                    "%s_%s_%s" % (prefix, alt, suffix)
                    for (prefix, alt, suffix)
                    in protein_fragment.cdna_sequence_tuples
                ]))
            column_dict["cdna_sequence_length"].append(
                sum(len(prefix) + len(alt) + len(suffix)
                    for (prefix, alt, suffix)
                    in protein_fragment.cdna_sequence_tuples) / len(
                    protein_fragment.cdna_sequence_tuples))
            column_dict["number_supporting_reads"].append(
                protein_fragment.number_supporting_reads)
    df = pd.DataFrame(column_dict)
    print(df)
    df.to_csv(args.output)
