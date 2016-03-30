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

import varcode
from skbio import DNA
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
    default="hg19")

parser.add_argument(
    "--min-weight", type=float, default=3.0)

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
        if len(variant_reads) == 0:
            continue
        print("VARIANT: %s" % variant)
        result = sequence_counts(variant_reads, context_size=args.context_size)
        for ((prefix, suffix), weight) in sorted(
                result.sequence_weights.items(),
                key=lambda x: -x[1]):
            if weight < args.min_weight:
                continue
            variant_seq = result.variant_nucleotides
            # translate in three reading frames:
            print("%s_%s_%s: %f" % (
                prefix,
                variant_seq,
                suffix,
                weight))
            for transcript in variant.transcripts:
                if not transcript.contains(
                        variant.contig,
                        variant.start,
                        variant.end):
                    print("Skipping %s" % transcript)
                    continue
                variant_transcript_idx = transcript.spliced_offset(
                    variant.start)
                start_idx = variant_transcript_idx - args.context_size
                transcript_sequence_before_variant = transcript.sequence[
                    start_idx:variant_transcript_idx]

                if transcript.strand == "-":
                    cdna_prefix = str(DNA(suffix).reverse_complement())
                    cdna_suffix = str(DNA(prefix).reverse_complement())
                    cdna_variant = str(
                        DNA(variant_seq).reverse_complement())
                else:
                    cdna_prefix = prefix
                    cdna_suffix = suffix
                    cdna_variant = variant_seq

                length_difference = (
                    len(transcript_sequence_before_variant) - len(cdna_prefix))
                assert length_difference >= 0
                aligned_transcript_seq = transcript_sequence_before_variant[
                    length_difference:]
                n_mismatch = sum(xi != yi for (xi, yi) in zip(
                    aligned_transcript_seq, cdna_prefix))
                if transcript_sequence_before_variant.endswith(cdna_prefix):
                    print("Match! %s" % transcript)
                else:
                    print("%s: %d mismatches" % (transcript, n_mismatch))
                    print(aligned_transcript_seq)
                    print(cdna_prefix)
                    print("".join([
                        "_" if xi == yi else "!"
                        for (xi, yi) in
                        zip(aligned_transcript_seq, cdna_prefix)
                    ]))
