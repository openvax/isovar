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
Helper functions which it easier to write shell scripts
"""

from __future__ import print_function, division, absolute_import
import varcode
from pysam import AlignmentFile

from .protein_sequence import variants_to_protein_sequences_dataframe

def variants_to_protein_sequences_dataframe_from_args(args):
    variants = varcode.load_vcf(args.vcf, genome=args.genome)
    samfile = AlignmentFile(args.bam)
    return variants_to_protein_sequences_dataframe(
        variants=variants,
        samfile=samfile,
        protein_sequence_length=args.protein_sequence_length,
        min_reads_supporting_rna_sequence=args.min_reads,
        min_transcript_prefix_length=args.min_transcript_prefix_length,
        max_transcript_mismatches=args.max_reference_transcript_mismatches,
        max_protein_sequences_per_variant=args.max_protein_sequences_per_variant,
        min_mapping_quality=args.min_mapping_quality)
