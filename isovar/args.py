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
Helper functions for building ArgumentParser instances. Keeps us from
having to copy the same arguments across many scripts.
"""

from __future__ import print_function, division, absolute_import

from .default_parameters import (
    MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
    MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
    MIN_TRANSCRIPT_PREFIX_LENGTH,
    MIN_READ_MAPPING_QUALITY,
    PROTEIN_SEQUENCE_LENGTH,
    MAX_PROTEIN_SEQUENCES_PER_VARIANT,
)




def add_rna_consensus_sequence_args(
        parser,
        min_rna_reads_arg="--min-reads",
        min_rna_reads_default=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE):
    rna_sequence_group = parser.add_argument_group("Consensus expressed coding sequence")
    rna_sequence_group.add_argument(
        min_rna_reads_arg,
        type=int,
        default=min_rna_reads_default,
        help="Minimum number of reads supporting a variant sequence")
    return parser

def add_protein_sequence_args(
        parser,
        protein_sequence_length_arg="--protein-sequence-length",
        protein_sequence_length_default=PROTEIN_SEQUENCE_LENGTH,
        max_protein_sequences_per_variant_arg="--max-protein-sequences-per-variant",
        max_protein_sequences_per_variant_default=MAX_PROTEIN_SEQUENCES_PER_VARIANT):
    """
    Extends an ArgumentParser instance with the following args:
        --protein-sequence-length
        --max-protein-sequences-per-variant
    """
    protein_sequence_group = parser.add_argument_group("Protein Sequence")
    protein_sequence_group.add_argument(
        protein_sequence_length_arg,
        default=protein_sequence_length_default,
        type=int)
    protein_sequence_group.add_argument(
        max_protein_sequences_per_variant_arg,
        type=int,
        default=max_protein_sequences_per_variant_default)
    return parser
