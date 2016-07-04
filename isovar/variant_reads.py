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
Functions for getting AlleleReads which support variant alleles.
"""

from .variant_helpers import trim_variant
from .allele_reads import reads_overlapping_variant, reads_overlapping_variants

def filter_non_alt_reads_for_variant(variant, allele_reads):
    _, _, alt = trim_variant(variant)
    return [read for read in allele_reads if read.allele == alt]

def filter_non_alt_reads_for_variants(variants_and_allele_reads_sequence):
    """
    Given a sequence of variants paired with all of their overlapping reads,
    yields a sequence of variants paired only with reads which contain their
    mutated nucleotide sequence.
    """
    for variant, allele_reads in variants_and_allele_reads_sequence:
        yield variant, filter_non_alt_reads_for_variant(variant, allele_reads)

def reads_supporting_variant(variant, samfile, **kwargs):
    allele_reads = reads_overlapping_variant(
        variant=variant,
        samfile=samfile,
        **kwargs)
    return filter_non_alt_reads_for_variant(
        variant=variant,
        allele_reads=allele_reads)

def reads_supporting_variants(variants, samfile, **kwargs):
    """
    Given a SAM/BAM file and a collection of variants, generates a sequence
    of variants paired with reads which support each variant.
    """
    for variant, allele_reads in reads_overlapping_variants(
            variants=variants,
            samfile=samfile,
            **kwargs):
        yield variant, filter_non_alt_reads_for_variant(variant, allele_reads)
