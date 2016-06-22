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
Summarizes which alleles are found at a locus overlapping a variant.
"""

from __future__ import print_function, division, absolute_import
from collections import namedtuple, defaultdict

from .default_parameters import (
    MIN_READ_MAPPING_QUALITY,
    USE_SECONDARY_ALIGNMENTS,
    USE_DUPLICATE_READS,
)
from .variant_helpers import trim_variant
from .dataframe_builder import DataFrameBuilder
from .allele_read import allele_reads_for_variants

AlleleCount = namedtuple(
    "AlleleCount", [
        "n_ref",
        "n_alt",
        "n_other",
    ])

def group_reads_by_allele(allele_reads):
    """
    Returns dictionary mapping each allele's nucleotide sequence to a list of
    supporting AlleleRead objects.
    """
    allele_to_reads_dict = defaultdict(list)
    for allele_read in allele_reads:
        allele_to_reads_dict[allele_read.allele].append(allele_read)
    return allele_to_reads_dict

def count_alleles_at_variant_locus(variant, allele_reads):
    allele_reads = list(allele_reads)
    n_total = len(allele_reads)
    allele_to_reads_dict = group_reads_by_allele(allele_reads)
    _, ref, alt = trim_variant(variant)
    n_ref = len(allele_to_reads_dict[ref])
    n_alt = len(allele_to_reads_dict[alt])
    n_other = n_total - (n_ref + n_alt)
    return AlleleCount(n_ref=n_ref, n_alt=n_alt, n_other=n_other)

def allele_counts_dataframe(
        variants,
        samfile,
        use_duplicate_reads=USE_DUPLICATE_READS,
        use_secondary_alignments=USE_SECONDARY_ALIGNMENTS,
        min_mapping_quality=MIN_READ_MAPPING_QUALITY):
    """
    Creates a DataFrame containing number of reads supporting the
    ref vs. alt alleles for each variant.

    Parameters
    ----------
    variants : varcode.VariantCollection

    samfile : pysam.AlignmentFile

    use_duplicate_reads : bool
        Should we use reads that have been marked as PCR duplicates

    use_secondary_alignments : bool
        Should we use reads at locations other than their best alignment

    min_mapping_quality : int
        Drop reads below this mapping quality
    """
    df_builder = DataFrameBuilder(AlleleCount)
    for variant, allele_reads in allele_reads_for_variants(
            variants,
            samfile,
            use_duplicate_reads=use_duplicate_reads,
            use_secondary_alignments=use_secondary_alignments,
            min_mapping_quality=min_mapping_quality):
        counts = count_alleles_at_variant_locus(variant, allele_reads)
        df_builder.add(variant, counts)
    return df_builder.to_dataframe()
