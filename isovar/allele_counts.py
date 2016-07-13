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
from collections import namedtuple

from .variant_helpers import trim_variant
from .dataframe_builder import DataFrameBuilder
from .read_helpers import group_reads_by_allele

AlleleCount = namedtuple(
    "AlleleCount", [
        "n_ref",
        "n_alt",
        "n_other",
    ])

def count_alleles_at_variant_locus(variant, allele_reads):
    allele_reads = list(allele_reads)
    n_total = len(allele_reads)
    allele_to_reads_dict = group_reads_by_allele(allele_reads)
    _, ref, alt = trim_variant(variant)
    n_ref = len(allele_to_reads_dict[ref])
    n_alt = len(allele_to_reads_dict[alt])
    n_other = n_total - (n_ref + n_alt)
    return AlleleCount(n_ref=n_ref, n_alt=n_alt, n_other=n_other)

def allele_counts_dataframe(variant_and_allele_reads_generator):
    """
    Creates a DataFrame containing number of reads supporting the
    ref vs. alt alleles for each variant.
    """
    df_builder = DataFrameBuilder(
        AlleleCount,
        extra_column_fns={
            "gene": lambda variant, _: ";".join(variant.gene_names),
        })
    for variant, allele_reads in variant_and_allele_reads_generator:
        counts = count_alleles_at_variant_locus(variant, allele_reads)
        df_builder.add(variant, counts)
    return df_builder.to_dataframe()
