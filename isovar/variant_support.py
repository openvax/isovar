# Copyright (c) 2018. Mount Sinai School of Medicine
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
Collect AlleleReads overlapping each variant grouped by whether
they support the reference, somatic allele, or some other allele.

Allows for filtering variants by RNA evidence.
"""

from .allele_read_helpers import group_reads_by_allele
from .variant_helpers import trim_variant

from collections import namedtuple

VariantSupport = namedtuple("VariantSupport", [
    "trimemd_start",
    "trimmed_ref",
    "trimmed_alt",
    "ref_reads",
    "alt_reads",
    "other_reads"
])

VariantSupportStats = namedtuple("VariantSupportStats", [
    "n_reads",
    "n_fragments",
    "n_ref_reads",
    "n_ref_fragments",
    "n_alt_reads",
    "n_alt_fragments",
    "n_other_reads",
    "n_other_fragments",
])


def split_reads_into_ref_alt_other(ref, alt, overlapping_reads):
    """
    Returns three lists of AlleleRead objects
        - reads which support the reference allele
        - reads which support the variant's alt allele
        - reads which support other alleles
    """
    # convert to list in case it's a generator since
    # we want to traverse the sequence repeatedly
    overlapping_reads = list(overlapping_reads)

    reads_grouped_by_allele = group_reads_by_allele(overlapping_reads)
    ref_reads = reads_grouped_by_allele.get(ref, [])
    alt_reads = reads_grouped_by_allele.get(alt, [])
    other_reads = []
    for allele, allele_reads in reads_grouped_by_allele.items():
        if allele in {ref, alt}:
            continue
        other_reads.extend(allele_reads)
    return ref_reads, alt_reads, other_reads


def gather_variant_support(variant, overlapping_reads):
    """
    Returns VariantSupport object for given variant and the reads
    which align to its location.
    """
    start, ref, alt = trim_variant(variant)
    ref_reads, alt_reads, other_reads = split_reads_into_ref_alt_other(
        ref=ref, alt=alt, overlapping_reads=overlapping_reads)
    return VariantSupport(
        trimmed_start=start,
        trimmed_ref=ref,
        trimmed_alt=alt,
        ref_reads=ref_reads,
        alt_reads=alt_reads,
        other_reads=other_reads)


def compute_variant_support_stats(variant_support):
    n_ref_reads = len(variant_support.ref_reads)
    n_alt_reads = len(variant_support.alt_reads)
    n_other_reads = len(variant_support.other_reads)
    n_total_reads = n_ref_reads + n_alt_reads + n_other_reads

    n_ref_fragments = len({r.name for r in variant_support.ref_reads})
    n_alt_fragments = len({r.name for r in variant_support.alt_reads})
    n_other_fragments = len({r.name for r in variant_support.other_reads})
    n_total_fragments = n_ref_fragments + n_alt_fragments + n_other_fragments

    return VariantSupportStats(
        n_reads=n_total_reads,
        n_ref_reads=n_ref_reads,
        n_alt_reads=n_alt_reads,
        n_other_reads=n_other_reads,
        n_fragments=n_total_fragments,
        n_ref_fragments=n_ref_fragments,
        n_alt_fragments=n_alt_fragments,
        n_other_fragments=n_other_fragments)


def variant_support_generator(variant_and_overlapping_read_generator):
    for variant, overlapping_reads in variant_and_overlapping_read_generator:
        pass
"""
def filter_reads():
    if n_alt_reads < self.min_alt_rna_reads:
        logger.info(
            "Skipping %s because only %d alt RNA reads (min=%d)",
            variant,
            n_alt_reads,
            min_alt_rna_reads)
        return []
    elif n_alt_fragments < self.min_alt_rna_fragments:
        pass
"""

