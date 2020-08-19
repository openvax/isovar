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
Functions for filtering, grouping, and summarizing collections of
AlleleRead objects.
"""

from collections import defaultdict

from .common import groupby
from .logging import get_logger
from .allele_read import AlleleRead

logger = get_logger(__name__)


def group_reads_by_allele(allele_reads):
    """
    Returns dictionary mapping each allele's nucleotide sequence to a list of
    supporting AlleleRead objects.
    """
    return groupby(allele_reads, lambda read: read.allele)


def get_single_allele_from_reads(allele_reads):
    """
    Given a sequence of AlleleRead objects, which are expected to all have
    the same allele, return that allele.
    """
    allele_reads = list(allele_reads)

    if len(allele_reads) == 0:
        raise ValueError("Expected non-empty list of AlleleRead objects")

    seq = allele_reads[0].allele
    if any(read.allele != seq for read in allele_reads):
        raise ValueError("Expected all AlleleRead objects to have same allele '%s', got %s" % (
            seq, allele_reads))
    return seq


def group_unique_sequences(
        allele_reads,
        max_prefix_size=None,
        max_suffix_size=None):
    """
    Given a list of AlleleRead objects, extracts all unique
    (prefix, allele, suffix) sequences and associate each with a list
    of reads that contained that sequence.
    """
    groups = defaultdict(set)
    for r in allele_reads:
        prefix = r.prefix
        allele = r.allele
        suffix = r.suffix
        if max_prefix_size and len(prefix) > max_prefix_size:
            prefix = prefix[-max_prefix_size:]
        if max_suffix_size and len(suffix) > max_suffix_size:
            suffix = suffix[:max_suffix_size]
        key = (prefix, allele, suffix)
        groups[key].add(r)
    return groups


def allele_reads_from_locus_reads(locus_reads):
    """
    Attempt to convert each LocusRead object to an AlleleRead and return
    the successfully converted objects.

    Parameters
    ----------
    locus_reads : list of LocusRead

    Returns list of AlleleRead
    -------

    """
    allele_reads = []
    for locus_read in locus_reads:
        allele_read = AlleleRead.from_locus_read(locus_read)
        if allele_read is None:
            continue
        else:
            allele_reads.append(allele_read)
    return allele_reads


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