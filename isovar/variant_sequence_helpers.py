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
Helper functions for constructing and filtering VariantSequence objects
from reads overlapping a variant locus.
"""

from __future__ import print_function, division, absolute_import

from .allele_read_helpers import group_unique_sequences
from .assembly import collapse_substrings
from .logging import get_logger
from .variant_sequence import VariantSequence

logger = get_logger(__name__)


def initial_variant_sequences_from_reads(
        variant_reads,
        max_nucleotides_before_variant=None,
        max_nucleotides_after_variant=None):
    """
    Get all unique sequences from reads spanning a variant locus. This will
    include partial sequences due to reads starting in the middle of the
    sequence around around a variant.

    Parameters
    ----------
    variant_reads : list of AlleleRead objects

    max_nucleotides_before_variant : int or None

    max_nucleotides_after_variant : int or None

    Returns
    -------
    list of VariantSequence
    """
    unique_sequence_groups = group_unique_sequences(
        variant_reads,
        max_prefix_size=max_nucleotides_before_variant,
        max_suffix_size=max_nucleotides_after_variant)

    return [
        VariantSequence(
            prefix=prefix,
            alt=alt,
            suffix=suffix,
            reads=reads)
        for ((prefix, alt, suffix), reads)
        in unique_sequence_groups.items()
    ]


def filter_variant_sequences_by_read_support(
        variant_sequences,
        min_variant_sequence_coverage):
    """
    Filter VariantSequences to only keep those with at least the desired
    level of coverage.

    Parameters
    ----------
    variant_sequences : list of VariantSequence

    min_variant_sequence_coverage : int
        Minimum number of reads which must cover each
        base of a VariantSequence

    Returns
    -------
    list of VariantSequence
    """
    n_total = len(variant_sequences)
    variant_sequences = [
        s
        for s in variant_sequences
        if s.min_coverage() >= min_variant_sequence_coverage
    ]
    n_dropped = n_total - len(variant_sequences)
    if n_dropped > 0:
        logger.info(
            "Dropped %d/%d variant sequences less than %d supporting reads",
            n_dropped,
            n_total,
            min_variant_sequence_coverage)
    return variant_sequences


def filter_variant_sequences_by_length(
        variant_sequences,
        preferred_sequence_length):
    """
    Parameters
    ----------
    variant_sequences : list of VariantSequence

    preferred_sequence_length : int
        If we get some sequences which are at least this long and others
        which are shorter, then drop the shorter ones.

    Returns
    -------
    list of VariantSequence
    """
    n_total = len(variant_sequences)
    if n_total == 0:
        return []
    # since we might have gotten some shorter fragments,
    # keep only the longest spanning sequence
    max_observed_sequence_length = max(len(s) for s in variant_sequences)

    # if we get back a sequence that's longer than the preferred length
    # then that doesn't mean we should necessarily drop the other sequences
    min_required_sequence_length = min(
        max_observed_sequence_length,
        preferred_sequence_length)

    variant_sequences = [
        s for s in variant_sequences
        if len(s.sequence) >= min_required_sequence_length
    ]
    n_dropped = n_total - len(variant_sequences)
    if n_dropped > 0:
        logger.info(
            "Dropped %d/%d variant sequences shorter than %d",
            n_dropped,
            n_total,
            min_required_sequence_length)
    return variant_sequences


def trim_variant_sequences(variant_sequences, min_variant_sequence_coverage):
    """
    Trim VariantSequences to desired coverage and then combine any
    subsequences which get generated.

    Parameters
    ----------
    variant_sequences : list of VariantSequence

    min_variant_sequence_coverage : int

    Returns list of VariantSequence
    """
    n_total = len(variant_sequences)
    trimmed_variant_sequences = [
        variant_sequence.trim_by_coverage(min_variant_sequence_coverage)
        for variant_sequence in variant_sequences
    ]
    collapsed_variant_sequences = collapse_substrings(trimmed_variant_sequences)
    n_after_trimming = len(collapsed_variant_sequences)
    logger.info(
        "Kept %d/%d variant sequences after read coverage trimming to >=%dx",
        n_after_trimming,
        n_total,
        min_variant_sequence_coverage)
    return collapsed_variant_sequences


def filter_variant_sequences(
        variant_sequences,
        preferred_sequence_length,
        min_variant_sequence_coverage):
    """
    Drop variant sequences which are shorter than request or don't have
    enough supporting reads.

    Parameters
    ----------
    variant_sequences : list of VariantSequence

    preferred_sequence_length : int

    min_variant_sequence_coverage : int

    Returns list of VariantSequence
    """
    variant_sequences = trim_variant_sequences(
        variant_sequences, min_variant_sequence_coverage)

    return filter_variant_sequences_by_length(
        variant_sequences=variant_sequences,
        preferred_sequence_length=preferred_sequence_length)

