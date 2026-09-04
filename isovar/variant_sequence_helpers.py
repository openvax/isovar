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

from .allele_read_helpers import group_unique_sequences
from .assembly import merge_identical_sequences
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
    Retain nonempty variant sequences for downstream translation.

    This function historically dropped every sequence shorter than the
    longest observed sequence (capped at ``preferred_sequence_length``). That
    decision cannot be made safely before translation: a longer sequence may
    fail to match every reference context, raw read counts are only the second
    part of the downstream support ordering, and multiple cDNA sequences may
    be grouped into one protein sequence whose combined support beats every
    individual candidate.

    Keep the legacy function name and ``preferred_sequence_length`` argument
    for API compatibility, but defer all candidate ranking until after
    reference matching and translation. Remove zero-length sequences here as
    a defense in depth; ``trim_variant_sequences`` normally removes them where
    they are produced.

    Parameters
    ----------
    variant_sequences : list of VariantSequence

    preferred_sequence_length : int
        Retained for API compatibility. Initial sequence construction already
        uses this value to bound assembled sequence length.

    Returns
    -------
    list of VariantSequence
    """
    n_total = len(variant_sequences)
    nonempty_variant_sequences = [
        s for s in variant_sequences
        if len(s) > 0
    ]
    n_dropped = n_total - len(nonempty_variant_sequences)
    if n_dropped > 0:
        logger.info(
            "Dropped %d/%d empty variant sequences before translation",
            n_dropped,
            n_total)
    return nonempty_variant_sequences


def trim_variant_sequences(variant_sequences, min_variant_sequence_coverage):
    """
    Trim VariantSequences to desired coverage and combine exact duplicates.

    Contained but non-identical candidates remain separate because reference
    compatibility is not known at this stage.

    Parameters
    ----------
    variant_sequences : list of VariantSequence

    min_variant_sequence_coverage : int

    Returns list of VariantSequence
    """
    n_total = len(variant_sequences)
    # when no base of a sequence has sufficient coverage, trim_by_coverage
    # returns an empty sequence which still carries every one of the original
    # reads. Those sequences are useless downstream and their read counts are
    # actively misleading when comparing candidates by support, so drop them
    # here rather than relying on a later validity check to sweep them up.
    trimmed_variant_sequences = [
        trimmed
        for trimmed in (
            variant_sequence.trim_by_coverage(min_variant_sequence_coverage)
            for variant_sequence in variant_sequences
        )
        if len(trimmed) > 0
    ]
    trimmed_variant_sequences = merge_identical_sequences(
        trimmed_variant_sequences)
    n_after_trimming = len(trimmed_variant_sequences)
    logger.info(
        "Kept %d/%d variant sequences after read coverage trimming to >=%dx",
        n_after_trimming,
        n_total,
        min_variant_sequence_coverage)
    return trimmed_variant_sequences


def filter_variant_sequences(
        variant_sequences,
        preferred_sequence_length,
        min_variant_sequence_coverage):
    """
    Trim variant sequences to the requested coverage and discard any which
    contain no usable sequence afterward.

    Candidate ranking is deliberately deferred until after reference matching
    and translation, when transcript mismatches and aggregate protein support
    are available.

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
