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
    Drop variant sequences shorter than the longest one we managed to
    construct (capped at preferred_sequence_length), except for shorter
    sequences with strictly better fragment support.

    Length on its own is a poor way to choose between candidates once the
    requested sequence is longer than a single read, since assembled length
    then reflects how far a chain of overlapping reads happened to extend
    rather than how well supported that chain is. A long sequence built from
    two reads can evict a shorter one built from ten.

    Isovar's actual preference order is spelled out by
    ProteinSequence.ascending_sort_key: supporting fragments, then supporting
    reads, then mismatches against the reference transcript, and only then
    length as a tie-break. We can't apply all of it here because mismatch
    counts don't exist until a VariantSequence has been matched against a
    ReferenceContext. So instead of committing to length at this stage, keep
    any better-supported shorter candidate and let the downstream ranking in
    protein_sequence_helpers.sort_protein_sequences decide.

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

    long_enough = [
        s for s in variant_sequences
        if len(s.sequence) >= min_required_sequence_length
    ]
    # number of distinct fragments (not reads) behind the best long sequence,
    # matching the unit used by ProteinSequence.num_supporting_fragments
    best_long_support = max(len(s.read_names) for s in long_enough)
    # a zero length sequence can't be rescued no matter how many reads it
    # claims, since trim_by_coverage hands back all of the original reads
    # along with an empty sequence
    better_supported = [
        s for s in variant_sequences
        if len(s.sequence) < min_required_sequence_length
        and len(s.sequence) > 0
        and len(s.read_names) > best_long_support
    ]
    variant_sequences = long_enough + better_supported

    n_dropped = n_total - len(variant_sequences)
    if n_dropped > 0:
        logger.info(
            "Dropped %d/%d variant sequences shorter than %d",
            n_dropped,
            n_total,
            min_required_sequence_length)
    for s in better_supported:
        logger.info(
            ("Keeping variant sequence shorter than %d (len=%d) since its %d "
             "supporting fragments beat the %d supporting the longest sequence"),
            min_required_sequence_length,
            len(s.sequence),
            len(s.read_names),
            best_long_support)
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

