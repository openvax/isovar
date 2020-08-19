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

from __future__ import print_function, division, absolute_import

from .logging import get_logger
from .variant_orf import VariantORF

logger = get_logger(__name__)


def match_variant_sequence_to_reference_context(
        variant_sequence,
        reference_context,
        min_transcript_prefix_length,
        max_transcript_mismatches,
        count_mismatches_after_variant=False,
        max_trimming_attempts=2):
    """
    Iteratively trim low-coverage subsequences of a variant sequence
    until it either matches the given reference context or there
    are too few nucleotides left in the variant sequence.

    Parameters
    ----------
    variant_sequence : VariantSequence
        Assembled sequence from RNA reads, will need to be to be reverse
        complemented if matching against a reference transcript on the
        negative strand.

    reference_context : ReferenceContext
        Sequence of reference transcript before the variant and associated
        metadata.

    min_transcript_prefix_length : int
        Minimum number of nucleotides we try to match against a reference
        transcript.

    max_transcript_mismatches : int
        Maximum number of nucleotide differences between reference transcript
        sequence and the variant sequence.

    count_mismatches_after_variant : bool
        Set to true if the number of mismatches after the variant locus should
        count toward the total max_transcript_mismatches, which by default
        only counts mismatches before the variant locus.

    max_trimming_attempts : int
        How many times do we try trimming the VariantSequence to higher
        levels of coverage before giving up?

    Returns VariantORF or None
    """
    # if we can't get the variant sequence to match this reference
    # context then keep trimming it by coverage until either
    for i in range(max_trimming_attempts + 1):
        # check the reverse-complemented prefix if the reference context is
        # on the negative strand since variant sequence is aligned to
        # genomic DNA (positive strand)
        variant_sequence_too_short = (
            (reference_context.strand == "+" and
                len(variant_sequence.prefix) < min_transcript_prefix_length) or
            (reference_context.strand == "-" and
                len(variant_sequence.suffix) < min_transcript_prefix_length)
        )
        if variant_sequence_too_short:
            logger.info(
                "Prefix of variant sequence %s shorter than min allowed %d (iter=%d)",
                variant_sequence,
                min_transcript_prefix_length,
                i + 1)
            return None

        variant_orf = \
            VariantORF.from_variant_sequence_and_reference_context(
                variant_sequence=variant_sequence,
                reference_context=reference_context)

        if variant_orf is None:
            return None

        n_mismatch_before_variant = (
            variant_orf.num_mismatches_before_variant)
        n_mismatch_after_variant = (
            variant_orf.num_mismatches_after_variant)

        logger.info("Iter #%d/%d: %s (len=%d)" % (
            i + 1,
            max_trimming_attempts + 1,
            variant_orf,
            len(variant_orf.cdna_sequence)))

        total_mismatches = n_mismatch_before_variant
        if count_mismatches_after_variant:
            total_mismatches += n_mismatch_after_variant
        if total_mismatches <= max_transcript_mismatches:
            # if we got a variant sequence + reading frame with sufficiently
            # few mismatches then call it a day
            return variant_orf

        logger.info(
            ("Too many mismatches (%d) between variant sequence %s and "
             "reference context %s (attempt=%d/%d)"),
            n_mismatch_before_variant,
            variant_sequence,
            reference_context,
            i + 1,
            max_trimming_attempts + 1)
        # if portions of the sequence are supported by only 1 read
        # then try trimming to 2 to see if the better supported
        # subsequence can be better matched against the reference
        current_min_coverage = variant_sequence.min_coverage()
        logger.info(
            "Trimming to subsequence covered by at least %d reads",
            current_min_coverage + 1)
        variant_sequence = variant_sequence.trim_by_coverage(
            current_min_coverage + 1)
    return None
