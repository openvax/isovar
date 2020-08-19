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

from .allele_read_helpers import get_single_allele_from_reads
from .assembly import iterative_overlap_assembly
from .default_parameters import (
    MIN_VARIANT_SEQUENCE_COVERAGE,
    VARIANT_SEQUENCE_LENGTH,
    VARIANT_SEQUENCE_ASSEMBLY,
    MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE
)
from .logging import get_logger
from .variant_sequence_helpers import (
    filter_variant_sequences,
    initial_variant_sequences_from_reads
)

logger = get_logger(__name__)


class VariantSequenceCreator(object):
    """
    Assembler is used to assemble a set of AlleleReads into a smaller set of
    VariantSequence objects based on overlap of read sequences.
    """
    def __init__(
            self,
            min_variant_sequence_coverage=MIN_VARIANT_SEQUENCE_COVERAGE,
            preferred_sequence_length=VARIANT_SEQUENCE_LENGTH,
            variant_sequence_assembly=VARIANT_SEQUENCE_ASSEMBLY,
            min_assembly_overlap_size=MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE):
        """
        Parameters
        ----------
        min_variant_sequence_coverage : int
            Minimum number of RNA reads supporting each nucleotide of the
            variant cDNA sequence

        preferred_sequence_length : int
            Total number of nucleotides in the assembled sequences, including
            variant nucleotides.

        variant_sequence_assembly : bool
            Construct variant sequences by merging overlapping reads. If False
            then variant sequences must be fully spanned by cDNA reads.

        min_assembly_overlap_size : int
            Minimum number of nucleotides shared by two sequences before they
            can be merged into a single VariantSequence object.

        """
        self.min_variant_sequence_coverage = min_variant_sequence_coverage
        self.preferred_sequence_length = preferred_sequence_length
        self.variant_sequence_assembly = variant_sequence_assembly
        self.min_assembly_overlap_size = min_assembly_overlap_size

    def reads_to_variant_sequences(
            self,
            variant,
            reads):
        """
        Collapse variant-supporting RNA reads into consensus sequences of
        approximately the preferred length (may differ at the ends of transcripts),
        filter consensus sequences by length and number of supporting RNA reads.

        Parameters
        ----------
        variant : varcode.Variant

        reads : list of AlleleRead objects
            Reads which support the variant allele

        Returns
        -------
        list of VariantSequence
        """
        # convert to list in case it's a generator
        variant_reads = list(reads)

        if len(variant_reads) == 0:
            return []

        alt_seq = get_single_allele_from_reads(variant_reads)

        # the number of context nucleotides on either side of the variant
        # is half the desired length (minus the number of variant nucleotides)
        n_alt_nucleotides = len(alt_seq)

        n_surrounding_nucleotides = self.preferred_sequence_length - n_alt_nucleotides
        max_nucleotides_after_variant = n_surrounding_nucleotides // 2

        # if the number of nucleotides we need isn't divisible by 2 then
        # prefer to have one more *before* the variant since we need the
        # prefix sequence to match against reference transcripts
        max_nucleotides_before_variant = (
                n_surrounding_nucleotides - max_nucleotides_after_variant)

        variant_sequences = initial_variant_sequences_from_reads(
            variant_reads=variant_reads,
            max_nucleotides_before_variant=max_nucleotides_before_variant,
            max_nucleotides_after_variant=max_nucleotides_after_variant)

        logger.info(
            "Initial pool of %d variant sequences (min length=%d, max length=%d)",
            len(variant_sequences),
            min(len(s) for s in variant_sequences),
            max(len(s) for s in variant_sequences))

        if self.variant_sequence_assembly:
            # this is a tricky parameter to set correctly:
            # by how many bases should two sequences overlap before
            # we merge, currently defaulting to either half the non-variant
            # nucleotides or the specified min_assembly_overlap_size
            # (whichever is smaller)
            min_overlap_size = min(
                self.min_assembly_overlap_size,
                n_surrounding_nucleotides // 2)
            variant_sequences = iterative_overlap_assembly(
                variant_sequences,
                min_overlap_size=min_overlap_size)

        if variant_sequences:
            logger.info(
                "After overlap assembly: %d variant sequences (min length=%d, max length=%d)",
                len(variant_sequences),
                min(len(s) for s in variant_sequences),
                max(len(s) for s in variant_sequences))
        else:
            logger.info("After overlap assembly: 0 variant sequences")
            return []

        variant_sequences = filter_variant_sequences(
            variant_sequences=variant_sequences,
            preferred_sequence_length=self.preferred_sequence_length,
            min_variant_sequence_coverage=self.min_variant_sequence_coverage)

        if variant_sequences:
            logger.info(
                ("After coverage & length filtering: %d variant sequences "
                 "(min length=%d, max length=%d)"),
                len(variant_sequences),
                min(len(s) for s in variant_sequences),
                max(len(s) for s in variant_sequences))
        else:
            logger.info("After coverage & length filtering: 0 variant sequences")
            return []

        # sort VariantSequence objects by decreasing order of supporting read
        # counts
        variant_sequences.sort(key=lambda vs: -len(vs.reads))
        return variant_sequences

    def sequences_from_alt_reads_generator(self, variant_and_reads_generator):
        """
        For each (variant, [AlleleRead]) pair in the input generator,
        collapse the reads into a list of VariantSequence objects.

        Parameters
        ----------
        variant_and_reads_generator : generator
            Sequence of Variant objects paired with a list of reads which
            support that variant.

        Yields pairs with the following fields:
            - Variant
            - list of VariantSequence objects
        """
        for variant, reads in variant_and_reads_generator:
            variant_sequences = self.reads_to_variant_sequences(
                variant=variant,
                reads=reads)
            yield variant, variant_sequences

    def sequences_from_read_evidence_generator(
            self, variant_and_read_evidence_generator):
        """
        Given a generator of (Variant, ReadEvidence) pairs, generate a
        sequence of (Variant, [VariantSequence]) pairs.
        """
        reads_gen = (
            (variant, read_evidence.alt_reads)
            for (variant, read_evidence) in
            variant_and_read_evidence_generator
        )
        return self.sequences_from_alt_reads_generator(reads_gen)
