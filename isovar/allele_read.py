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
Reads overlapping a locus of interest split into prefix,
allele (ref, alt, or otherwise), and suffix portions
"""

from __future__ import print_function, division, absolute_import
import logging

from .string_helpers import convert_from_bytes_if_necessary, trim_N_nucleotides
from .value_object import ValueObject

logger = logging.getLogger(__name__)

class AlleleRead(ValueObject):
    __slots__ = ["prefix", "allele", "suffix", "name", "sequence"]

    def __init__(self, prefix, allele, suffix, name):
        self.prefix = prefix
        self.allele = allele
        self.suffix = suffix
        self.name = name
        self.sequence = prefix + allele + suffix

    def __len__(self):
        return len(self.sequence)

    @classmethod
    def from_locus_read(cls, locus_read, n_ref):
        """
        Given a single LocusRead object, return either an AlleleRead or None

        Parameters
        ----------
        locus_read : LocusRead
            Read which overlaps a variant locus but doesn't necessarily contain the
            alternate nucleotides

        n_ref : int
            Number of reference positions we are expecting to be modified or
            deleted (for insertions this should be 0)
        """
        sequence = locus_read.sequence
        reference_positions = locus_read.reference_positions

        # positions of the nucleotides before and after the variant within
        # the read sequence
        read_pos_before = locus_read.base0_read_position_before_variant
        read_pos_after = locus_read.base0_read_position_after_variant

        # positions of the nucleotides before and after the variant on the
        # reference genome
        ref_pos_before = reference_positions[read_pos_before]

        if ref_pos_before is None:
            logger.warn(
                "Missing reference pos for nucleotide before variant on read: %s",
                locus_read)
            return None

        ref_pos_after = reference_positions[read_pos_after]

        if ref_pos_after is None:
            logger.warn(
                "Missing reference pos for nucleotide after variant on read: %s",
                locus_read)
            return None

        if n_ref == 0:
            if ref_pos_after - ref_pos_before != 1:
                # if the number of nucleotides skipped isn't the same
                # as the number of reference nucleotides in the variant then
                # don't use this read
                logger.debug(
                    "Positions before (%d) and after (%d) variant should be adjacent on read %s",
                    ref_pos_before,
                    ref_pos_after,
                    locus_read)
                return None

            # insertions require a sequence of non-aligned bases
            # followed by the subsequence reference position
            ref_positions_for_inserted = reference_positions[
                read_pos_before + 1:read_pos_after]
            if any(insert_pos is not None for insert_pos in ref_positions_for_inserted):
                # all these inserted nucleotides should *not* align to the
                # reference
                logger.debug(
                    "Skipping read, inserted nucleotides shouldn't map to reference")
                return None
        else:
            # substitutions and deletions
            if ref_pos_after - ref_pos_before != n_ref + 1:
                # if the number of nucleotides skipped isn't the same
                # as the number of reference nucleotides in the variant then
                # don't use this read
                logger.debug(
                    ("Positions before (%d) and after (%d) variant should be "
                     "adjacent on read %s"),
                    ref_pos_before,
                    ref_pos_after,
                    locus_read)
                return None

        nucleotides_at_variant_locus = sequence[read_pos_before + 1:read_pos_after]

        prefix = sequence[:read_pos_before + 1]
        suffix = sequence[read_pos_after:]

        prefix, suffix = convert_from_bytes_if_necessary(prefix, suffix)
        prefix, suffix = trim_N_nucleotides(prefix, suffix)

        return cls(
            prefix,
            nucleotides_at_variant_locus,
            suffix,
            name=locus_read.name)

