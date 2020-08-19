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
    """
    Extremely simplified representation of a read at a locus: just the allele
    at the locus and sequence before/after. We're ignoring the base qualities
    and any additional information about splicing, clipping or alignment.
    """
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
    def from_locus_read(cls, locus_read):
        """
        Given a single LocusRead object, return either an AlleleRead or None

        Parameters
        ----------
        locus_read : LocusRead
            Read which overlaps a variant locus but doesn't necessarily contain the
            alternate nucleotides
        """
        sequence = locus_read.sequence
        read_name = locus_read.name

        reference_base0_start_inclusive = locus_read.reference_base0_start_inclusive
        reference_base0_end_exclusive = locus_read.reference_base0_end_exclusive

        read_base0_start_inclusive = locus_read.read_base0_start_inclusive
        read_base0_end_exclusive = locus_read.read_base0_end_exclusive

        if read_base0_start_inclusive is None or read_base0_end_exclusive is None:
            logger.debug(
                "Skipping read '%s' because required bases in reference interval %s:%s aren't mapped",
                read_name,
                reference_base0_start_inclusive,
                reference_base0_end_exclusive)
            return None

        reference_positions = locus_read.reference_positions

        n_ref_bases = reference_base0_end_exclusive - reference_base0_start_inclusive

        insertion = (n_ref_bases == 0)

        if insertion:
            # insertions require a sequence of non-aligned bases
            # followed by the subsequence reference position
            for read_index in range(read_base0_start_inclusive, read_base0_end_exclusive):
                # all the inserted nucleotides should *not* align to the reference
                if reference_positions[read_index] is not None:
                    logger.debug(
                        "Skipping read '%s', inserted nucleotides shouldn't map to reference",
                        read_name)
                    return None

        nucleotides_at_variant_locus = convert_from_bytes_if_necessary(
            sequence[read_base0_start_inclusive:read_base0_end_exclusive])

        if "N" in nucleotides_at_variant_locus:
            logger.debug(
                "Skipping read '%s', found N nucleotides at variant locus",
                read_name)
        prefix = convert_from_bytes_if_necessary(sequence[:read_base0_start_inclusive])
        suffix = convert_from_bytes_if_necessary(sequence[read_base0_end_exclusive:])

        prefix, suffix = trim_N_nucleotides(prefix, suffix)

        return AlleleRead(
            prefix,
            nucleotides_at_variant_locus,
            suffix,
            name=read_name)
