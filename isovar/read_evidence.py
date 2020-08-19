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
"""

from __future__ import print_function, division, absolute_import

from .allele_read_helpers import split_reads_into_ref_alt_other
from .variant_helpers import trim_variant
from .value_object import ValueObject


class ReadEvidence(ValueObject):
    """
    This class represents the reads at a variant locus partitioned
    by allele (ref/alt/other) relative to a variant.
    """

    __slots__ = [
        "trimmed_base1_start",
        "trimmed_ref",
        "trimmed_alt",
        "ref_reads",
        "alt_reads",
        "other_reads"
    ]

    @classmethod
    def from_variant_and_allele_reads(
            cls,
            variant,
            allele_reads):
        """
        Create a GroupedAlleleReads object from a variant and the set of reads overlapping
        the location of that variant.

        Parameters
        ----------
        variant : varcode.Variant

        allele_reads : list of AlleleRead

        Returns GroupedAlleleReads

        """
        trimmed_base1_start, trimmed_ref, trimmed_alt = \
            trim_variant(variant)
        ref_reads, alt_reads, other_reads = split_reads_into_ref_alt_other(
            ref=trimmed_ref,
            alt=trimmed_alt,
            overlapping_reads=allele_reads)
        return cls(
            trimmed_base1_start=trimmed_base1_start,
            trimmed_ref=trimmed_ref,
            trimmed_alt=trimmed_alt,
            ref_reads=ref_reads,
            alt_reads=alt_reads,
            other_reads=other_reads)

    def __init__(
            self,
            trimmed_base1_start,
            trimmed_ref,
            trimmed_alt,
            ref_reads,
            alt_reads,
            other_reads):
        """
        Parameters
        ----------
        trimmed_base1_start : int
            Base-1 inclusive genomic position indicating start of variant,
            trimmed to avoid any shared prefix between the ref and alt sequences

        trimmed_ref : str
            Reference allele, trimmed to avoid any shared prefix with the alt
            sequence.

        trimmed_alt : str
            Alternate allele, trimmed to avoid any shared prefix with the ref
            sequence.

        ref_reads : list of AlleleRead
            Reads supporting the reference allele

        alt_reads : list of AlleleRead
            Reads supporting the alt allele

        other_reads : list of AlleleRead
            Reads supporting some allele other than ref or alt.
        """
        self.trimmed_base1_start = trimmed_base1_start
        self.trimmed_ref = trimmed_ref
        self.trimmed_alt = trimmed_alt
        self.ref_reads = ref_reads
        self.alt_reads = alt_reads
        self.other_reads = other_reads

    @property
    def ref_read_names(self):
        """
        Names of reads which match the ref allele

        Returns set of str
        """
        return {r.name for r in self.ref_reads}

    @property
    def alt_read_names(self):
        """
        Names of reads which match the alt allele

        Returns set of str
        """
        return {r.name for r in self.alt_reads}

    @property
    def other_read_names(self):
        """
        Names of reads which match non-ref/non-alt alleles

        Returns set of str
        """
        return {r.name for r in self.other_reads}
