# Copyright (c) 2018-2019. Mount Sinai School of Medicine
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
"""



from .allele_read_helpers import  split_reads_into_ref_alt_other
from .coverage_stats import CoverageStats
from .variant_helpers import trim_variant
from .value_object import ValueObject


class VariantSupport(ValueObject):
    """
    This class represents the reads at a variant locus partitioned
    by allele (ref/alt/other).
    """

    __slots__ = [
        "variant",
        "reads",
        "trimemd_start",
        "trimmed_ref",
        "trimmed_alt",
        "ref_reads",
        "alt_reads",
        "other_reads"
    ]

    def __init__(
            self,
            variant,
            reads):
        self.variant = variant
        self.reads = reads

        trimmed_start, trimmed_ref, trimmed_alt = \
            trim_variant(variant)
        ref_reads, alt_reads, other_reads = split_reads_into_ref_alt_other(
            ref=trimmed_ref,
            alt=trimmed_alt,
            overlapping_reads=reads)
        self.trimmed_start = trimmed_start
        self.trimmed_ref = trimmed_ref
        self.trimmed_alt = trimmed_alt
        self.ref_reads = ref_reads
        self.alt_reads = alt_reads
        self.other_reads = other_reads

    def coverage_stats(self):
        """
        Create a CoverageStats object with counts of
        the number of reads and distinct fragments covering
        the ref, alt, and other alleles. This object can then
        be used for further filtering.

        Returns CoverageStats
        """
        return CoverageStats.from_reads(
            ref_reads=self.ref_reads,
            alt_reads=self.alt_reads,
            other_reads=self.other_reads)

