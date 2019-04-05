# Copyright (c) 2016-2018. Mount Sinai School of Medicine
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
IsovarResult is a collection of all information gathered about a variant
and any protein sequences which were successfully translated for it.
"""

from __future__ import print_function, division, absolute_import

from .value_object import  ValueObject

class VariantResult(ValueObject):
    """
    This object represents all information gathered about a variant,
    which includes the AlleleReads supporting any allele at this variant's
    locus and any protein sequences generated from an alt-allele cDNA
    assembly.
    """
    __slots__ = [
        "variant",
        "grouped_allele_reads",
        "sorted_protein_sequences",
    ]

    def __init__(self, variant, grouped_allele_reads, sorted_protein_sequences):
        self.variant = variant
        self.grouped_allele_reads = grouped_allele_reads
        self.sorted_protein_sequences = sorted_protein_sequences

    @property
    def ref_reads(self):
        return self.grouped_allele_reads.ref_reads

    @property
    def alt_reads(self):
        return self.grouped_allele_reads.alt_reads

    @property
    def other_reads(self):
        return self.grouped_allele_reads.other_reads

    @property
    def ref_read_names(self):
        return {r.name for r in self.ref_reads}

    @property
    def alt_read_names(self):
        return {r.name for r in self.alt_reads}

    @property
    def other_read_names(self):
        return {r.name for r in self.other_reads}

    @property
    def all_read_names(self):
        return self.ref_read_names.union(self.alt_read_names).union(self.other_read_names)

    @property
    def num_reads(self):
        return self.num_ref_reads + self.num_alt_reads + self.num_other_nonref_reads

    @property
    def num_fragments(self):
        return len(self.all_read_names)

    @property
    def ref_and_alt_read_names(self):
        return self.ref_read_names.union(self.alt_read_names)

    @property
    def num_ref_reads(self):
        return len(self.ref_reads)

    @property
    def num_ref_fragments(self):
        return len(self.ref_read_names)

    @property
    def num_alt_reads(self):
        return len(self.alt_reads)

    @property
    def num_alt_fragments(self):
        return len(self.alt_read_names)

    @property
    def num_other_reads(self):
        num_nonref = self.num_overlapping_reads - self.num_ref_reads
        return num_nonref - self.num_alt_reads

    @property
    def num_other_fragments(self):
        num_nonref = self.num_overlapping_fragments - self.num_ref_fragments
        return num_nonref - self.num_alt_fragments
