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
Helper functions for working with RNA reads
"""

from __future__ import print_function, division, absolute_import

def get_single_allele_from_reads(allele_reads):
    """
    Given a sequence of AlleleRead objects, which are expected to all have
    the same allele, return that allele.
    """
    allele_reads = list(allele_reads)

    if len(allele_reads) == 0:
        raise ValueError("Expected non-empty list of AlleleRead objects")

    seq = allele_reads[0].allele
    if any(read.allele != seq for read in allele_reads):
        raise ValueError("Expected all AlleleRead objects to have same allele '%s', got %s" % (
            seq, allele_reads))
    return seq


def make_prefix_suffix_pairs(allele_reads):
    return [(r.prefix, r.suffix) for r in allele_reads]
