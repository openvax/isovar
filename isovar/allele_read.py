# Copyright (c) 2016-2019. Mount Sinai School of Medicine
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
