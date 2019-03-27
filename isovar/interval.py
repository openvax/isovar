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


from __future__ import print_function, division, absolute_import

from .value_object import ValueObject
from .variant_helpers import base0_interval_for_variant


class Interval(ValueObject):
    """
    Representation of a half-open interval on a chromosome with helpers
    for converting to other coordinate systems.
    """
    __slots__ = [
        "chromosome",
        "base0_start_inclusive",
        "base0_end_exclusive",
    ]

    def __len__(self):
        return self.base0_end_exclusive - self.base0_start_inclusive

    @property
    def base0_half_open_coordinates(self):
        return (self.base0_start_inclusive,  self.base0_end_exclusive)

    @property
    def base1_half_open_coordinates(self):
        return (self.base1_start_inclusive, self.base1_end_exclusive)

    @property
    def base0_open_coorindates(self):
        if self.base0_start_inclusive == self.base0_end_exclusive:
            return (self.base0_start_inclusive, self.base0_start_inclusive + 1)
        return (self.base0_start_exclusive, self.base0_end_exclusive)

    @property
    def base1_open_coorindates(self):
        return (self.base1_start_exclusive, self.base1_end_exclusive)

    @property
    def base0_closed_coorindates(self):
        return (self.base0_start_exclusive, self.base0_end_exclusive)

    @property
    def base1_closed_coorindates(self):
        return (self.base1_start_exclusive, self.base1_end_exclusive)


    @property
    def base0_start_exclusive(self):
        return self.base0_position_before

    @property
    def base0_position_before(self):
        return self.base0_start_inclusive - 1

    @property
    def base0_position_after(self):
        return self.base0_end_exclusive

    @property
    def base1_position_before(self):
        return self.base0_start_inclusive

    @property
    def base1_position_after(self):
        return self.base0_end_exclusive + 1

    @property
    def base1_start_inclusive(self):
        return self.base0_start_inclusive + 1

    @property
    def base1_end_exclusive(self):
        return self.base0_end_exclusive + 1

    @classmethod
    def from_variant(cls, variant):
        base0_start_inclusive, base0_end_exclusive = base0_interval_for_variant(variant)
        return cls(
            chromosome=variant.contig,
            base0_start_inclusive=base0_start_inclusive,
            base0_end_exclusive=base0_end_exclusive)
