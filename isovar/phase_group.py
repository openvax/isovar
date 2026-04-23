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
Explicit representation of a phased set of variants observed together in RNA.
"""

from .value_object import ValueObject


class PhaseGroup(ValueObject):
    """
    A connected set of somatic variants that co-occur in RNA evidence.

    This is a group-level object rather than a focal-variant annotation. A
    variant can have pairwise phasing support with only a subset of the other
    variants in the group, while still belonging to the same connected
    component of the phasing graph.

    Germline variants are not populated yet, but the field is included so the
    public model can grow into that use case without changing shape again.
    """

    __slots__ = [
        "somatic_variants",
        "germline_variants",
        "supporting_read_names",
    ]

    def __init__(self, somatic_variants, germline_variants=(), supporting_read_names=()):
        self.somatic_variants = tuple(somatic_variants)
        self.germline_variants = tuple(germline_variants)
        self.supporting_read_names = frozenset(supporting_read_names)
