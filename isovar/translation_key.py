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
This module combines variant cDNA sequences collected from a BAM file with
the reading frames of annotated reference transcripts to create candidate
translations.
"""


from __future__ import print_function, division, absolute_import

from .logging import get_logger
from .value_object import ValueObject


logger = get_logger(__name__)


class TranslationKey(ValueObject):
    """
    TranslationKey contains fields related to a translated protein sequence
    which should be used to combine multiple equivalent mutated amino acid
    sequences.
    """
    __slots__ = [
        # translated sequence of a variant sequence in the ORF established
        # by a reference context
        "amino_acids",
        # is there a mutation in the amino acid sequence?
        "contains_mutation",
        # half-open interval coordinates for variant amino acids
        # in the translated sequence
        "mutation_start_idx",
        "mutation_end_idx",
        # did the amino acid sequence end due to a stop codon or did we
        # just run out of sequence context around the variant?
        "ends_with_stop_codon",
        # was the variant a frameshift relative to the reference sequence?
        "frameshift"
    ]
