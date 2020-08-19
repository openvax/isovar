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
Since multiple variant sequences can translate to the same amino acid sequence,
this module aggregates equivalent Translation objects into a single
ProteinSequence.
"""

from __future__ import print_function, division, absolute_import

from .common import groupby
from .logging import get_logger
from .protein_sequence import ProteinSequence
from .translation import Translation

logger = get_logger(__name__)


def sort_protein_sequences(protein_sequences):
    """
    Sort protein sequences in decreasing order of priority
    """
    return list(
        sorted(
            protein_sequences,
            key=ProteinSequence.ascending_sort_key,
            reverse=True))


def group_equivalent_translations(translations):
    """
    Convert a list of Translation objects into a (potentially smaller) list
    of ProteinSequence objects by grouping the equivalent amino acid sequences.

    Parameters
    ----------
    translations : list of Translation objects

    Returns list of ProteinSequence objects
    """
    protein_sequences = []
    translation_groups = groupby(
        translations,
        key_fn=Translation.as_translation_key)
    for equivalent_translations in translation_groups.values():
        protein_sequences.append(ProteinSequence.from_translations(equivalent_translations))
    return protein_sequences
