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

"""
This module implements basic DNA functionality in Python strings to
to avoid having to depend on a bigger library such as BioPython.
"""

dna_complement_dictionary = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
}

dna_nucleotides = list(sorted(dna_complement_dictionary.keys()))

dna_nucleotide_to_index = {c: i for (i, c) in enumerate(dna_nucleotides)}

index_to_dna_nucleotide = {i: c for (i, c) in enumerate(dna_nucleotides)}


def complement_dna(seq):
    """
    Convert every A->T, T->A, C->G, G->C in a DNA sequence

    Parameters
    ----------
    seq : str

    Returns str
    """
    return "".join(dna_complement_dictionary[nt] for nt in seq)


def reverse_complement_dna(seq):
    """
    Reverse complement of a DNA sequence

    Parameters
    ----------
    seq : str

    Returns str
    """
    return complement_dna(seq)[::-1]
