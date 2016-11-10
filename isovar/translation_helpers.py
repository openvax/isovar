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
This module combines variant cDNA sequences collected from a BAM file with
the reading frames of annotated reference transcripts to create candidate
translations.
"""

from __future__ import print_function, division, absolute_import

from .genetic_code import (
    standard_genetic_code_with_extra_start_codons,
    vertebrate_mitochondrial_genetic_code,
)

def compute_offset_to_first_complete_codon(
        offset_to_first_complete_reference_codon,
        n_trimmed_from_reference_sequence):
    """
    Once we've aligned the variant sequence to the ReferenceContext, we need
    to transfer reading frame from the reference transcripts to the variant
    sequences.

    Parameters
    ----------
    offset_to_first_complete_reference_codon : ReferenceContext

    n_trimmed_from_reference_sequence : int

    Returns an offset into the variant sequence that starts from a complete
    codon.
    """
    if n_trimmed_from_reference_sequence <= offset_to_first_complete_reference_codon:
        return (
            offset_to_first_complete_reference_codon -
            n_trimmed_from_reference_sequence)
    else:
        n_nucleotides_trimmed_after_first_codon = (
            n_trimmed_from_reference_sequence -
            offset_to_first_complete_reference_codon)
        frame = n_nucleotides_trimmed_after_first_codon % 3
        return (3 - frame) % 3

def translate_cdna(
        cdna_sequence,
        first_codon_is_start=False,
        mitochondrial=False):
    """
    Given a cDNA sequence which is aligned to a reading frame, returns
    the translated protein sequence and a boolean flag indicating whether
    the translated sequence ended on a stop codon (or just ran out of codons).

    Parameters
    ----------
    cdna_sequence : str
        cDNA sequence which is expected to start and end on complete codons.

    first_codon_is_start : bool

    mitochondrial : bool
        Use the mitochondrial codon table instead of standard
        codon to amino acid table.
    """
    # once we drop some of the prefix nucleotides, we should be in a reading frame
    # which allows us to translate this protein
    if mitochondrial:
        genetic_code = vertebrate_mitochondrial_genetic_code
    else:
        genetic_code = standard_genetic_code_with_extra_start_codons

    return genetic_code.translate(
        cdna_sequence=cdna_sequence,
        first_codon_is_start=first_codon_is_start)
