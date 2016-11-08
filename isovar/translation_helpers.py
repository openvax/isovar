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

human_standard_codon_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

# Non-canonical start sites based on figure 2 of
#   "Global mapping of translation initiation sites in mammalian
#   cells at single-nucleotide resolution"

human_standard_start_codons = {
    'ATG',
    'CTG',  # how often are CTG translated as leucine? skimming of papers
            # makes it seem rare and dependent on upstream sequence
    'TTG',
    'GTG',
    'AGG',
    'ACG',
    'AAG',
    'ATC',
    'ATA',
    'ATT',
}

def translate_cdna(
        cdna_sequence_from_first_codon,
        first_codon_is_start=False,
        mitochondrial=False):
    """
    Given a cDNA sequence which is aligned to a reading frame, returns
    the translated protein sequence and a boolean flag indicating whether
    the translated sequence ended on a stop codon (or just ran out of codons).

    Parameters
    ----------
    cdna_sequence_from_first_codon : skbio.DNA
        cDNA sequence which is expected to start and end on complete codons.
    """
    # once we drop some of the prefix nucleotides, we should be in a reading frame
    # which allows us to translate this protein
    if mitochondrial:
        raise ValueError("Translation of mitochondrial cDNA not yet supported")
    if not isinstance(cdna_sequence_from_first_codon, str):
        cdna_sequence_from_first_codon = str(cdna_sequence_from_first_codon)
    n = len(cdna_sequence_from_first_codon)
    if len(cdna_sequence_from_first_codon) < 3:
        return ''

    if first_codon_is_start and (
            cdna_sequence_from_first_codon[:3] in human_standard_start_codons):
        amino_acid_list = ['M']
        start_index = 3
    else:
        start_index = 0
        amino_acid_list = []

    # trim to multiple of 3 length
    end_idx = 3 * (n // 3)
    ends_with_stop_codon = False
    for i in range(start_index, end_idx, 3):
        codon = cdna_sequence_from_first_codon[i:i + 3]
        assert codon, (i, start_index, end_idx, cdna_sequence_from_first_codon)
        aa = human_standard_codon_table[codon]

        if aa == "*":
            ends_with_stop_codon = True
            break
        amino_acid_list.append(aa)

    amino_acids = "".join(amino_acid_list)

    return amino_acids, ends_with_stop_codon
