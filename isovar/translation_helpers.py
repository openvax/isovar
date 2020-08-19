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
Helper functions used for creating translating a variant's cDNA sequence
into a particular reading frame.
"""

from __future__ import print_function, division, absolute_import

import math

def find_mutant_amino_acid_interval(
        cdna_sequence,
        cdna_first_codon_offset,
        cdna_variant_start_offset,
        cdna_variant_end_offset,
        n_ref,
        n_amino_acids):
    """
    Parameters
    ----------
    cdna_sequence : skbio.DNA or str
        cDNA sequence found in RNAseq data

    cdna_first_codon_offset : int
        Offset into cDNA sequence to first complete codon, lets us skip
        past UTR region and incomplete codons.

    cdna_variant_start_offset : int
        Interbase start offset into cDNA sequence for selecting mutant
        nucleotides.

    cdna_variant_end_offset : int
        Interbase end offset into cDNA sequence for selecting mutant
        nucleotides.

    n_ref : int
        Number of reference nucleotides

    n_amino_acids : int
        Number of translated amino acids

    Returns
    -------
    tuple with three fields:
        1) Start offset for interval of mutant amino acids in translated sequence
        2) End offset for interval of mutant amino acids in translated sequence
        3) Boolean flag indicating whether the variant was a frameshift.
    """
    cdna_alt_nucleotides = cdna_sequence[
        cdna_variant_start_offset:cdna_variant_end_offset]

    n_alt = len(cdna_alt_nucleotides)

    # sequence of nucleotides before the variant starting from the first codon
    cdna_coding_prefix = cdna_sequence[cdna_first_codon_offset:cdna_variant_start_offset]

    # rounding down since a change in the middle of a codon should count
    # toward the variant codons
    n_coding_nucleotides_before_variant = len(cdna_coding_prefix)

    n_complete_prefix_codons = n_coding_nucleotides_before_variant // 3

    frame_of_variant_nucleotides = n_coding_nucleotides_before_variant % 3
    frameshift = abs(n_ref - n_alt) % 3 != 0
    indel = n_ref != n_alt

    variant_aa_interval_start = n_complete_prefix_codons

    if frameshift:
        # if mutation is a frame shift then every amino acid from the
        # first affected codon to the stop is considered mutant
        #
        # TODO: what if the first k amino acids are synonymous with the reference sequence?
        variant_aa_interval_end = n_amino_acids
    else:
        n_alt_codons = int(math.ceil(n_alt / 3.0))
        if indel:
            # We need to adjust the number of affected codons by whether the
            # variant is aligned with codon boundaries, since in-frame indels
            # may still be split across multiple codons.
            #
            # Example of in-frame deletion of 3 nucleotides which leaves
            # 0 variant codons in the sequence (interval = 1:1)
            #   ref = CCC|AAA|GGG|TTT
            #   alt = CCC|GGG|TTT
            #
            # Example of in-frame deletion of 3 nucleotides which leaves
            # 1 variant codon in the sequence (interval = 1:2)
            #   ref = CCC|AAA|GGG|TTT
            #   alt = CCC|AGG|TTT
            #
            # Example of in-frame insertion of 3 nucleotides which
            # yields two variant codons:
            #   ref = CCC|AAA|GGG|TTT
            #   alt = CTT|TCC|AAA|GGG|TTT
            extra_affected_codon = int(frame_of_variant_nucleotides != 0)
            variant_aa_interval_end = (
                variant_aa_interval_start + n_alt_codons + extra_affected_codon)
        else:
            # if the variant is a simple substitution then it only affects
            # as many codons as are in the alternate sequence
            variant_aa_interval_end = variant_aa_interval_start + n_alt_codons
    return variant_aa_interval_start, variant_aa_interval_end, frameshift
