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
GeneticCode objects contain the rules for translating cDNA into a protein
sequence: the set of valid start and stop codons, as well as which
amino acid each DNA triplet is translated into.
"""


class GeneticCode(object):
    """
    Represents distinct translation tables to go from cDNA triplets to amino
    acids.
    """
    def __init__(self, name, start_codons, stop_codons, codon_table):
        self.name = name
        self.start_codons = set(start_codons)
        self.stop_codons = set(stop_codons)
        self.codon_table = dict(codon_table)
        self._check_codons()

    def _check_codons(self):
        """
        If codon table is missing stop codons, then add them.
        """
        for stop_codon in self.stop_codons:
            if stop_codon in self.codon_table:
                if self.codon_table[stop_codon] != "*":
                    raise ValueError(
                        ("Codon '%s' not found in stop_codons, but codon table "
                         "indicates that it should be") % (stop_codon,))
            else:
                self.codon_table[stop_codon] = "*"

        for start_codon in self.start_codons:
            if start_codon not in self.codon_table:
                raise ValueError(
                    "Start codon '%s' missing from codon table" % (
                        start_codon,))

        for codon, amino_acid in self.codon_table.items():
            if amino_acid == "*" and codon not in self.stop_codons:
                raise ValueError(
                    "Non-stop codon '%s' can't translate to '*'" % (
                        codon,))

        if len(self.codon_table) != 64:
            raise ValueError(
                "Expected 64 codons but found %d in codon table" % (
                    len(self.codon_table,)))

    def translate(self, cdna_sequence, first_codon_is_start=False):
        """
        Given a cDNA sequence which is aligned to a reading frame, returns
        the translated protein sequence and a boolean flag indicating whether
        the translated sequence ended on a stop codon (or just ran out of codons).

        Parameters
        ----------
        cdna_sequence : str
            cDNA sequence which is expected to start and end on complete codons.

        first_codon_is_start : bool
            Is the first codon of the sequence a start codon?
        """
        if not isinstance(cdna_sequence, str):
            cdna_sequence = str(cdna_sequence)
        n = len(cdna_sequence)

        # trim to multiple of 3 length, if there are 1 or 2 nucleotides
        # dangling at the end of an mRNA they will not affect translation
        # since ribosome will fall off at that point
        end_idx = 3 * (n // 3)

        codon_table = self.codon_table
        if first_codon_is_start and cdna_sequence[:3] in self.start_codons:
            amino_acid_list = ['M']
            start_index = 3
        else:
            start_index = 0
            amino_acid_list = []

        ends_with_stop_codon = False
        for i in range(start_index, end_idx, 3):
            codon = cdna_sequence[i:i + 3]
            aa = codon_table[codon]

            if aa == "*":
                ends_with_stop_codon = True
                break
            amino_acid_list.append(aa)

        amino_acids = "".join(amino_acid_list)
        return amino_acids, ends_with_stop_codon

    def copy(
            self,
            name,
            start_codons=None,
            stop_codons=None,
            codon_table=None,
            codon_table_changes=None):
        """
        Make copy of this GeneticCode object with optional replacement
        values for all fields.
        """
        new_start_codons = (
            self.start_codons.copy()
            if start_codons is None
            else start_codons)

        new_stop_codons = (
            self.stop_codons.copy()
            if stop_codons is None
            else stop_codons)

        new_codon_table = (
            self.codon_table.copy()
            if codon_table is None
            else codon_table)

        if codon_table_changes is not None:
            new_codon_table.update(codon_table_changes)

        return GeneticCode(
            name=name,
            start_codons=new_start_codons,
            stop_codons=new_stop_codons,
            codon_table=new_codon_table)

standard_genetic_code = GeneticCode(
    name="standard",
    start_codons={'ATG', 'CTG', 'TTG'},
    stop_codons={'TAA', 'TAG', 'TGA'},
    codon_table={
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
)

# Non-canonical start sites based on figure 2 of
#   "Global mapping of translation initiation sites in mammalian
#   cells at single-nucleotide resolution"
standard_genetic_code_with_extra_start_codons = standard_genetic_code.copy(
    name="standard-with-extra-start-codons",
    start_codons=standard_genetic_code.start_codons.union({
        'GTG',
        'AGG',
        'ACG',
        'AAG',
        'ATC',
        'ATA',
        'ATT'}))

vertebrate_mitochondrial_genetic_code = standard_genetic_code.copy(
    name="verterbrate-mitochondrial",
    # "For thirty years AGA and AGG were considered terminators instead
    #  of coding for arginine. However, Temperley (2010) has recently shown
    #  that human mitochondria use only UAA and UAG stop codons."
    # (http://mitomap.org/bin/view.pl/MITOMAP/HumanMitoCode)
    stop_codons={'TAA', 'TAG'},
    # "AUU codes for isoleucine during elongation but can code for
    #  methionine for initiation (ND2) See Fearnley & Walker (1987) and
    #  Peabody (1989)."
    # (http://mitomap.org/bin/view.pl/MITOMAP/HumanMitoCode)
    start_codons=['ATT', 'ATC', 'ATA', 'ATG', 'GTG'],
    # "UGA codes for tryptophan instead of termination and AUA codes for
    #  methionine instead of isoleucine."
    # (http://mitomap.org/bin/view.pl/MITOMAP/HumanMitoCode)
    codon_table_changes={'TGA': 'W', 'ATA': 'M'},
)


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
