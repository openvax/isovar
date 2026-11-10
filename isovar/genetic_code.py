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
GeneticCode objects contain the rules for translating cDNA into a protein
sequence: the set of valid start and stop codons, as well as which
amino acid each DNA triplet is translated into.
"""

from functools import lru_cache


@lru_cache(maxsize=None)
def genetic_code_by_id(translation_table):
    """Use Varcode's offline NCBI definitions; unsupported/ambiguous codes raise."""
    from varcode.effects.codon_tables import codon_table_by_id

    table = codon_table_by_id(translation_table)
    return GeneticCode(table.name, table.start_codons, table.stop_codons,
                       table.forward_table)


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

        # Only complete codons in the supplied sequence can be decoded.
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

# Generic nuclear/fusion consumers use the same NCBI registry as explicit codes.
standard_genetic_code = genetic_code_by_id(1)


def translate_cdna(
        cdna_sequence,
        first_codon_is_start=False,
        *, translation_table=1):
    """
    Given a cDNA sequence which is aligned to a reading frame, returns
    the translated protein sequence and a boolean flag indicating whether
    the translated sequence ended on a stop codon (or just ran out of codons).

    Parameters
    ----------
    cdna_sequence : str
        cDNA sequence which is expected to start and end on complete codons.

    first_codon_is_start : bool
        Decode a recognized initiation codon as methionine at the first
        position only; internal codons retain their normal assignments.
    translation_table : int or None, optional
        NCBI table ID, defaulting to 1 (standard, including mammalian nuclear
        translation). Use 2 for vertebrate mitochondrial or 3 for yeast
        mitochondrial translation. None retains the existing
        table-1 decoding convention for undeclared annotated-ORF codes; it
        does not change their source metadata. Unsupported explicit codes
        and codes requiring context-dependent stop recoding fail explicitly.
    """
    genetic_code = genetic_code_by_id(1 if translation_table is None else translation_table)

    return genetic_code.translate(
        cdna_sequence=cdna_sequence,
        first_codon_is_start=first_codon_is_start)
