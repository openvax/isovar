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

from itertools import permutations

import pytest
from varcode import Variant

from isovar.allele_read import AlleleRead
from isovar.dna import reverse_complement_dna
from isovar.protein_sequence_creator import ProteinSequenceCreator
from isovar.reference_context import ReferenceContext
from isovar.variant_sequence_creator import VariantSequenceCreator


def _genomic_allele(cdna, strand):
    return cdna if strand == "+" else reverse_complement_dna(cdna)


def _read(prefix, alt, suffix, name, strand):
    if strand == "-":
        prefix, alt, suffix = (
            reverse_complement_dna(suffix),
            reverse_complement_dna(alt),
            reverse_complement_dna(prefix),
        )
    return AlleleRead(prefix=prefix, allele=alt, suffix=suffix, name=name)


def _context(variant, strand, prefix, suffix):
    return ReferenceContext(
        strand=strand,
        sequence_before_variant_locus=prefix,
        sequence_at_variant_locus=_genomic_allele(variant.ref, strand),
        sequence_after_variant_locus=suffix,
        offset_to_first_complete_codon=0,
        contains_start_codon=False,
        overlaps_start_codon=False,
        contains_five_prime_utr=False,
        amino_acids_before_variant="",
        variant=variant,
        transcripts=(),
    )


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("assembly", [False, True])
def test_multibase_substitution_marks_both_changed_amino_acids(strand, assembly):
    """#208: RNA-backed translation must expose both changes to vaxrank."""
    variant = Variant("1", 100, _genomic_allele("GAA", strand),
                      _genomic_allele("TCC", strand), "GRCh38")
    prefix, suffix = "AAA" * 4 + "AT", "AGGG"
    reads = [_read(prefix, "TCC", suffix, str(i), strand) for i in range(2)]
    sequences = VariantSequenceCreator(
        variant_sequence_assembly=assembly).reads_to_variant_sequences(variant, reads)
    context = _context(variant, strand, prefix, suffix)
    translation, = ProteinSequenceCreator().all_pairs_translations(sequences, [context])

    # NCBI standard code: ATG|AAA -> ATT|CCA, i.e. MK -> IP.
    # https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1
    assert translation.amino_acids == "KKKKIPG"
    assert (translation.mutation_start_idx, translation.mutation_end_idx) == (4, 6)
    assert translation.contains_mutation
    assert not translation.frameshift


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("assembly", [False, True])
def test_nested_single_read_candidates_keep_shared_support(strand, assembly):
    """#209: retained candidates need support, not just their original reads."""
    variant = Variant("1", 100, _genomic_allele("G", strand),
                      _genomic_allele("C", strand), "GRCh38")
    prefixes = [p + "A" * 12 for p in (
        "", "GGG", "CCCGGG", "TTTCCCGGG", "GGGTTTCCCGGG")]
    reads = [_read(prefix, "C", "A" * 30, str(i), strand)
             for i, prefix in enumerate(prefixes)]
    context = _context(variant, strand, "A" * 24, "A" * 30)
    creator = VariantSequenceCreator(variant_sequence_assembly=assembly)
    protein_creator = ProteinSequenceCreator()

    for ordered_reads in (reads, list(reversed(reads))):
        sequences = creator.reads_to_variant_sequences(variant, ordered_reads)
        core, = [sequence for sequence in sequences
                 if (sequence.prefix if strand == "+" else
                     reverse_complement_dna(sequence.suffix)) == "A" * 12]
        assert core.reads == frozenset(reads)
        assert core.min_coverage() == 5
        translations = protein_creator.all_pairs_translations(sequences, [context])
        assert translations
        assert {translation.amino_acids for translation in translations} == {"KKKKQKKKKKKKKK"}


@pytest.mark.parametrize("strand", ["+", "-"])
def test_competing_assembly_branches_translate_independently_of_read_order(strand):
    """#210: equal merge scores must not choose a peptide by BAM order."""
    variant = Variant("1", 100, _genomic_allele("G", strand),
                      _genomic_allele("C", strand), "GRCh38")
    groups = [
        [_read(prefix, "C", suffix, "%d_%d" % (i, j), strand) for j in range(2)]
        for i, (prefix, suffix) in enumerate([
            ("A" * 12, "G" * 21),
            ("A" * 9, "G" * 21 + "A" * 6),
            ("A" * 9, "G" * 21 + "T" * 6),
        ])
    ]
    context = _context(variant, strand, "A" * 12, "G" * 21 + "A" * 6)
    creator = VariantSequenceCreator(variant_sequence_assembly=True)
    protein_creator = ProteinSequenceCreator()
    expected = None
    for order in permutations(groups):
        reads = [read for group in order for read in group]
        sequences = creator.reads_to_variant_sequences(variant, reads)
        translations = protein_creator.all_pairs_translations(sequences, [context])
        observed = sorted(
            (t.amino_acids, tuple(sorted(r.name for r in t.reads))) for t in translations)
        assert any(len(amino_acids) == 13 for amino_acids, _ in observed)
        if expected is None:
            expected = observed
        assert observed == expected
