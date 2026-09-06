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

import pytest

from isovar.translation_helpers import find_mutant_amino_acid_interval


@pytest.mark.parametrize("codon_offset", [0, 1, 7])
@pytest.mark.parametrize("frame", range(3))
@pytest.mark.parametrize("n_ref,n_alt", [
    (n_ref, n_alt)
    for n_ref in range(10)
    for n_alt in range(10)
    if n_ref or n_alt
])
def test_mutation_interval_covers_affected_codons(codon_offset, frame, n_ref, n_alt):
    """#208: enumerate codons, independently of the interval arithmetic.

    Cover substitutions, insertions, deletions and delins at every codon
    position, including incomplete leading codons and untranslated prefixes.
    """
    coding_prefix_length = 12 + frame
    variant_start = codon_offset + coding_prefix_length
    cdna = "T" * codon_offset + "A" * coding_prefix_length + "C" * n_alt + "G" * 36
    n_amino_acids = (len(cdna) - codon_offset) // 3
    frameshift = (n_alt - n_ref) % 3 != 0

    if frameshift:
        expected_end = n_amino_acids
    else:
        affected_codons = {
            (coding_prefix_length + i) // 3 for i in range(n_alt)
        }
        # An in-frame deletion either joins two partial codons into one
        # affected codon, or leaves a zero-width junction between codons.
        if not n_alt and frame:
            affected_codons.add(coding_prefix_length // 3)
        expected_end = max(affected_codons) + 1 if affected_codons else 4

    assert find_mutant_amino_acid_interval(
        cdna_sequence=cdna,
        cdna_first_codon_offset=codon_offset,
        cdna_variant_start_offset=variant_start,
        cdna_variant_end_offset=variant_start + n_alt,
        n_ref=n_ref,
        n_amino_acids=n_amino_acids,
    ) == (4, expected_end, frameshift)
