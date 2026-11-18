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

from isovar.dna import reverse_complement_dna, complement_dna
from .common import eq_


def test_reverse_complement_dna():
    eq_("ATGCAATTGGCC", reverse_complement_dna("GGCCAATTGCAT"))


def test_complement_dna():
    eq_("ATGC", complement_dna("TACG"))


def test_reverse_complement_handles_n_case_and_ambiguity_codes():
    eq_("NnRyAC", reverse_complement_dna("GTrYnN"))


def test_shared_complement_covers_iupac_and_preserves_case():
    from isovar.sv_rna_path_support import reverse_complement
    # INSDC Feature Table Appendix IV 7.4.1; ambiguity complements follow
    # the complement of each represented base set (e.g. R=A/G -> Y=T/C).
    sequence = 'ACGTRYKMSWBDHVNacgtrykmswbdhvn'
    complement = 'TGCAYRMKSWVHDBNtgcayrmkswvhdbn'
    assert complement_dna(sequence) == complement
    assert reverse_complement_dna(sequence) == complement[::-1]
    assert reverse_complement(sequence) == complement[::-1]
    assert reverse_complement(reverse_complement(sequence)) == sequence
