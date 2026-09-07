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

from numbers import Integral, Real

from .common import groupby
from .logging import get_logger
from .protein_sequence import ProteinSequence
from .translation import Translation

logger = get_logger(__name__)


def validate_protein_sequence_preference(preference, peptide_length, min_support_fraction):
    if preference not in ("balanced", "support", "context"):
        raise ValueError("protein_sequence_preference must be balanced, support or context")
    if isinstance(peptide_length, bool) or not isinstance(peptide_length, Integral) or peptide_length < 1:
        raise ValueError("protein_context_peptide_length must be a positive integer")
    if (isinstance(min_support_fraction, bool) or not isinstance(min_support_fraction, Real)
            or not 0 <= min_support_fraction <= 1):
        raise ValueError("min_protein_sequence_support_fraction must be finite and between 0 and 1")


def mutant_peptide_window_count(protein_sequence, peptide_length):
    """Number of full windows overlapping the edit, not the number of epitopes.

    A zero-width deletion junction must be strictly inside a window. Windows
    touching only its boundary contain no novel adjacency. Distinct start
    positions are counted even if repeated sequence makes their strings equal.
    """
    if isinstance(peptide_length, bool) or not isinstance(peptide_length, Integral) or peptide_length < 1:
        raise ValueError("peptide_length must be a positive integer")
    if not protein_sequence.contains_mutation:
        return 0
    start, end = protein_sequence.mutation_start_idx, protein_sequence.mutation_end_idx
    first_start = max(0, start - peptide_length + 1)
    last_start = min(len(protein_sequence.amino_acids) - peptide_length,
                     (end if start < end else start) - 1)
    return max(0, last_start - first_start + 1)


def sort_protein_sequences(protein_sequences, preference="support", peptide_length=25,
                          min_support_fraction=0.9):
    """
    Sort candidates without dropping reads or changing their support counts.

    The historical helper default is support-first; ProteinSequenceCreator
    explicitly requests its configured preference (balanced by default).
    Balanced maximizes useful peptide windows within a relative compatible
    read-name support budget. Context ignores that relative budget. Neither
    mode treats the budget as a calibrated confidence or expression estimate.
    """
    validate_protein_sequence_preference(preference, peptide_length, min_support_fraction)
    protein_sequences = list(protein_sequences)
    if preference == "support":
        return sorted(protein_sequences, key=ProteinSequence.ascending_sort_key, reverse=True)
    best_support = max((p.num_supporting_fragments for p in protein_sequences if p.contains_mutation), default=0)

    def key(p):
        # Compare the reported ratio directly: multiplying e.g. .14 * 50
        # rounds above 7 and would incorrectly reject exactly 7/50 support.
        within_budget = (best_support == 0
                         or p.num_supporting_fragments / best_support >= min_support_fraction)
        return (
            p.contains_mutation,
            within_budget if preference == "balanced" else True,
            mutant_peptide_window_count(p, peptide_length),
            # If no full window is available, retain useful partial context
            # within the same budget; never pretend it is a complete peptide.
            min(len(p.amino_acids), peptide_length),
            p.ascending_sort_key(),
            # Deterministic content ties, independent of input/hash order.
            p.amino_acids, p.mutation_start_idx, p.mutation_end_idx,
            p.ends_with_stop_codon, p.frameshift,
        )

    return sorted(protein_sequences, key=key, reverse=True)


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
