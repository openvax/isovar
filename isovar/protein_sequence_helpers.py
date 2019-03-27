# Copyright (c) 2019. Mount Sinai School of Medicine
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
Since multiple variant sequences can translate to the same amino acid sequence,
this module aggregates equivalent Translation objects into a single
ProteinSequence.
"""

from __future__ import print_function, division, absolute_import


def protein_sequence_ascending_sort_key(protein_sequence):
    """
    Sort protein sequences lexicographically by three criteria:
        - number of unique supporting reads
        - minimum mismatch versus a supporting reference transcript before variant
        - minimum mismatch versus a supporting reference transcript after variant
        - number of supporting reference transcripts
    """
    return (
        len(protein_sequence.alt_reads_supporting_protein_sequence),
        min(t.number_mismatches_before_variant for t in protein_sequence.translations),
        min(t.number_mismatches_after_variant for t in protein_sequence.translations),
        len(protein_sequence.transcripts_supporting_protein_sequence)
    )

def collapse_translations(translations):
    """

    Parameters
    ----------
    translations : list of Translation objects

    Returns list of ProteinSequence objects
    """
    protein_sequences = []
    for (key, equivalent_translations) in groupby(
            translations, key_fn=Translation.as_translation_key).items():
        # get the variant read names, transcript IDs and gene names for
        # protein sequence we're about to construct
        alt_reads_supporting_protein_sequence, group_transcript_ids, group_gene_names = \
            ProteinSequence._summarize_translations(equivalent_translations)

        logger.info(
            "%s: %s alt reads supporting protein sequence (gene names = %s)",
            key,
            len(alt_reads_supporting_protein_sequence),
            group_gene_names)

        protein_sequence = ProteinSequence.from_translation_key(
            translation_key=key,
            translations=equivalent_translations,
            overlapping_reads=overlapping_reads,
            alt_reads=variant_support.alt_reads,
            ref_reads=variant_support.ref_reads,
            alt_reads_supporting_protein_sequence=alt_reads_supporting_protein_sequence,
            transcripts_supporting_protein_sequence=group_transcript_ids,
            transcripts_overlapping_variant=overlapping_transcript_ids,
            gene=list(group_gene_names))
        logger.info("%s: protein sequence = %s" % (key, protein_sequence.amino_acids))
        protein_sequences.append(protein_sequence)


def sort_protein_sequences(protein_sequences):
    """
    Sort protein sequences in decreasing order of priority
    """
    return list(
        sorted(
            protein_sequences,
            key=protein_sequence_ascending_sort_key,
            reverse=True))
