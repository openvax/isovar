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

from .common import groupby
from .logging import get_logger
from .protein_sequence import ProteinSequence
from .translation import  Translation

logger = get_logger(__name__)


def sort_protein_sequences(protein_sequences):
    """
    Sort protein sequences in decreasing order of priority
    """
    return list(
        sorted(
            protein_sequences,
            key=ProteinSequence.ascending_sort_key,
            reverse=True))


def summarize_translations(translations):
    """
    Summarize a collection of Translation objects into three values:
        1) List of unique reads supporting underlying variant sequences
        2) Set of unique transcript names for establishing reading frames of the
           translations.
        3) Set of unique gene names for all transcripts used by translations.
    """
    read_name_to_reads = {}
    gene_names = set([])
    transcript_ids = set([])
    for translation in translations:
        for read in translation.reads:
            read_name_to_reads[read.name] = read
        for transcript in translation.reference_context.transcripts:
            transcript_ids.add(transcript.id)
            gene_names.add(transcript.gene.name)
    unique_reads = list(read_name_to_reads.values())
    return unique_reads, transcript_ids, gene_names


def collapse_translations(translations):
    """
    Convert a list of Translation objects into a (potentially smaller) list
    of ProteinSequence objects by grouping the equivalent amino acid sequences.

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
            summarize_translations(equivalent_translations)

        logger.info(
            "%s: %s alt reads supporting protein sequence (gene names = %s)",
            key,
            len(alt_reads_supporting_protein_sequence),
            group_gene_names)

        protein_sequence = ProteinSequence.from_translation_key(
            translation_key=key,
            translations=equivalent_translations,
            transcript_ids_supporting_protein_sequence=group_transcript_ids,
            reads_supporting_protein_sequence=alt_reads_supporting_protein_sequence)
        logger.info("%s: protein sequence = %s" % (key, protein_sequence.amino_acids))
        protein_sequences.append(protein_sequence)
    return protein_sequences

