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


from __future__ import print_function, division, absolute_import

from .protein_sequence_helpers import sort_protein_sequences, collapse_translations
from .protein_sequence_creator import ProteinSequenceCreator
from .read_collector import ReadCollector
from .logging import get_logger
from .isovar_result import IsovarResult
from .value_object import ValueObject

logger = get_logger(__name__)

class Isovar(ValueObject):
    """
    This is the main entrypoint into the Isovar library, which collects
    RNA reads supporting variants and translates their coding sequence
    into amino acid sequences.
    """
    def __init__(
            self,
            read_collector=None,
            protein_sequence_creator=None):

        """
        read_collector : ReadCollector or None
            Object used to collect ReadEvidence for each variant, created
            with default settings if not supplied.

        protein_sequence_creator : ProteinSequenceCreator or None
            Object used to turn (Variant, ReadEvidence) into one or more
            ProteinSequence objects. Created with default settings if not
            supplied.
        """
        if read_collector is None:
            read_collector = ReadCollector()

        self.read_collector = read_collector

        if protein_sequence_creator is None:
            self.protein_sequence_creator = ProteinSequenceCreator()

    def process_variants(
            self,
            variants,
            alignment_file,
            transcript_id_whitelist=None):
        """
        Main entry-point into Isovar, which attempts to translate each
        variant in a protein sequence using supporting reads from the
        alignment_file.

        Parameters
        ----------
        variants : varcode.VariantCollection
            Somatic variants

        alignment_file : pysam.AlignmentFile
            Aligned tumor RNA reads

        transcript_id_whitelist : set of str or None
            Which transcripts should be considered when predicting DNA-only
            coding effects of mutations and also when trying to establish a
            reading frame for identified cDNA sequences.

        Generator of IsovarResult objects, one for each variant.
        The `protein_sequences` field of the IsovarVar result will be empty
        if no sequences could be determined.
        """

        # create generator which returns (Variant, ReadEvidence) pairs
        read_evidence_gen = \
            self.read_collector.read_evidence_generator(
               variants=variants,
               alignment_file=alignment_file)

        for variant, read_evidence in read_evidence_gen:
            # generate protein sequences by assembling variant reads
            protein_sequences = \
                self.sorted_protein_sequences_for_variant(
                    variant=variant,
                    read_evidence=read_evidence,
                    transcript_id_whitelist=transcript_id_whitelist)

            # if we're only keeping the top-k, get rid of the rest
            if self.max_protein_sequences_per_variant:
                protein_sequences = protein_sequences[:self.max_protein_sequences_per_variant]

            yield IsovarResult(
                variant=variant,
                transcript_id_whitelist=transcript_id_whitelist,
                read_evidence=read_evidence,
                protein_sequences=protein_sequences)