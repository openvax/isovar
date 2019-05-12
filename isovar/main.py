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

from six import string_types
import pandas as pd
from varcode import load_vcf
from pysam import AlignmentFile

from .filtering import DEFAULT_FILTER_THRESHOLDS, apply_filters
from .protein_sequence_creator import ProteinSequenceCreator
from .read_collector import ReadCollector
from .logging import get_logger
from .isovar_result import IsovarResult
from .value_object import ValueObject

logger = get_logger(__name__)


def run_isovar(
        variants,
        alignment_file,
        transcript_id_whitelist=None,
        read_collector=None,
        protein_sequence_creator=None,
        filter_thresholds=DEFAULT_FILTER_THRESHOLDS):
    """
    This is the main entrypoint into the Isovar library, which collects
    RNA reads supporting variants and translates their coding sequence
    into amino acid sequences. Collects both the read evidence and
    protein sequences into IsovarResult objects. The values of any filters
    which are supplied in the filter_thresholds argument are attached to
    each IsovarResult's filter_values_dict field.

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

    read_collector : ReadCollector or None
        Object used to collect ReadEvidence for each variant, created
        with default settings if not supplied.

    protein_sequence_creator : ProteinSequenceCreator or None
        Object used to turn (Variant, ReadEvidence) into one or more
        ProteinSequence objects. Created with default settings if not
        supplied.

    filter_thresholds : dict or OrderedDict
        Dictionary whose entries have names like "min_num_alt_reads"
        mapping to a numerical threshold value. In general, the keys
        must start with either "min_" or "max_" followed by a property
        of the IsovarResult class.

    Generator of IsovarResult objects, one for each variant. The
    `protein_sequences` field of the IsovarVar result will be empty
    if no sequences could be determined.
    """
    if isinstance(variants, string_types):
        variants = load_vcf(variants)

    if isinstance(alignment_file, string_types):
        alignment_file = AlignmentFile(alignment_file)

    if read_collector is None:
        read_collector = ReadCollector()

    if protein_sequence_creator is None:
        protein_sequence_creator = ProteinSequenceCreator()

    # create generator which returns (Variant, ReadEvidence) pairs
    read_evidence_gen = \
        read_collector.read_evidence_generator(
           variants=variants,
           alignment_file=alignment_file)

    for variant, read_evidence in read_evidence_gen:
        # generate protein sequences by assembling variant reads
        protein_sequences = \
            protein_sequence_creator.sorted_protein_sequences_for_variant(
                variant=variant,
                read_evidence=read_evidence,
                transcript_id_whitelist=transcript_id_whitelist)

        isovar_result =  IsovarResult(
            variant=variant,
            transcript_id_whitelist=transcript_id_whitelist,
            read_evidence=read_evidence,
            protein_sequences=protein_sequences)
        isovar_result = isovar_result.clone_with_extra_filters(filter_thresholds)
        yield isovar_result

