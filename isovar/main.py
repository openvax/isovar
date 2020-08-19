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
from varcode import load_vcf
from pysam import AlignmentFile
from collections import OrderedDict

from .protein_sequence_creator import ProteinSequenceCreator
from .read_collector import ReadCollector
from .logging import get_logger
from .isovar_result import IsovarResult
from .default_parameters import (
    MIN_NUM_RNA_ALT_READS,
    MIN_NUM_RNA_ALT_FRAGMENTS,
    MIN_FRACTION_RNA_ALT_READS,
    MIN_FRACTION_RNA_ALT_FRAGMENTS,
    MAX_NUM_RNA_REF_READS,
    MAX_NUM_RNA_REF_FRAGMENTS,
    MAX_FRACTION_RNA_REF_READS,
    MAX_FRACTION_RNA_REF_FRAGMENTS,
    MAX_NUM_RNA_OTHER_READS,
    MAX_NUM_RNA_OTHER_FRAGMENTS,
    MAX_FRACTION_RNA_OTHER_READS,
    MAX_FRACTION_RNA_OTHER_FRAGMENTS,
    MIN_RATIO_RNA_ALT_TO_OTHER_FRAGMENTS,
    MIN_SHARED_FRAGMENTS_FOR_PHASING
)
from .effect_prediction import top_varcode_effect
from .filtering import apply_filters
from .phasing import annotate_phased_variants

logger = get_logger(__name__)


DEFAULT_FILTER_THRESHOLDS =  OrderedDict([
    # alt allele
    ("min_num_alt_reads", MIN_NUM_RNA_ALT_READS),
    ("min_num_alt_fragments", MIN_NUM_RNA_ALT_FRAGMENTS),
    ("min_fraction_alt_reads", MIN_FRACTION_RNA_ALT_READS),
    ("min_fraction_alt_fragments", MIN_FRACTION_RNA_ALT_FRAGMENTS),

    # ref allele coverage and VAF
    ("max_num_ref_reads", MAX_NUM_RNA_REF_READS),
    ("max_num_ref_fragments", MAX_NUM_RNA_REF_FRAGMENTS),
    ("max_fraction_ref_reads", MAX_FRACTION_RNA_REF_READS),
    ("max_fraction_ref_fragments", MAX_FRACTION_RNA_REF_FRAGMENTS),

    # other alleles
    ("max_num_other_reads", MAX_NUM_RNA_OTHER_READS),
    ("max_num_other_fragments", MAX_NUM_RNA_OTHER_FRAGMENTS),
    ("max_fraction_other_reads", MAX_FRACTION_RNA_OTHER_READS),
    ("max_fraction_other_fragments", MAX_FRACTION_RNA_OTHER_FRAGMENTS),

    # misc. filters
    ("min_ratio_alt_to_other_fragments", MIN_RATIO_RNA_ALT_TO_OTHER_FRAGMENTS)
])

DEFAULT_FILTER_FLAGS = [
    "predicted_effect_modifies_protein_sequence",
    "has_mutant_protein_sequence_from_rna",
    "protein_sequence_contains_mutation",
]

def run_isovar(
        variants,
        alignment_file,
        transcript_id_whitelist=None,
        read_collector=None,
        protein_sequence_creator=None,
        filter_thresholds=DEFAULT_FILTER_THRESHOLDS,
        filter_flags=DEFAULT_FILTER_FLAGS,
        min_shared_fragments_for_phasing=MIN_SHARED_FRAGMENTS_FOR_PHASING,
        decompression_threads=1):
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

    filter_flags : list of str
        List of boolean fields of IsovarResult used for filtering,
        they can also be negated by prepending "not_",
        such as "not_has_protein_sequence".

    decompress_threads : int
        Number of threads used by htslib to decompress BAM/CRAM
        files.

    Generator of IsovarResult objects, one for each variant. The
    `protein_sequences` field of the IsovarVar result will be empty
    if no sequences could be determined.
    """
    if isinstance(variants, string_types):
        variants = load_vcf(variants)

    if isinstance(alignment_file, string_types):
        alignment_file = AlignmentFile(
            alignment_file,
            threads=decompression_threads)

    if read_collector is None:
        read_collector = ReadCollector()

    if protein_sequence_creator is None:
        protein_sequence_creator = ProteinSequenceCreator()

    # create generator which returns (Variant, ReadEvidence) pairs
    read_evidence_gen = \
        read_collector.read_evidence_generator(
           variants=variants,
           alignment_file=alignment_file)

    results = []
    for variant, read_evidence in read_evidence_gen:
        # generate protein sequences by assembling variant reads
        protein_sequences = \
            protein_sequence_creator.sorted_protein_sequences_for_variant(
                variant=variant,
                read_evidence=read_evidence,
                transcript_id_whitelist=transcript_id_whitelist)
        predicted_effect = top_varcode_effect(
            variant=variant,
            transcript_id_whitelist=transcript_id_whitelist)
        isovar_result = IsovarResult(
            variant=variant,
            predicted_effect=predicted_effect,
            read_evidence=read_evidence,
            sorted_protein_sequences=protein_sequences)
        isovar_result = apply_filters(
            isovar_result,
            filter_thresholds=filter_thresholds,
            filter_flags=filter_flags)
        results.append(isovar_result)
    results = annotate_phased_variants(
        results,
        min_shared_fragments_for_phasing)
    return results
