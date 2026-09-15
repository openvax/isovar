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


from varcode import load_vcf
from pysam import AlignmentFile
from collections import OrderedDict

from .protein_sequence_creator import ProteinSequenceCreator
from .read_collector import ReadCollector
from .logging import get_logger
from .isovar_result import IsovarResult
from .default_parameters import (
    DEFAULT_FILTER_THRESHOLDS as DEFAULT_FILTER_THRESHOLDS,
    DEFAULT_FILTER_FLAGS as DEFAULT_FILTER_FLAGS,
    MIN_SHARED_FRAGMENTS_FOR_PHASING,
    NUM_RNA_DECOMPRESSION_THREADS,
)
from .effect_prediction import top_varcode_effect
from .filtering import apply_filters
from .phasing import annotate_phased_variants

logger = get_logger(__name__)


def run_isovar(
        variants,
        alignment_file,
        transcript_id_whitelist=None,
        read_collector=None,
        protein_sequence_creator=None,
        filter_thresholds=None,
        filter_flags=None,
        min_shared_fragments_for_phasing=MIN_SHARED_FRAGMENTS_FOR_PHASING,
        decompression_threads=NUM_RNA_DECOMPRESSION_THREADS):
    """
    This is the main entrypoint into the Isovar library, which collects
    RNA reads supporting variants and translates their coding sequence
    into amino acid sequences. Collects both the read evidence and
    protein sequences into IsovarResult objects. The values of any filters
    which are supplied in the filter_thresholds argument are attached to
    each IsovarResult's filter_values field.

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
        such as "not_has_mutant_protein_sequence_from_rna".

    decompression_threads : int
        Number of threads used by htslib to decompress BAM/CRAM
        files.

    Returns
    -------
    list of IsovarResult
        One per variant. The `protein_sequences` field will be empty
        if no sequences could be determined.
    """
    if filter_thresholds is None:
        filter_thresholds = OrderedDict(DEFAULT_FILTER_THRESHOLDS)

    if filter_flags is None:
        filter_flags = list(DEFAULT_FILTER_FLAGS)

    if isinstance(variants, str):
        variants = load_vcf(variants)

    if isinstance(alignment_file, str):
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
