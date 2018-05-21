# Copyright (c) 2016-2018. Mount Sinai School of Medicine
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
This module combines variant cDNA sequences collected from a BAM file with
the reading frames of annotated reference transcripts to create candidate
translations.
"""


from __future__ import print_function, division, absolute_import
import math


from .reference_context import reference_contexts_for_variant
from .variant_sequences import reads_to_variant_sequences
from .genetic_code import translate_cdna
from .variant_sequence_in_reading_frame import (
    match_variant_sequence_to_reference_context,
)
from .default_parameters import (
    MIN_TRANSCRIPT_PREFIX_LENGTH,
    MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
    INCLUDE_MISMATCHES_AFTER_VARIANT,
    PROTEIN_SEQUENCE_LENGTH,
    MIN_ALT_RNA_READS,
    MIN_VARIANT_SEQUENCE_COVERAGE,
    VARIANT_SEQUENCE_ASSEMBLY,
)
from .dataframe_builder import dataframe_from_generator
from .value_object import ValueObject
from .logging import get_logger

logger = get_logger(__name__)


class TranslationKey(ValueObject):
    """
    TranslationKey contains fields related to a translated protein sequence
    which should be used to combine multiple equivalent mutated amino acid
    sequences.
    """
    __slots__ = [
        # translated sequence of a variant sequence in the ORF established
        # by a reference context
        "amino_acids",
        # half-open interval coordinates for variant amino acids
        # in the translated sequence
        "variant_aa_interval_start",
        "variant_aa_interval_end",
        # did the amino acid sequence end due to a stop codon or did we
        # just run out of sequence context around the variant?
        "ends_with_stop_codon",
        # was the variant a frameshift relative to the reference sequence?
        "frameshift"
    ]


class Translation(TranslationKey):
    """
    Translated amino acid sequence of a VariantSequenceInReadingFrame for a
    particular ReferenceContext and VariantSequence.
    """
    __slots__ = [
        "untrimmed_variant_sequence",
        "reference_context",
        "variant_sequence_in_reading_frame"
    ]

    def __init__(
            self,
            amino_acids,
            variant_aa_interval_start,
            variant_aa_interval_end,
            ends_with_stop_codon,
            frameshift,
            untrimmed_variant_sequence,
            reference_context,
            variant_sequence_in_reading_frame):
        # TODO: get rid of untrimmed_variant_sequence by making
        # VariantSequenceInReadingFrame keep track of its inputs
        self.amino_acids = amino_acids
        self.variant_aa_interval_start = variant_aa_interval_start
        self.variant_aa_interval_end = variant_aa_interval_end
        self.ends_with_stop_codon = ends_with_stop_codon
        self.frameshift = frameshift
        # this variant sequence might differ from the one
        # in variant_sequence_in_reading_frame due to trimming
        # required to match the reference
        self.untrimmed_variant_sequence = untrimmed_variant_sequence
        self.reference_context = reference_context
        self.variant_sequence_in_reading_frame = variant_sequence_in_reading_frame

    @property
    def reads(self):
        """
        RNA reads which were used to construct the coding sequence
        from which we translated these amino acids.
        """
        return self.untrimmed_variant_sequence.reads

    @property
    def reference_cdna_sequence_before_variant(self):
        return (
            self.
            variant_sequence_in_reading_frame.
            reference_cdna_sequence_before_variant)

    @property
    def number_mismatches(self):
        """Only counting number of mismatches before the variant locus.
        """
        return self.number_mismatches_before_variant

    @property
    def number_mismatches_before_variant(self):
        return self.variant_sequence_in_reading_frame.number_mismatches_before_variant

    @property
    def number_mismatches_after_variant(self):
        return self.variant_sequence_in_reading_frame.number_mismatches_after_variant

    @property
    def cdna_sequence(self):
        return self.variant_sequence_in_reading_frame.cdna_sequence

    @property
    def offset_to_first_complete_codon(self):
        return self.variant_sequence_in_reading_frame.offset_to_first_complete_codon

    @property
    def variant_cdna_interval_start(self):
        return self.variant_sequence_in_reading_frame.variant_cdna_interval_start

    @property
    def variant_cdna_interval_end(self):
        return self.variant_sequence_in_reading_frame.variant_cdna_interval_end

    def as_translation_key(self):
        """
        Project Translation object or any other derived class into just a
        TranslationKey, which has fewer fields and can be used as a
        dictionary key.
        """
        return TranslationKey(**{
            name: getattr(self, name)
            for name in TranslationKey._fields})

    @classmethod
    def from_variant_sequence_and_reference_context(
            cls,
            variant_sequence,
            reference_context,
            min_transcript_prefix_length,
            max_transcript_mismatches,
            include_mismatches_after_variant,
            protein_sequence_length=None):
        """
        Attempt to translate a single VariantSequence using the reading frame
        from a single ReferenceContext.

        Parameters
        ----------
        variant_sequence : VariantSequence

        reference_context : ReferenceContext

        min_transcript_prefix_length : int
            Minimum number of nucleotides before the variant to test whether
            our variant sequence can use the reading frame from a reference
            transcript.

        max_transcript_mismatches : int
            Don't use the reading frame from a context where the cDNA variant
            sequences disagrees at more than this number of positions before the
            variant nucleotides.

        include_mismatches_after_variant : bool
            If true, mismatches after the variant nucleotides will also count
            against max_transcript_mismatches filtering.

        protein_sequence_length : int, optional
            Truncate protein to be at most this long

        Returns either a ProteinSequence object or None if the number of
        mismatches between the RNA and reference transcript sequences exceeds
        given threshold.
        """
        variant_sequence_in_reading_frame = match_variant_sequence_to_reference_context(
            variant_sequence,
            reference_context,
            min_transcript_prefix_length=min_transcript_prefix_length,
            max_transcript_mismatches=max_transcript_mismatches,
            include_mismatches_after_variant=include_mismatches_after_variant)

        if variant_sequence_in_reading_frame is None:
            logger.info("Unable to determine reading frame for %s", variant_sequence)
            return None

        cdna_sequence = variant_sequence_in_reading_frame.cdna_sequence
        cdna_codon_offset = variant_sequence_in_reading_frame.offset_to_first_complete_codon

        # get the offsets into the cDNA sequence which pick out the variant nucleotides
        cdna_variant_start_offset = variant_sequence_in_reading_frame.variant_cdna_interval_start
        cdna_variant_end_offset = variant_sequence_in_reading_frame.variant_cdna_interval_end

        # TODO: determine if the first codon is the start codon of a
        # transcript, for now any of the unusual start codons like CTG
        # will translate to leucine instead of methionine.
        variant_amino_acids, ends_with_stop_codon = translate_cdna(
            cdna_sequence[cdna_codon_offset:],
            first_codon_is_start=False,
            mitochondrial=reference_context.mitochondrial)

        variant_aa_interval_start, variant_aa_interval_end, frameshift = \
            find_mutant_amino_acid_interval(
                cdna_sequence=cdna_sequence,
                cdna_first_codon_offset=cdna_codon_offset,
                cdna_variant_start_offset=cdna_variant_start_offset,
                cdna_variant_end_offset=cdna_variant_end_offset,
                n_ref=len(reference_context.sequence_at_variant_locus),
                n_amino_acids=len(variant_amino_acids))

        if protein_sequence_length and len(variant_amino_acids) > protein_sequence_length:
            if protein_sequence_length <= variant_aa_interval_start:
                logger.warn(
                    ("Truncating amino acid sequence %s "
                     "to only %d elements loses all variant residues"),
                    variant_amino_acids,
                    protein_sequence_length)
                return None
            # if the protein is too long then shorten it, which implies
            # we're no longer stopping due to a stop codon and that the variant
            # amino acids might need a new stop index
            variant_amino_acids = variant_amino_acids[:protein_sequence_length]
            variant_aa_interval_end = min(variant_aa_interval_end, protein_sequence_length)
            ends_with_stop_codon = False

        return Translation(
            amino_acids=variant_amino_acids,
            frameshift=frameshift,
            ends_with_stop_codon=ends_with_stop_codon,
            variant_aa_interval_start=variant_aa_interval_start,
            variant_aa_interval_end=variant_aa_interval_end,
            untrimmed_variant_sequence=variant_sequence,
            reference_context=reference_context,
            variant_sequence_in_reading_frame=variant_sequence_in_reading_frame)


def find_mutant_amino_acid_interval(
        cdna_sequence,
        cdna_first_codon_offset,
        cdna_variant_start_offset,
        cdna_variant_end_offset,
        n_ref,
        n_amino_acids):
    """
    Parameters
    ----------
    cdna_sequence : skbio.DNA or str
        cDNA sequence found in RNAseq data

    cdna_first_codon_offset : int
        Offset into cDNA sequence to first complete codon, lets us skip
        past UTR region and incomplete codons.

    cdna_variant_start_offset : int
        Interbase start offset into cDNA sequence for selecting mutant
        nucleotides.

    cdna_variant_end_offset : int
        Interbase end offset into cDNA sequence for selecting mutant
        nucleotides.

    n_ref : int
        Number of reference nucleotides

    n_amino_acids : int
        Number of translated amino acids

    Returns tuple with three fields:
        1) Start offset for interval of mutant amino acids in translated sequence
        2) End offset for interval of mutant amino acids in translated sequence
        3) Boolean flag indicating whether the variant was a frameshift.
    """
    cdna_alt_nucleotides = cdna_sequence[
        cdna_variant_start_offset:cdna_variant_end_offset]

    n_alt = len(cdna_alt_nucleotides)

    # sequence of nucleotides before the variant starting from the first codon
    cdna_coding_prefix = cdna_sequence[cdna_first_codon_offset:cdna_variant_start_offset]

    # rounding down since a change in the middle of a codon should count
    # toward the variant codons
    n_coding_nucleotides_before_variant = len(cdna_coding_prefix)

    n_complete_prefix_codons = n_coding_nucleotides_before_variant // 3

    frame_of_variant_nucleotides = n_coding_nucleotides_before_variant % 3
    frameshift = abs(n_ref - n_alt) % 3 != 0
    indel = n_ref != n_alt

    variant_aa_interval_start = n_complete_prefix_codons

    if frameshift:
        # if mutation is a frame shift then every amino acid from the
        # first affected codon to the stop is considered mutant
        #
        # TODO: what if the first k amino acids are synonymous with the
        # reference sequence?
        variant_aa_interval_end = n_amino_acids
    else:
        n_alt_codons = int(math.ceil(n_alt / 3.0))
        if indel:
            # We need to adjust the number of affected codons by whether the
            # variant is aligned with codon boundaries, since in-frame indels
            # may still be split across multiple codons.
            #
            # Example of in-frame deletion of 3 nucleotides which leaves
            # 0 variant codons in the sequence (interval = 1:1)
            #   ref = CCC|AAA|GGG|TTT
            #   alt = CCC|GGG|TTT
            #
            # Example of in-frame deletion of 3 nucleotides which leaves
            # 1 variant codon in the sequence (interval = 1:2)
            #   ref = CCC|AAA|GGG|TTT
            #   alt = CCC|AGG|TTT
            #
            # Example of in-frame insertion of 3 nucleotides which
            # yields two variant codons:
            #   ref = CCC|AAA|GGG|TTT
            #   alt = CTT|TCC|AAA|GGG|TTT
            extra_affected_codon = int(frame_of_variant_nucleotides != 0)
            variant_aa_interval_end = (
                variant_aa_interval_start + n_alt_codons + extra_affected_codon)
        else:
            # if the variant is a simple substitution then it only affects
            # as many codons as are in the alternate sequence
            variant_aa_interval_end = variant_aa_interval_start + n_alt_codons
    return variant_aa_interval_start, variant_aa_interval_end, frameshift


def translation_generator(
        variant_sequences,
        reference_contexts,
        min_transcript_prefix_length,
        max_transcript_mismatches,
        include_mismatches_after_variant,
        protein_sequence_length=None):
    """
    Given all detected VariantSequence objects for a particular variant
    and all the ReferenceContext objects for that locus, translate
    multiple protein sequences, up to the number specified by the argument
    max_protein_sequences_per_variant.

    Parameters
    ----------
    variant_sequences : list of VariantSequence objects
        Variant sequences overlapping a single original variant

    reference_contexts : list of ReferenceContext objects
        Reference sequence contexts from the same variant as the variant_sequences

    min_transcript_prefix_length : int
        Minimum number of nucleotides before the variant to test whether
        our variant sequence can use the reading frame from a reference
        transcript.

    max_transcript_mismatches : int
        Maximum number of mismatches between coding sequence before variant
        and reference transcript we're considering for determing the reading
        frame.

    include_mismatches_after_variant : bool
        If true, mismatches occurring after the variant locus will also count
        toward max_transcript_mismatches filtering.

    protein_sequence_length : int, optional
        Truncate protein to be at most this long.

    Yields a sequence of Translation objects.
    """
    for reference_context in reference_contexts:
        for variant_sequence in variant_sequences:
            translation = Translation.from_variant_sequence_and_reference_context(
                variant_sequence=variant_sequence,
                reference_context=reference_context,
                min_transcript_prefix_length=min_transcript_prefix_length,
                max_transcript_mismatches=max_transcript_mismatches,
                include_mismatches_after_variant=include_mismatches_after_variant,
                protein_sequence_length=protein_sequence_length)
            if translation is not None:
                yield translation


def translate_variant_reads(
        variant,
        variant_reads,
        protein_sequence_length,
        transcript_id_whitelist=None,
        min_alt_rna_reads=MIN_ALT_RNA_READS,
        min_variant_sequence_coverage=MIN_VARIANT_SEQUENCE_COVERAGE,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
        max_transcript_mismatches=MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
        include_mismatches_after_variant=INCLUDE_MISMATCHES_AFTER_VARIANT,
        variant_sequence_assembly=VARIANT_SEQUENCE_ASSEMBLY):
    """
    Given a variant and its associated alt reads, construct variant sequences
    and translate them into Translation objects.

    Returns 0 or more Translation objects.

    Parameters
    ----------
    variant : varcode.Variant

    variant_reads : sequence or generator
        AlleleRead objects supporting the variant

    protein_sequence_length : int
        Try to translate protein sequences of this length, though sometimes
        we'll have to return something shorter (depending on the RNAseq data,
        and presence of stop codons).

    transcript_id_whitelist : set, optional
        If given, expected to be a set of transcript IDs which we should use
        for determining the reading frame around a variant. If omitted, then
        try to use all overlapping reference transcripts.

    min_alt_rna_reads : int
        Drop variant sequences from loci with fewer than this number of
        RNA reads supporting the alt allele.

    min_variant_sequence_coverage : int
        Trim variant sequences to nucleotides covered by at least this many
        reads.

    min_transcript_prefix_length : int
        Minimum number of bases we need to try matching between the reference
        context and variant sequence.

    max_transcript_mismatches : int
        Don't try to determine the reading frame for a transcript if more
        than this number of bases differ.

    include_mismatches_after_variant : bool
        Include mismatches after the variant locus in the count compared
        against max_transcript_mismatches.

    variant_sequence_assembly : bool
        Use overlap assembly to construct longer variant cDNA sequences.
    """
    if len(variant_reads) == 0:
        logger.info("No supporting reads for variant %s", variant)
        return []

    # Adding an extra codon to the desired RNA sequence length in case we
    # need to clip nucleotides at the start/end of the sequence
    cdna_sequence_length = (protein_sequence_length + 1) * 3

    variant_sequences = reads_to_variant_sequences(
        variant=variant,
        reads=variant_reads,
        preferred_sequence_length=cdna_sequence_length,
        min_alt_rna_reads=min_alt_rna_reads,
        min_variant_sequence_coverage=min_variant_sequence_coverage,
        variant_sequence_assembly=variant_sequence_assembly)

    if not variant_sequences:
        logger.info("No spanning cDNA sequences for variant %s", variant)
        return []

    # try translating the variant sequences from the same set of
    # ReferenceContext objects, which requires using the longest
    # context_size to be compatible with all of the sequences. Some
    # sequences maybe have fewer nucleotides than this before the variant
    # and will thus have to be trimmed.
    context_size = max(
        len(variant_sequence.prefix)
        for variant_sequence in variant_sequences)

    reference_contexts = reference_contexts_for_variant(
        variant,
        context_size=context_size,
        transcript_id_whitelist=transcript_id_whitelist)

    return list(translation_generator(
        variant_sequences=variant_sequences,
        reference_contexts=reference_contexts,
        min_transcript_prefix_length=min_transcript_prefix_length,
        max_transcript_mismatches=max_transcript_mismatches,
        include_mismatches_after_variant=include_mismatches_after_variant,
        protein_sequence_length=protein_sequence_length))


def translate_variants(
        variants_with_supporting_reads,
        transcript_id_whitelist=None,
        protein_sequence_length=PROTEIN_SEQUENCE_LENGTH,
        min_alt_rna_reads=MIN_ALT_RNA_READS,
        min_variant_sequence_coverage=MIN_VARIANT_SEQUENCE_COVERAGE,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
        max_transcript_mismatches=MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
        include_mismatches_after_variant=INCLUDE_MISMATCHES_AFTER_VARIANT,
        variant_sequence_assembly=VARIANT_SEQUENCE_ASSEMBLY):
    """
    Translates each coding variant in a collection to one or more protein
    fragment sequences (if the variant is not filtered and its spanning RNA
    sequences can be given a reading frame).

    Parameters
    ----------
    variants_with_reads : sequence or generator
        Each item of this sequence should be a pair containing a varcode.Variant
        and a list of AlleleRead objects supporting that variant.

    transcript_id_whitelist : set, optional
        If given, expected to be a set of transcript IDs which we should use
        for determining the reading frame around a variant. If omitted, then
        try to use all overlapping reference transcripts.

    protein_sequence_length : int
        Try to translate protein sequences of this length, though sometimes
        we'll have to return something shorter (depending on the RNAseq data,
        and presence of stop codons).

    min_alt_rna_reads : int
        Drop variant sequences from loci with fewer than this number of
        RNA reads supporting the alt allele.

    min_variant_sequence_coverage : int
        Trim variant sequences to nucleotides covered by at least this many
        reads.

    min_transcript_prefix_length : int
        Minimum number of bases we need to try matching between the reference
        context and variant sequence.

    max_transcript_mismatches : int
        Don't try to determine the reading frame for a transcript if more
        than this number of bases differ.

    include_mismatches_after_variant : bool
        Include mismatches after the variant locus in the count compared
        against max_transcript_mismatches.

    variant_sequence_assembly : bool
        Use overlap assembly to construct longer variant cDNA sequences.

    Yields pairs of a Variant and a sequence of all its candidate
    Translation objects.
    """
    for variant, variant_reads in variants_with_supporting_reads:
        translations = translate_variant_reads(
            variant=variant,
            variant_reads=variant_reads,
            protein_sequence_length=protein_sequence_length,
            transcript_id_whitelist=transcript_id_whitelist,
            min_alt_rna_reads=min_alt_rna_reads,
            min_variant_sequence_coverage=min_variant_sequence_coverage,
            min_transcript_prefix_length=min_transcript_prefix_length,
            max_transcript_mismatches=max_transcript_mismatches,
            include_mismatches_after_variant=include_mismatches_after_variant,
            variant_sequence_assembly=variant_sequence_assembly)
        yield variant, translations


def translations_generator_to_dataframe(translations_generator):
    """
    Given a generator of (Variant, [Translation]) pairs,
    returns a DataFrame of translated protein fragments with columns
    for each field of a Translation object (and chr/pos/ref/alt per variant).
    """
    return dataframe_from_generator(
        element_class=Translation,
        variant_and_elements_generator=translations_generator,
        exclude=[],
        converters={
            "untrimmed_variant_sequence": lambda vs: vs.sequence,
            "variant_sequence_in_reading_frame": (
                lambda vs: vs.in_frame_cdna_sequence),
            "reference_context": (
                lambda rc: ";".join([
                    transcript.name for
                    transcript in rc.transcripts]))
        },
        extra_column_fns={
            "untrimmed_variant_sequence_read_count": (
                lambda _, t: len(t.untrimmed_variant_sequence.reads)),
        })
