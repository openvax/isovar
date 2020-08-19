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
IsovarResult is a collection of all information gathered about a variant
and any protein sequences which were successfully translated for it.
"""

from __future__ import print_function, division, absolute_import

from collections import OrderedDict

from cached_property import cached_property

from .common import safediv
from .alignment_score import alignment_score

class IsovarResult(object):
    """
    This object represents all information gathered about a variant,
    which includes the AlleleReads supporting any allele at this variant's
    locus and any protein sequences generated from an alt-allele cDNA
    assembly.
    """

    def __init__(
            self,
            variant,
            read_evidence,
            predicted_effect,
            sorted_protein_sequences=None,
            filter_values=None,
            phased_variants_in_supporting_reads=None,
            phased_variants_in_protein_sequence=None):
        """
        Parameters
        ----------
        variant : varcode.Variant
            Somatic mutation.

        read_evidence : ReadEvidence
            Reads grouped by allele at the locus of this variant.

        predicted_effect : varcode effect
            Top predicted impact on protein sequence from DNA.

        sorted_protein_sequences : list of ProteinSequence
            Protein sequences sorted by degree of support.

        filter_values : OrderedDict
            Dictionary from names of filters to whether this IsovarResult
            passed that filter. 

        phased_variants_in_supporting_reads : set of varcode.Variant
            Other somatic variants which occur in the alt reads supporting the 
            variant associated with this IsovarResult.

        phased_variants_in_protein_sequence : set of varcode.Variant
            Other somatic variants which occur in the reads used to construct
            the top protein sequence associated with this IsovarResult.
        """
        self.variant = variant
        self.read_evidence = read_evidence
        self.predicted_effect = predicted_effect

        if sorted_protein_sequences is None:
            sorted_protein_sequences = []

        self.sorted_protein_sequences = sorted_protein_sequences

        if filter_values is None:
            self.filter_values = OrderedDict()
        else:
            self.filter_values = filter_values

        if phased_variants_in_supporting_reads is None:
            self.phased_variants_in_supporting_reads = set()
        else:
            self.phased_variants_in_supporting_reads = \
                phased_variants_in_supporting_reads

        if phased_variants_in_protein_sequence is None:
            self.phased_variants_in_protein_sequence = set()
        else:
            self.phased_variants_in_protein_sequence = \
                phased_variants_in_protein_sequence

    @property
    def fields(self):
        """
        List of field names used to construct an IsovarResult instance.
        """
        return [
            "variant",
            "predicted_effect",
            "read_evidence",
            "sorted_protein_sequences",
            "filter_values",
            "phased_variants_in_supporting_reads",
            "phased_variants_in_protein_sequence",
        ]

    def __str__(self):
        field_strings = ["%s=%s" % (k, v) for (k, v) in self.to_dict().items()]
        return "%s(%s)" % (
            self.__class__.__name__,
            ", ".join(field_strings)
        )

    def __eq__(self, other):
        for field_name in self.fields:
            my_value = getattr(self, field_name)
            other_value = getattr(other, field_name)
            if my_value is None and other_value is not None:
                return False
            elif other_value is None and my_value is not None:
                return False
            elif my_value != other_value:
                return False
        return True

    def __repr__(self):
        return str(self)

    def to_dict(self):
        """
        Dictionary representation of fields used to construct this IsovarResult

        Returns dict
        """
        return OrderedDict([
            (k, getattr(self, k))
            for k in self.fields
        ])

    def clone_with_updates(self, **kwargs):
        """
        Create a copy of this IsovarResult object including any new
        parameters in `kwargs`.

        Returns IsovarResult
        """
        for (k, v) in self.to_dict().items():
            if k not in kwargs:
                kwargs[k] = v
        return IsovarResult(**kwargs)

    def clone(self):
        """
        Make a copy of this IsovarResult without updating any fields.

        Returns IsovarResult
        """
        return self.clone_with_updates()

    def to_record(self):
        """
        Create an OrderedDict of essential information from
        this IsovarResult to be used for building a DataFrame across
        variants.

        Returns OrderedDict
        """
        d = OrderedDict([
            ("variant", self.variant.short_description),
            ("overlapping_gene_names",
                ";".join(self.overlapping_gene_names(only_coding=False))),
            ("overlapping_gene_ids",
                ";".join(self.overlapping_gene_ids(only_coding=False))),
            ("overlapping_coding_gene_names",
             ";".join(self.overlapping_gene_names(only_coding=True))),
            ("overlapping_coding_gene_ids",
             ";".join(self.overlapping_gene_ids(only_coding=True))),

        ])

        # get all quantitative fields from this object
        for key in dir(self):
            if key.startswith("num_") or key.startswith("fraction_") or key.startswith("ratio_"):
                d[key] = getattr(self, key)

        # get all boolean properties that start with "has_"
        for key in dir(self):
            if key.startswith("has_"):
                d[key] = getattr(self, key)

        ########################################################################
        # predicted protein changes without looking at RNA reads
        ########################################################################
        effect = self.predicted_effect

        d["predicted_effect"] = effect.short_description
        d["predicted_effect_class"] = effect.__class__.__name__
        d["predicted_effect_modifies_protein_sequence"] = \
            self.predicted_effect_modifies_protein_sequence

        # list of field names on varcode effect properties
        effect_properties = [
            "gene_name",
            "gene_id",
            "transcript_id",
            "transcript_name",
            "modifies_protein_sequence",
            "original_protein_sequence",
            "aa_mutation_start_offset",
            "aa_mutation_end_offset",
            "mutant_protein_sequence"
        ]
        for field_name in effect_properties:
            # store effect fields with prefix 'predicted_effect_' and use
            # getattr in case the field is not available for all effects
            d["predicted_effect_%s" % field_name] = getattr(
                effect,
                field_name,
                None)

        ########################################################################
        # get the top protein sequence, if one exists
        ########################################################################
        protein_sequence = self.top_protein_sequence

        # list of names we want to use in the result dictionary,
        # paired with names of fields on ProteinSequence
        protein_sequence_properties = [
            ("protein_sequence", "amino_acids"),
            ("protein_sequence_ends_with_stop_codon", "ends_with_stop_codon"),
            ("protein_sequence_gene_names", "gene_names"),
            ("protein_sequence_gene_ids", "gene_ids"),
            ("protein_sequence_transcript_names", "transcript_names"),
            ("protein_sequence_transcript_ids", "transcript_ids"),
        ]
        for (name, protein_sequence_field) in protein_sequence_properties:
            value = getattr(protein_sequence, protein_sequence_field, None)
            if isinstance(value, (list, set, tuple)):
                value = ";".join(value)
            d[name] = value
        d["protein_sequence_mutation_start_idx"] = \
            self.protein_sequence_mutation_start_idx
        d["protein_sequence_mutation_end_idx"] = \
            self.protein_sequence_mutation_end_idx
        d["num_mutant_amino_acids_in_protein_sequence"] = \
            self.num_mutant_amino_acids_in_protein_sequence

        d["trimmed_predicted_mutant_protein_sequence"] = self.trimmed_predicted_mutant_protein_sequence
        d["trimmed_reference_protein_sequence"] = self.trimmed_reference_protein_sequence
        d["protein_sequence_contains_mutation"] = self.protein_sequence_contains_mutation
        d["protein_sequence_contains_deletion"] = self.protein_sequence_contains_deletion
        d["protein_sequence_matches_predicted_mutation_effect"] = \
            self.protein_sequence_matches_predicted_mutation_effect

        ########################################################################
        # filters
        ########################################################################
        for filter_name, filter_value in self.filter_values.items():
            d["filter:%s" % filter_name] = filter_value
        d["passes_all_filters"] = self.passes_all_filters
        return d


    @cached_property
    def passes_all_filters(self):
        """
        Does this IsovarResult have True for all the filter values in
        self.filter_values?
        """
        if len(self.filter_values) == 0:
            return True
        else:
            return all(list(self.filter_values.values()))

    @cached_property
    def top_gene_name(self):
        """
        If there is a single best gene name associated with this IsovarResult's
        mutation, then return that gene name. Otherwise return all overlapping
        gene names.

        Returns
        -------
        str
        """
        if self.top_protein_sequence is not None:
            return self.top_protein_sequence.gene_name

        # if no protein sequence was determined from the RNA,
        # see if the top predicted effect is associated with a gene
        gene_name = getattr(self.predicted_effect, "gene_name", None)
        if gene_name:
            return gene_name
        else:
            # if nothing about the gene can be guessed from the predicted
            # effect, then use all overlapping gene names
            return ";".join(
                self.overlapping_gene_names(only_coding=False))

    @cached_property
    def top_protein_sequence(self):
        """
        If any protein sequences were assembled for this variant then
        return the best according to coverage, number of mismatches
        relative to the reference, number of reference transcripts
        which match sequence before the variant and protein
        sequence length.

        Returns ProteinSequence or None
        """
        if len(self.sorted_protein_sequences) > 0:
            return self.sorted_protein_sequences[0]
        else:
            return None

    @cached_property
    def protein_sequence_mutation_start_idx(self):
        """
        Interbase start coordinate for mutated amino acids in top protein
        sequence.

        Returns
        -------
        int or None
        """
        if self.has_mutant_protein_sequence_from_rna:
            return self.top_protein_sequence.mutation_start_idx
        else:
            return None

    @cached_property
    def protein_sequence_mutation_end_idx(self):
        """
        Interbase end coordinate for mutated amino acids in top protein
        sequence.

        Returns
        -------
        int or None
        """
        if self.has_mutant_protein_sequence_from_rna:
            return self.top_protein_sequence.mutation_end_idx
        else:
            return None

    @cached_property
    def predicted_effect_modifies_protein_sequence(self):
        """
        Does the predicted effect change the protein sequence?

        Returns bool
        """
        if self.predicted_effect is None:
            return None
        return self.predicted_effect.modifies_protein_sequence

    @cached_property
    def num_cdna_mismatches_in_top_protein_sequence(self):
        """
        How many nucleotide positions aside from the variant location don't
        match the reference transcript sequence in the cDNA sequence used
        to translate the top protein sequence. When multiple cDNA sequences
        were used then take the min across them.

        Returns int
        """
        if self.top_protein_sequence is None:
            return None
        return self.top_protein_sequence.num_mismatches

    @cached_property
    def num_cdna_mismatches_before_variant_in_top_protein_sequence(self):
        """
        How many nucleotide positions before the variant location don't
        match the reference transcript sequence in the cDNA sequence used
        to translate the top protein sequence. When multiple cDNA sequences
        were used then take the min across them.

        Returns int
        """
        if self.top_protein_sequence is None:
            return None
        return self.top_protein_sequence.num_mismatches_before_variant

    @cached_property
    def num_cdna_mismatches_after_variant_in_top_protein_sequence(self):
        """
        How many nucleotide positions after the variant location don't
        match the reference transcript sequence in the cDNA sequence used
        to translate the top protein sequence. When multiple cDNA sequences
        were used then take the min across them.

        Returns int
        """
        if self.top_protein_sequence is None:
            return None
        return self.top_protein_sequence.num_mismatches_after_variant

    @cached_property
    def num_reads_supporting_top_protein_sequence(self):
        """
        How many reads support the assembled protein sequence?

        Returns int
        """
        if self.top_protein_sequence is None:
            return 0
        return self.top_protein_sequence.num_supporting_reads

    @cached_property
    def num_fragments_supporting_top_protein_sequence(self):
        """
        How many fragments (unique read names) support the assembled
        protein sequence?

        Returns int
        """
        if self.top_protein_sequence is None:
            return 0
        return self.top_protein_sequence.num_supporting_fragments

    @cached_property
    def num_translations_from_top_protein_sequence(self):
        """
        How many distinct translations were used to create the top
        protein sequence? These can arise either from different cDNA sequences
        or different reading frames.

        Returns int
        """
        if self.top_protein_sequence is None:
            return 0
        return len(self.top_protein_sequence.translations)

    @cached_property
    def num_cdna_sequences_from_top_protein_sequence(self):
        """
        How many distinct cDNA sequences were used to create the top
        protein sequence?

        Returns int
        """
        if self.top_protein_sequence is None:
            return 0
        return self.top_protein_sequence.num_cdna_sequences

    @cached_property
    def trimmed_predicted_mutant_protein_sequence(self):
        """
        Trim the predicted mutant protein sequence from Varcode
        to match the length of the protein subsequence assembled from RNA.

        Returns str or None
        """
        p = self.top_protein_sequence
        e = self.predicted_effect

        if e is None or p is None:
            return None
        if e.mutant_protein_sequence is None:
            return None
        if e.aa_mutation_start_offset is None:
            return None
        n_before_mutation = p.mutation_start_idx
        n_after_mutation = len(p.amino_acids) - p.mutation_end_idx
        return e.mutant_protein_sequence[
             e.aa_mutation_start_offset - n_before_mutation:
             e.aa_mutation_start_offset + len(e.aa_alt) + n_after_mutation]

    @cached_property
    def has_mutant_protein_sequence_from_rna(self):
        """
        Does this IsovarResult have an associated protein sequence
        determined from RNA reads?

        Returns bool
        """
        return self.top_protein_sequence is not None

    @cached_property
    def has_reference_protein_sequence_from_predicted_effect(self):
        """
        Does the variant for this IsovarResult affect a protein
        whose reference sequence can be determined?

        Returns bool
        """
        if self.predicted_effect is None:
            return False
        return self.predicted_effect.original_protein_sequence is not None

    @cached_property
    def has_mutant_protein_sequence_from_predicted_effect(self):
        """
        Does the variant for this IsovarResult have a predicted
        mutant protein sequence?

        Returns bool
        """
        if self.predicted_effect is None:
            return False
        return self.predicted_effect.mutant_protein_sequence is not None

    @cached_property
    def trimmed_reference_protein_sequence(self):
        """
        Trim the reference protein sequence from the top Varcode effect
        to match the length of the protein sequence assembled from RNA.

        Returns str
        """
        p = self.top_protein_sequence
        e = self.predicted_effect
        if e is None or p is None:
            return None
        if e.original_protein_sequence is None:
            return None
        if e.aa_mutation_start_offset is None:
            return None
        n_before_mutation = p.mutation_start_idx
        n_total = len(p.amino_acids)
        start_index_in_original_protein = (
               e.aa_mutation_start_offset - n_before_mutation)
        return e.original_protein_sequence[
               start_index_in_original_protein:
               start_index_in_original_protein + n_total]

    @cached_property
    def num_amino_acid_mismatches_from_predicted_effect(self):
        """
        Compute the number of mismatches between the mutant protein sequence
        predicted by Varcode and the best supported sequence translated
        from assembled RNA reads by Isovar. We're not allowing any
        insertions or deletions in the middle of the sequences but do
        allow a shorter sequence to start anywhere within a longer one.

        Returns int
        """
        predicted_sequence = self.trimmed_predicted_mutant_protein_sequence
        if predicted_sequence is None:
            return None

        protein_sequence_object = self.top_protein_sequence
        if protein_sequence_object is None:
            return None

        assembled_protein_sequence = protein_sequence_object.amino_acids

        return alignment_score(
            assembled_protein_sequence,
            predicted_sequence)

    @cached_property
    def protein_sequence_matches_predicted_mutation_effect(self):
        """
        Does the top protein sequence translated from RNA reads
        match the predicted protein change determined by Varcode?

        Returns bool
        """
        n_mismatches = self.num_amino_acid_mismatches_from_predicted_effect
        if n_mismatches is None:
            return None
        return n_mismatches == 0


    @cached_property
    def num_amino_acid_mismatches_from_reference(self):
        """
        Compute the number of mismatches between the original protein sequence
        and the best supported sequence translated from assembled RNA reads by
        Isovar. We're not allowing any insertions or deletions in the middle of
        the sequences but do allow a shorter sequence to start anywhere
        within a longer one.

        Returns int
        """
        original_sequence = self.trimmed_reference_protein_sequence
        if original_sequence is None:
            return None

        protein_sequence_object = self.top_protein_sequence
        if protein_sequence_object is None:
            return None
        assembled_protein_sequence = protein_sequence_object.amino_acids

        return alignment_score(
            assembled_protein_sequence,
            original_sequence)

    @cached_property
    def protein_sequence_contains_mutation(self):
        """
        Does the protein sequence assembled from RNA contain a mutation?

        Returns bool
        """
        if self.top_protein_sequence is None:
            return None
        return self.top_protein_sequence.contains_mutation

    @cached_property
    def protein_sequence_contains_deletion(self):
        """
        Does the protein sequence assembled from RNA span both sides of a
        deletion mutation?

        Returns bool
        """
        if self.top_protein_sequence is None:
            return None
        return self.top_protein_sequence.contains_deletion

    @cached_property
    def num_mutant_amino_acids_in_protein_sequence(self):
        """
        Number of amino acids modified by mutation in the top protein sequence
        from RNA.

        Returns
        -------
        int or None
        """
        if self.top_protein_sequence is None:
            return None
        return self.top_protein_sequence.num_mutant_amino_acids

    @cached_property
    def protein_sequence_matches_reference(self):
        """
        Does the top protein sequence translated from RNA reads
        match the reference protein sequence?

        Returns bool
        """
        n_mismatches = self.num_amino_acid_mismatches_from_reference
        if n_mismatches is None:
            return None
        return n_mismatches == 0

    @cached_property
    def protein_sequence_different_from_reference(self):
        """
        Does the top protein sequence translated from RNA reads
        differ from the reference?

        Returns bool
        """
        if self.protein_sequence_matches_reference is None:
            return None
        return not self.protein_sequence_matches_reference

    @cached_property
    def num_protein_sequences(self):
        """
        Number of distinct protein sequences which were translated from
        assembled RNA reads.

        Returns int
        """
        return len(self.sorted_protein_sequences)

    def transcripts_from_protein_sequences(self, max_num_protein_sequences=None):
        """
        Ensembl transcript IDs of all transcripts which support the reading
        frame used by protein sequences in this IsovarResult.

        Parameters
        ----------
        max_num_protein_sequences : int or None
            If supplied then only consider the top protein sequences up to
            this number.

        Returns list of pyensembl.Transcript
        """
        transcript_set = set([])
        for p in self.sorted_protein_sequences[:max_num_protein_sequences]:
            transcript_set.update(p.transcripts)
        return sorted(transcript_set)

    def transcript_ids_from_protein_sequences(self, max_num_protein_sequences=None):
        """
        Ensembl transcripts IDs which support the reading frame used by protein
        sequences in this IsovarResult.

        Parameters
        ----------
        max_num_protein_sequences : int or None
            If supplied then only consider the top protein sequences up to
            this number.

        Returns list of str
        """
        return sorted({t.id for t in self.transcripts_from_protein_sequences(
            max_num_protein_sequences=max_num_protein_sequences)})

    @cached_property
    def num_transcripts_from_protein_sequences(self):
        """
        Number of genes used by any translated protein sequence associated
        with this IsovarResult.

        Returns int
        """
        return len(self.transcript_ids_from_protein_sequences(
            max_num_protein_sequences=None))

    @cached_property
    def num_transcripts_from_top_protein_sequence(self):
        """
        Number of genes used by any translated protein sequence associated
        with this IsovarResult.

        Returns int
        """
        return len(self.transcript_ids_from_protein_sequences(
            max_num_protein_sequences=1))

    def genes_from_protein_sequences(self, max_num_protein_sequences=None):
        """
        Ensembl genes which support the reading frame used by protein
        sequences in this IsovarResult.

        Parameters
        ----------
        max_num_protein_sequences : int or None
            If supplied then only consider the top protein sequences up to
            this number.

        Returns list of pyensembl.Gene
        """
        transcripts = self.transcripts_from_protein_sequences(
            max_num_protein_sequences=max_num_protein_sequences)
        genes = [t.gene for t in transcripts]
        return sorted(genes)

    def gene_ids_from_protein_sequences(self, max_num_protein_sequences=None):
        """
        Ensembl genes IDs which support the reading frame used by protein
        sequences in this IsovarResult.

        Parameters
        ----------
        max_num_protein_sequences : int or None
            If supplied then only consider the top protein sequences up to
            this number.

        Returns list of str
        """
        return sorted({
            g.id
            for g
            in
            self.genes_from_protein_sequences(
                max_num_protein_sequences=max_num_protein_sequences)
        })

    @cached_property
    def num_genes_from_protein_sequences(self):
        """
        Number of genes used by any translated protein sequence associated
        with this IsovarResult.

        Returns int
        """
        return len(self.gene_ids_from_protein_sequences(
            max_num_protein_sequences=None))

    @cached_property
    def num_genes_from_top_protein_sequence(self):
        """
        Number of genes used by any translated protein sequence associated
        with this IsovarResult.

        Returns int
        """
        return len(self.gene_ids_from_protein_sequences(
            max_num_protein_sequences=1))

    def overlapping_transcripts(self, only_coding=True):
        """
        Transcripts which this variant overlaps.

        Parameters
        ----------
        only_coding : bool
            Only return transcripts which are annotated as coding for a
            protein (default=True)

        Returns set of pyensembl.Transcript objects
        """
        return {
            t
            for t in self.variant.transcripts
            if not only_coding or t.is_protein_coding
        }

    def overlapping_transcript_ids(self, only_coding=True):
        """
        Transcript IDs which this variant overlaps.

        Parameters
        ----------
        only_coding : bool
            Only return transcripts which are annotated as coding for a
            protein (default=True)
        Returns set of str
        """
        return {
            t.id
            for t in self.variant.transcripts
            if not only_coding or t.is_protein_coding
        }

    @cached_property
    def num_overlapping_transcripts(self):
        """
        Number of transcripts overlapped by the variant

        Returns int
        """
        return len(self.overlapping_transcript_ids(only_coding=False))

    @cached_property
    def num_overlapping_coding_transcripts(self):
        """
        Number of coding transcripts overlapped by the variant

        Returns int
        """
        return len(self.overlapping_transcript_ids(only_coding=True))

    def overlapping_genes(self, only_coding=True):
        """
        Genes which this variant overlaps.

        Parameters
        ----------
        only_coding : bool
            Only return genes which are annotated as coding for a
            protein (default=True)

        Returns list of pyensembl.Gene objects
        """
        return sorted({
            g
            for g in self.variant.genes
            if not only_coding or g.is_protein_coding
        })

    def overlapping_gene_names(self, only_coding=True):
        """
        Names of genes which this variant overlaps.

        Parameters
        ----------
        only_coding : bool
           Only return genes which are annotated as coding for a
           protein (default=True)

        Returns list of str
        """
        return [
            g.name for g in self.overlapping_genes(only_coding=only_coding)
        ]

    def overlapping_gene_ids(self, only_coding=True):
        """
        Gene IDs which this variant overlaps.

        Parameters
        ----------
        only_coding : bool
            Only return genes which are annotated as coding for a
            protein (default=True)

        Returns set of str
        """
        return {
            g.id
            for g in self.variant.genes
            if not only_coding or g.is_protein_coding
        }

    @cached_property
    def num_overlapping_genes(self):
        """
        Number of genes overlapped by the variant

        Returns int
        """
        return len(self.overlapping_gene_ids(only_coding=False))

    @cached_property
    def num_overlapping_coding_genes(self):
        """
        Number of coding genes overlapped by the variant

        Returns int
        """
        return len(self.overlapping_gene_ids(only_coding=True))

    @cached_property
    def ref_reads(self):
        """
        AlleleRead objects at this locus which support the reference allele
        """
        return self.read_evidence.ref_reads

    @cached_property
    def alt_reads(self):
        """
        AlleleRead objects at this locus which support the mutant allele
        """
        return self.read_evidence.alt_reads

    @cached_property
    def other_reads(self):
        """
        AlleleRead objects at this locus which support some allele other than
        either the reference or alternate.
        """
        return self.read_evidence.other_reads

    @cached_property
    def ref_read_names(self):
        """
        Names of reference reads at this locus.
        """
        return {r.name for r in self.ref_reads}

    @cached_property
    def alt_read_names(self):
        """
        Names of alt reads at this locus.
        """
        return {r.name for r in self.alt_reads}

    @cached_property
    def ref_and_alt_read_names(self):
        """
        Names of reads which support either the ref or alt alleles.
        """
        return self.ref_read_names.union(self.alt_read_names)

    @cached_property
    def other_read_names(self):
        """
        Names of other (non-alt, non-ref) reads at this locus.
        """
        return {r.name for r in self.other_reads}

    @cached_property
    def all_read_names(self):
        """
        Names of all reads at this locus.
        """
        return self.ref_read_names.union(self.alt_read_names).union(self.other_read_names)

    @cached_property
    def num_total_reads(self):
        """
        Total number of reads at this locus, regardless of allele.
        """
        return self.num_ref_reads + self.num_alt_reads + self.num_other_reads

    @cached_property
    def num_total_fragments(self):
        """
        Total number of distinct fragments at this locus, which also corresponds
        to the total number of read names.
        """
        return len(self.all_read_names)

    @cached_property
    def num_ref_reads(self):
        """
        Number of reads which support the reference allele.
        """
        return len(self.ref_reads)

    @cached_property
    def num_ref_fragments(self):
        """
        Number of distinct fragments which support the reference allele.
        """
        return len(self.ref_read_names)

    @cached_property
    def num_alt_reads(self):
        """
        Number of reads which support the alt allele.
        """
        return len(self.alt_reads)

    @cached_property
    def num_alt_fragments(self):
        """
        Number of distinct fragments which support the alt allele.
        """
        return len(self.alt_read_names)

    @cached_property
    def num_other_reads(self):
        """
        Number of reads which support neither the reference nor alt alleles.
        """
        return len(self.other_reads)

    @cached_property
    def num_other_fragments(self):
        """
        Number of distinct fragments which support neither the reference nor
        alt alleles.
        """
        return len(self.other_read_names)

    @cached_property
    def fraction_ref_reads(self):
        """
        Allelic fraction of the reference allele among all reads at this site.
        """
        return safediv(self.num_ref_reads, self.num_total_reads)

    @cached_property
    def fraction_ref_fragments(self):
        """
        Allelic fraction of the reference allele among all fragments at this site.
        """
        return safediv(self.num_ref_fragments, self.num_total_fragments)

    @cached_property
    def fraction_alt_reads(self):
        """
        Allelic fraction of the variant allele among all reads at this site.
        """
        return safediv(self.num_alt_reads, self.num_total_reads)

    @cached_property
    def fraction_alt_fragments(self):
        """
        Allelic fraction of the variant allele among all fragments at this site.
        """
        return safediv(self.num_alt_fragments, self.num_total_fragments)

    @cached_property
    def fraction_other_reads(self):
        """
        Allelic fraction of the "other" (non-ref, non-alt) alleles among all
        reads at this site.
        """
        return safediv(self.num_other_reads, self.num_total_reads)

    @cached_property
    def fraction_other_fragments(self):
        """
        Allelic fraction of the "other" (non-ref, non-alt) alleles among all
        fragments at this site.
        """
        return safediv(self.num_other_fragments, self.num_total_fragments)

    @cached_property
    def ratio_other_to_ref_reads(self):
        """
        Ratio of the number of reads which support alleles which are neither
        ref/alt to the number of ref reads.
        """
        return safediv(self.num_other_reads, self.num_ref_reads)

    @cached_property
    def ratio_other_to_ref_fragments(self):
        """
        Ratio of the number of fragments which support alleles which are neither
        ref/alt to the number of ref fragments.
        """
        return safediv(self.num_other_fragments, self.num_ref_fragments)

    @cached_property
    def ratio_other_to_alt_reads(self):
        """
        Ratio of the number of reads which support alleles which are neither
        ref/alt to the number of alt reads.
        """
        return safediv(self.num_other_reads, self.num_alt_reads)

    @cached_property
    def ratio_other_to_alt_fragments(self):
        """
        Ratio of the number of fragments which support alleles which are neither
        ref/alt to the number of alt fragments.
        """
        return safediv(self.num_other_fragments, self.num_alt_fragments)

    @cached_property
    def ratio_ref_to_other_reads(self):
        """
        Ratio of the number of reference reads to non-ref/non-alt reads
        """
        return safediv(self.num_ref_reads, self.num_other_reads)

    @cached_property
    def ratio_ref_to_other_fragments(self):
        """
        Ratio of the number of reference fragments to non-ref/non-alt fragments
        """
        return safediv(self.num_ref_fragments, self.num_other_fragments)

    @cached_property
    def ratio_alt_to_other_reads(self):
        """
        Ratio of alt allele reads to non-ref/non-alt reads
        """
        return safediv(self.num_alt_reads, self.num_other_reads)

    @cached_property
    def ratio_alt_to_other_fragments(self):
        """
        Ratio of the number of fragments which support the alt allele
        to the number of non-alt/non-ref allele fragments.
        """
        return safediv(self.num_alt_fragments, self.num_other_fragments)


    @cached_property
    def num_phased_variants_from_reads(self):
        """
        Number of variants explicitly phased with the variant in this
        IsovarResult by shared RNA fragments.

        Returns int
        """
        return len(self.phased_variants_in_supporting_reads)

    @cached_property
    def num_phased_variants_in_protein_sequence(self):
        """
        Number of variants explicitly phased with the variant in this
        IsovarResult by shared RNA fragments which were used to create
        the top protein sequences for this variant and whatever it
        gets phased with.

        Returns int
        """
        return len(self.phased_variants_in_protein_sequence)