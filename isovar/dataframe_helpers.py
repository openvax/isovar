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

import pandas as pd
from varcode.errors import InvalidContigError

from .allele_read import AlleleRead
from .read_identity import allele_read_mate_label, count_reads, fragment_ids
from .common import list_to_string
from .dataframe_builder import DataFrameBuilder
from .isovar_result import IsovarResult
from .locus_read import LocusRead
from .protein_sequence import ProteinSequence
from .read_collector import ReadCollector
from .read_evidence import ReadEvidence
from .reference_context import ReferenceContext
from .translation import Translation
from .variant_sequence import VariantSequence


def variant_gene_names(variant):
    """Return semicolon-joined gene names overlapping ``variant``, or ``""``.

    Variants on synthetic or unannotated contigs return ``""`` rather than
    raising ``InvalidContigError``.
    """
    try:
        return ";".join(variant.gene_names)
    except InvalidContigError:
        return ""


def dataframe_from_generator(
        element_class,
        variant_and_elements_generator,
        **kwargs):
    """
    Creates a DataFrame from a generator whose elements
    are varcode.Variant objects paired with objects
    of the `element_class` type.

    Parameters
    ----------
    element_class : type

    variant_and_elements_generator : generator
        Elements are (varcode.Variant, element_class)

    **kwargs : dict
        Additional arguments to pass to DataFrameBuilder.

    Returns
    -------
    pandas.DataFrame
    """
    builder = DataFrameBuilder(element_class, **kwargs)
    for variant, x in variant_and_elements_generator:
        if isinstance(x, (list, tuple)):
            builder.add_many(variant, x)
        else:
            builder.add(variant, x)
    return builder.to_dataframe()


def protein_sequences_generator_to_dataframe(variant_and_protein_sequences_generator):
    """
    Given a generator which yields (Variant, [ProteinSequence]) elements,
    returns a pandas.DataFrame
    """
    return dataframe_from_generator(
        element_class=ProteinSequence,
        variant_and_elements_generator=variant_and_protein_sequences_generator,
        exclude={"known_variants"},
        rename_dict={"translations": "num_translations"},
        extra_column_fns=dict(
            gene_names=lambda _, x: ";".join(x.gene_names),
            transcript_ids=lambda _, x: ";".join(x.transcript_ids),
            transcript_names=lambda _, x: ";".join(x.transcript_names),
            num_supporting_reads=lambda _, x: x.num_supporting_reads,
            num_supporting_fragments=lambda _, x: x.num_supporting_fragments))


def allele_counts_dataframe(read_evidence_generator):
    """
    Creates a DataFrame containing number of reads supporting the
    ref vs. alt alleles for each variant.
    """
    return dataframe_from_generator(
        element_class=ReadEvidence,
        variant_and_elements_generator=read_evidence_generator,
        converters={
            "ref_reads": count_reads,
            "alt_reads": count_reads,
            "other_reads": count_reads,
        },
        rename_dict={
            "ref_reads": "num_ref_reads",
            "alt_reads": "num_alt_reads",
            "other_reads": "num_other_reads",
        },
        extra_column_fns={
            "gene_names": lambda variant, _: variant_gene_names(variant),
            "num_ref_fragments": lambda _, x: len(fragment_ids(x.ref_reads)),
            "num_alt_fragments": lambda _, x: len(fragment_ids(x.alt_reads)),
            "num_other_fragments": lambda _, x: len(fragment_ids(x.other_reads))
        })


def allele_reads_to_dataframe(variants_and_allele_reads, allele_class=None):
    """
    Parameters
    ----------
    variants_and_allele_reads : sequence
        List or generator of pairs whose first element is a Variant and
        whose second element is a sequence of AlleleRead objects.
    allele_class : str, optional
        If provided (e.g. ``"alt"``), include an ``allele_class`` column set to
        this constant for every row.
    """
    extra_column_fns = {
        "gene_names": lambda v, _: variant_gene_names(v),
        "mate": lambda _, r: allele_read_mate_label(r),
        "is_secondary": lambda _, r: getattr(r, "is_secondary", None),
    }
    if allele_class is not None:
        extra_column_fns = {
            "allele_class": lambda _, _r: allele_class,
            **extra_column_fns,
        }
    df_builder = DataFrameBuilder(
        AlleleRead,
        exclude={"source_read_count", "source_alignments", "source_alignment_paths"},
        extra_column_fns=extra_column_fns)
    for variant, allele_reads in variants_and_allele_reads:
        df_builder.add_many(variant, allele_reads)
    return df_builder.to_dataframe()


def read_evidence_to_allele_reads_dataframe(
        variant_and_read_evidence_generator,
        allele_classes=("ref", "alt", "other")):
    """Build a per-read DataFrame from ``(Variant, ReadEvidence)`` pairs.

    Preserves each observation's ``allele_class`` (``"ref"``, ``"alt"``, or
    ``"other"``), mate label (``"R1"``, ``"R2"``, ``"R1+R2"``, ``"unpaired"``,
    or ``None``), and ``is_secondary`` alignment flag so rows reconcile
    directly with ``allele-counts``.

    Parameters
    ----------
    variant_and_read_evidence_generator : iterable
        Sequence or generator of ``(Variant, ReadEvidence)`` pairs.
    allele_classes : sequence of str, optional
        Subset and order of allele buckets to export (default:
        ``("ref", "alt", "other")``).

    Returns
    -------
    pandas.DataFrame
    """
    valid = {"ref", "alt", "other"}
    invalid = [name for name in allele_classes if name not in valid]
    if invalid:
        raise ValueError("Unknown allele_classes %s; expected subset of %s" % (
            invalid, sorted(valid)))
    rows = []
    for variant, read_evidence in variant_and_read_evidence_generator:
        gene_names = variant_gene_names(variant)
        for allele_class in allele_classes:
            for read in getattr(read_evidence, allele_class + "_reads"):
                rows.append({
                    "variant": variant.short_description,
                    "chr": variant.original_contig,
                    "pos": variant.original_start,
                    "ref": variant.original_ref,
                    "alt": variant.original_alt,
                    "gene_names": gene_names,
                    "allele_class": allele_class,
                    "name": read.name,
                    "mate": allele_read_mate_label(read),
                    "is_secondary": getattr(read, "is_secondary", None),
                    "prefix": read.prefix,
                    "allele": read.allele,
                    "suffix": read.suffix,
                    "sequence": read.sequence,
                })
    columns = [
        "variant", "chr", "pos", "ref", "alt",
        "gene_names", "allele_class", "name", "mate", "is_secondary",
        "prefix", "allele", "suffix", "sequence",
    ]
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame.from_records(rows, columns=columns)


def locus_reads_dataframe(alignments, chromosome, base0_start, base0_end, *args, **kwargs):
    """
    Traverse an alignment file (typically a BAM) to find all the reads
    overlapping a specified locus.

    Extra parameters are passed to ReadCollector.
    """
    df_builder = DataFrameBuilder(
        LocusRead,
        exclude={"source_read_count", "source_alignments", "source_alignment_paths", "is_primary", "source_read_views"},
        variant_columns=False,
        converters={
            "reference_positions": list_to_string,
            "quality_scores": list_to_string,
        })
    read_creator = ReadCollector(*args, **kwargs)
    for locus_read in read_creator.get_locus_reads(
            alignments, chromosome, base0_start, base0_end):
        df_builder.add(variant=None, element=locus_read)
    return df_builder.to_dataframe()


def variants_to_reference_contexts_dataframe(
        variant_and_reference_contexts_generator):
    """
    Given a generator of (Variant, [ReferenceContext]) pairs, create a
    DataFrame.

    Returns a DataFrame with {"variant", "chr", "pos", "ref", "alt"} columns
    for variants, as well as the fields of ReferenceContext and normalized
    ``gene_names``, ``transcript_ids``, and ``transcript_names`` columns.
    """
    df_builder = DataFrameBuilder(
        ReferenceContext,
        exclude=["variant", "transcripts"],
        extra_column_fns={
            "gene_names": lambda variant, _: variant_gene_names(variant),
            "transcript_ids": lambda _, rc: ";".join(t.id for t in rc.transcripts),
            "transcript_names": lambda _, rc: ";".join(t.name for t in rc.transcripts),
        })
    for variant, reference_contexts in variant_and_reference_contexts_generator:
        df_builder.add_many(variant, reference_contexts)
    return df_builder.to_dataframe()


def variant_sequences_generator_to_dataframe(variant_sequences_generator):
    """
    Creates a dataframe from a generator which yields
    (Variant, [VariantSequence]) pairs.

    Returns pandas.DataFrame
    """
    return dataframe_from_generator(
        VariantSequence,
        variant_sequences_generator,
        rename_dict={
            "alt": "allele",
            "reads": "num_supporting_reads",
            "compatible_transcript_ids": "transcript_ids",
        },
        converters={
            "reads": count_reads,
            "compatible_transcript_ids": lambda ids: None if ids is None else ";".join(sorted(ids)),
        },
        extra_column_fns={
            "num_supporting_fragments": lambda _, x: len(fragment_ids(x.reads)),
            "gene_names": lambda variant, _: variant_gene_names(variant),
        })


def translations_generator_to_dataframe(translations_generator):
    """
    Given a generator of (Variant, [Translation]) pairs,
    returns a DataFrame of translated protein fragments with columns
    for each field of a Translation object (and variant/chr/pos/ref/alt).
    """
    def translation_gene_names(variant, translation):
        names = sorted({
            transcript.gene.name
            for transcript in translation.reference_context.transcripts
            if getattr(transcript, "gene", None) is not None and getattr(transcript.gene, "name", None)
        })
        return ";".join(names) if names else variant_gene_names(variant)

    return dataframe_from_generator(
        element_class=Translation,
        variant_and_elements_generator=translations_generator,
        exclude=["reference_context"],
        converters={
            "untrimmed_variant_sequence": lambda vs: vs.sequence,
            "variant_orf": (
                lambda variant_orf: variant_orf.in_frame_cdna_sequence),
        },
        rename_dict={
            "variant_orf": "in_frame_cdna_sequence",
        },
        extra_column_fns={
            "gene_names": translation_gene_names,
            "transcript_ids": lambda _, t: ";".join(
                transcript.id for transcript in t.reference_context.transcripts),
            "transcript_names": lambda _, t: ";".join(
                transcript.name for transcript in t.reference_context.transcripts),
            "num_supporting_reads": lambda _, t: count_reads(t.reads),
            "num_supporting_fragments": lambda _, t: len(fragment_ids(t.reads)),
        })


def isovar_results_to_dataframe(isovar_results):
    """
    Create a DataFrame from a sequence of IsovarResult objects.

    Parameters
    ----------
    isovar_results : list or generator of IsovarResult

    Returns pandas.DataFrame
    """
    records = [isovar_result.to_record() for isovar_result in isovar_results]
    if not records:
        # Keep the header so that an empty result is still a valid table.
        return pd.DataFrame(columns=IsovarResult.record_columns())
    return pd.DataFrame.from_records(records)
