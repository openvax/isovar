"""Every protein hypothesis for small variants, with scoped RNA evidence.

`export_protein_hypotheses` turns `run_isovar` results into one JSON-ready
exchange format, ``isovar.protein_hypotheses.v1``. It keeps every protein
Isovar reconstructed for each variant, in Isovar's ranked order, and every
nucleotide translation behind each protein, so synonymous cDNAs stay distinct.
RNA support is given as counts plus an evidence set ID; each set's hashed
read identities (see `isovar.rna_evidence`) are stored once, so consumers can
combine evidence without counting a read twice. The export selects nothing:
ranking and filtering stay later policies, and recorded filter outcomes are
reported, not applied.
"""

import csv
import json
from pathlib import Path

from . import __version__
from .genetic_code import translate_cdna
from .protein_sequence_helpers import covered_protein_groups
from .read_identity import count_reads, fragment_ids, source_read_ids
from .read_metadata import unique_header_entries
from .rna_evidence import content_identifier, segment_support
from .transcript_edit_helpers import categorize_transcript_assembly_edits_from_translation

SCHEMA = "isovar.protein_hypotheses.v1"

_EFFECT_FIELDS = ("gene_name", "gene_id", "transcript_id", "transcript_name", "modifies_protein_sequence")


class _EvidenceSets:
    """Collects each distinct read set once, keyed by its evidence set ID."""

    def __init__(self, scope):
        self.scope, self.sets = scope, {}

    def support(self, reads, description):
        """Counts and evidence set ID; the ID is null when any read lacks an identity."""
        reads = list(reads)
        identities = [source_read_ids(read) for read in reads]
        if not all(identities):
            return dict(scope=description, segments=count_reads(reads), fragments=len(fragment_ids(reads)),
                        evidence_set_id=None, independent_molecules=None)
        evidence = segment_support(self.scope, [key for keys in identities for key in keys])
        self.sets[evidence["evidence_set_id"]] = evidence
        return dict(scope=description, segments=evidence["segments"], fragments=evidence["fragments"],
                    evidence_set_id=evidence["evidence_set_id"], independent_molecules=None)


def _counts(reads):
    reads = list(reads)
    return dict(segments=count_reads(reads), fragments=len(fragment_ids(reads)))


def _variant(variant):
    return dict(reference_name=variant.reference_name, contig=variant.contig, start=variant.start,
                ref=variant.ref, alt=variant.alt, description=variant.short_description)


def _event_id(variant):
    return content_identifier("variant", [variant.reference_name, variant.contig, variant.start,
                                          variant.ref, variant.alt])


def _starts_at_annotated_start_codon(translation):
    """Whether translation begins at the annotated start codon, as Isovar frames it."""
    context, orf = translation.reference_context, translation.variant_orf
    matched_from = len(context.sequence_before_variant_locus) - len(orf.reference_cdna_sequence_before_variant)
    return bool(context.contains_start_codon and matched_from <= context.offset_to_first_complete_codon)


def _check_translation(translation):
    """The exported cDNA, frame and protein must agree (a length cap may shorten the protein)."""
    orf = translation.variant_orf
    amino_acids, stop = translate_cdna(
        orf.cdna_sequence[orf.offset_to_first_complete_codon:], first_codon_is_start=False,
        mitochondrial=translation.reference_context.mitochondrial)
    if (amino_acids[:len(translation.amino_acids)] != translation.amino_acids
            or (translation.ends_with_stop_codon and not (stop and amino_acids == translation.amino_acids))):
        raise ValueError("Translation cDNA, frame and amino acids disagree")


def _edits(translation, known_variants):
    """Each difference from each transcript, with the supplied variant it is, if any."""
    focal = translation.reference_context.variant
    rows = set()
    for transcript in translation.reference_context.transcripts:
        categories = categorize_transcript_assembly_edits_from_translation(
            translation, transcript, known_variants)
        for category, edits in categories.items():
            for edit in edits:
                source = edit.source_variant
                origin = ("unexplained" if source is None else "nominated_variant" if source == focal
                          else "germline_variant" if category == "known_germline" else "co_somatic_variant")
                rows.add((edit.transcript_id, edit.cdna_start, edit.cdna_end, edit.alt_bases, origin,
                          None if source is None else _event_id(source)))
    return [dict(transcript_id=t, cdna_interval=[a, b], alt_bases=alt, origin=origin, source_event_id=source)
            for t, a, b, alt, origin, source in sorted(rows, key=lambda row: row[:5])]


def _translation(translation, event_id, evidence, known_variants):
    _check_translation(translation)
    orf, context = translation.variant_orf, translation.reference_context
    transcripts = sorted(context.transcripts, key=lambda t: t.id)
    offset, cdna = orf.offset_to_first_complete_codon, orf.cdna_sequence
    translated = 3 * (len(translation.amino_acids) + int(translation.ends_with_stop_codon))
    return dict(
        translation_id=content_identifier("translation", [
            event_id, cdna, offset, [t.id for t in transcripts], translation.amino_acids]),
        nucleotide_sequence=cdna, nucleotide_sequence_id=content_identifier("nt", cdna),
        translated_interval=[offset, offset + translated],
        variant_cdna_interval=[orf.variant_cdna_interval_start, orf.variant_cdna_interval_end],
        starts_at_annotated_start_codon=_starts_at_annotated_start_codon(translation),
        reference_context=dict(
            strand=context.strand, transcript_ids=[t.id for t in transcripts],
            transcript_names=[getattr(t, "name", None) for t in transcripts],
            contains_start_codon=context.contains_start_codon,
            overlaps_start_codon=context.overlaps_start_codon,
            contains_five_prime_utr=context.contains_five_prime_utr,
            amino_acids_before_variant=context.amino_acids_before_variant),
        mismatches_before_variant=orf.num_mismatches_before_variant,
        mismatches_after_variant=orf.num_mismatches_after_variant,
        observed_edits=_edits(translation, known_variants),
        rna_support=evidence.support(translation.reads, "reads_assembled_into_this_cdna"))


def _frame_contexts(protein):
    return {(t.reference_context.strand,
             (t.variant_orf.variant_cdna_interval_start - t.variant_orf.offset_to_first_complete_codon) % 3,
             transcript.id)
            for t in protein.translations for transcript in t.reference_context.transcripts}


def _hypothesis_id(protein, event_id):
    return content_identifier("protein_hypothesis", [
        event_id, protein.amino_acids, [protein.mutation_start_idx, protein.mutation_end_idx],
        protein.ends_with_stop_codon, protein.frameshift])


def _protein(protein, rank, event_id, evidence):
    translations = sorted((_translation(t, event_id, evidence, protein.known_variants)
                           for t in protein.translations),
                          key=lambda t: t["translation_id"])
    starts = {t["starts_at_annotated_start_codon"] for t in translations}
    interval = [protein.mutation_start_idx, protein.mutation_end_idx]
    return dict(
        hypothesis_id=_hypothesis_id(protein, event_id), isovar_rank=rank, amino_acids=protein.amino_acids,
        protein_sequence_id=content_identifier("aa", protein.amino_acids),
        mutation_interval=interval, mutant_amino_acids=protein.mutant_amino_acids,
        contains_mutation=protein.contains_mutation, frameshift=protein.frameshift,
        ends_with_stop_codon=protein.ends_with_stop_codon,
        n_terminus=("annotated_start_codon" if starts == {True} else "local_window" if starts == {False}
                    else "mixed"),
        c_terminus="stop_codon" if protein.ends_with_stop_codon else "no_stop_in_context",
        transcript_ids=list(protein.transcript_ids), transcript_names=list(protein.transcript_names),
        gene_ids=list(protein.gene_ids), gene_names=list(protein.gene_names),
        rna_support=evidence.support(protein.supporting_reads, "reads_assembled_into_this_protein"),
        translations=translations)


def _proteins(proteins, event_id, evidence):
    """Every protein, marking the shorter contexts that another one contains."""
    rows = [_protein(p, rank, event_id, evidence) for rank, p in enumerate(proteins, 1)]
    keys = [(_frame_contexts(p), p.frameshift, p.mutation_start_idx,
             p.amino_acids + ("*" if p.ends_with_stop_codon else "")) for p in proteins]
    groups = covered_protein_groups(keys)
    representatives = {i for i, _ in groups}
    for j, row in enumerate(rows):
        row["representative"] = j in representatives
        row["covered_by"] = [
            dict(hypothesis_id=rows[i]["hypothesis_id"], protein_interval=[
                proteins[i].mutation_start_idx - proteins[j].mutation_start_idx,
                proteins[i].mutation_start_idx - proteins[j].mutation_start_idx + len(proteins[j])])
            for i, covered in groups if i != j and j in covered]
    return rows


def _edit_attribution(proteins):
    """How many supplied variants the edits were checked against; null if none were."""
    known = {getattr(p, "known_variants", None) for p in proteins} - {None}
    if len(known) != 1:
        return None
    known, = known
    return dict(somatic_variants=len(known.somatic), germline_variants=len(known.germline))


def _completeness(result):
    """Whether the returned proteins are all that Isovar found, given its cap."""
    settings = result.protein_sequence_settings
    if settings is None:
        return None, None
    limit = settings["max_protein_sequences_per_variant"]
    return limit, not limit or len(result.sorted_protein_sequences) < limit


def _event(result, evidence):
    variant, reads = result.variant, result.read_evidence
    event_id = _event_id(variant)
    effect = result.predicted_effect
    limit, complete = _completeness(result)
    return dict(
        event_id=event_id, variant=_variant(variant),
        reference_prediction=None if effect is None else dict(
            effect_class=type(effect).__name__, description=getattr(effect, "short_description", None),
            **{name: getattr(effect, name, None) for name in _EFFECT_FIELDS}),
        allele_support=dict(
            alt=evidence.support(reads.alt_reads, "alt_allele_reads_at_locus"),
            ref=_counts(reads.ref_reads), other=_counts(reads.other_reads),
            total=_counts(list(reads.ref_reads) + list(reads.alt_reads) + list(reads.other_reads))),
        filters=dict(values=dict(result.filter_values), passes_all_filters=result.passes_all_filters),
        phased_variants=dict(
            supporting_reads=sorted(_event_id(v) for v in result.phased_variants_in_supporting_reads),
            top_protein_sequence=sorted(_event_id(v) for v in result.phased_variants_in_protein_sequence)),
        protein_sequence_limit=limit, protein_hypotheses_complete=complete,
        edit_attribution=_edit_attribution(result.sorted_protein_sequences),
        protein_sequence_settings=result.protein_sequence_settings,
        protein_hypotheses=_proteins(result.sorted_protein_sequences, event_id, evidence))


def read_group_metadata(alignment_header):
    """Sample, library and platform of each read group in a SAM header.

    Parameters
    ----------
    alignment_header : pysam.AlignmentHeader or dict
        The header of the alignment file, or its ``to_dict()``.

    Returns
    -------
    dict
        Read group ID to ``sample``, ``library`` and ``platform`` (the ``SM``,
        ``LB`` and ``PL`` fields; null when absent). IDs that appear more than
        once are omitted, since they cannot identify one read group.
    """
    header = alignment_header if isinstance(alignment_header, dict) else alignment_header.to_dict()
    return {group: dict(sample=entry.get("SM"), library=entry.get("LB"), platform=entry.get("PL"))
            for group, entry in unique_header_entries(header.get("RG", [])).items()}


def _read_groups(results, alignment_header):
    groups = set()
    for result in results:
        evidence = result.read_evidence
        for read in list(evidence.ref_reads) + list(evidence.alt_reads) + list(evidence.other_reads):
            groups.update(key[0] for key in source_read_ids(read))
    metadata = {} if alignment_header is None else read_group_metadata(alignment_header)
    return {group: metadata.get(group) for group in sorted(groups)}


def export_protein_hypotheses(isovar_results, *, sample_id, source, alignment_header=None):
    """Export every protein hypothesis and translation, with scoped RNA evidence.

    Parameters
    ----------
    isovar_results : iterable of IsovarResult
        Results from `run_isovar`. Proteins beyond the creator's
        ``max_protein_sequences_per_variant`` were never kept, so run with 0
        to export them all; each event says whether its list is complete.
    sample_id : str
        The sample the reads came from.
    source : str
        Identity of the read set, such as the BAM path or an accession. With
        `sample_id` it scopes the read identities: exports with the same
        scope give the same read the same ID.
    alignment_header : pysam.AlignmentHeader or dict, optional
        Header of the alignment file, used to report each read group's sample
        and library. Without it their metadata are null (unknown).

    Returns
    -------
    dict
        ``isovar.protein_hypotheses.v1``: one event per result, in input
        order, each with its protein hypotheses in Isovar's ranked order and
        their translations, plus ``evidence_sets``, the hashed read IDs of
        every support by evidence set ID. Intervals are 0-based and
        half-open. Read and fragment counts are not molecule counts.

    Raises
    ------
    ValueError
        If a translation's cDNA and frame do not give its amino acids.
    """
    results = list(isovar_results)
    scope = [sample_id, source]
    evidence = _EvidenceSets(scope)
    events = [_event(result, evidence) for result in results]
    return dict(
        schema=SCHEMA, isovar_version=__version__, sample_id=sample_id, source=source,
        evidence_scope=scope, interval_convention="zero_based_half_open",
        evidence_identity_policy="sha256_of_domain_and_canonical_json; sample/source scoped",
        read_group_metadata="alignment_header" if alignment_header is not None else "not_supplied",
        read_groups=_read_groups(results, alignment_header),
        interpretation="RNA-derived protein hypotheses in Isovar's ranked order; nothing is selected. "
                       "Counts are sequenced segments and fragments, not molecules or abundance.",
        events=events,
        evidence_sets={key: evidence.sets[key] for key in sorted(evidence.sets)})


_TSV_COLUMNS = (
    "schema", "sample_id", "source", "event_id", "variant", "hypothesis_id", "isovar_rank", "representative",
    "amino_acids", "mutation_start", "mutation_end", "ends_with_stop_codon", "frameshift", "n_terminus",
    "gene_names", "protein_transcript_ids", "protein_segments", "protein_fragments", "protein_evidence_set_id",
    "translation_id", "nucleotide_sequence", "translated_start", "translated_end",
    "starts_at_annotated_start_codon", "context_transcript_ids", "translation_segments",
    "translation_fragments", "translation_evidence_set_id", "protein_hypotheses_complete")


def _tsv_rows(export):
    for event in export["events"]:
        for protein in event["protein_hypotheses"]:
            for translation in protein["translations"]:
                yield dict(
                    schema=export["schema"], sample_id=export["sample_id"], source=export["source"],
                    event_id=event["event_id"], variant=event["variant"]["description"],
                    hypothesis_id=protein["hypothesis_id"], isovar_rank=protein["isovar_rank"],
                    representative=protein["representative"],
                    amino_acids=protein["amino_acids"], mutation_start=protein["mutation_interval"][0],
                    mutation_end=protein["mutation_interval"][1],
                    ends_with_stop_codon=protein["ends_with_stop_codon"], frameshift=protein["frameshift"],
                    n_terminus=protein["n_terminus"], gene_names=";".join(protein["gene_names"]),
                    protein_transcript_ids=";".join(protein["transcript_ids"]),
                    protein_segments=protein["rna_support"]["segments"],
                    protein_fragments=protein["rna_support"]["fragments"],
                    protein_evidence_set_id=protein["rna_support"]["evidence_set_id"] or "",
                    translation_id=translation["translation_id"],
                    nucleotide_sequence=translation["nucleotide_sequence"],
                    translated_start=translation["translated_interval"][0],
                    translated_end=translation["translated_interval"][1],
                    starts_at_annotated_start_codon=translation["starts_at_annotated_start_codon"],
                    context_transcript_ids=";".join(translation["reference_context"]["transcript_ids"]),
                    translation_segments=translation["rna_support"]["segments"],
                    translation_fragments=translation["rna_support"]["fragments"],
                    translation_evidence_set_id=translation["rna_support"]["evidence_set_id"] or "",
                    protein_hypotheses_complete=("" if event["protein_hypotheses_complete"] is None
                                                 else event["protein_hypotheses_complete"]))


def write_protein_hypotheses(export, path):
    """Write an export as JSON, plus a TSV with one row per translation.

    Parameters
    ----------
    export : dict
        An ``isovar.protein_hypotheses.v1`` export.
    path : str or Path
        The JSON path; the TSV is written beside it with ``.tsv`` in place of
        the ``.json`` suffix (or appended). Existing files are replaced.

    Returns
    -------
    dict
        The written ``json`` and ``tsv`` paths.
    """
    if export.get("schema") != SCHEMA:
        raise ValueError("Expected an %s export" % SCHEMA)
    path = Path(path).expanduser()
    tsv = path.with_suffix(".tsv") if path.suffix == ".json" else Path(str(path) + ".tsv")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(export, indent=2) + "\n")
    with tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_TSV_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(_tsv_rows(export))
    return dict(json=path, tsv=tsv)
