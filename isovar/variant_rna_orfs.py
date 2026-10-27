"""Mutation-supporting RNA interpreted in supplied, validated reference ORFs.

An ORFdata TranscriptORF can describe an ordinary CDS or a noncanonical ORF
on a noncoding transcript. Its reference observations supply a frame hypothesis,
not evidence that this sample initiates translation or produces a protein.
"""

from orfdata import ReferenceContext as ORFReferenceContext

from .dna import reverse_complement_dna
from .genetic_code import translate_cdna
from .protein_sequence_creator import ProteinSequenceCreator
from .reference_context import ReferenceContext
from .reference_sequence_key import ReferenceSequenceKey
from .rna_evidence import EvidenceSets, content_identifier
from .read_identity import fragment_ids
from .transcript_compatibility import annotate_reads_with_transcript_compatibility
from .translation_helpers import find_mutant_amino_acid_interval
from .variant_helpers import (
    NoExonicReferenceBases, interbase_range_affected_by_variant_on_transcript,
    require_literal_variant, trim_variant,
)
from .variant_orf_helpers import match_variant_sequence_to_reference_context


SCHEMA = 'isovar.variant_rna_orfs.v1'


def _context(variant, reference, size):
    """Use the selected ORF frame without changing its parent's annotation."""
    normalize = lambda c: str(c).removeprefix('chr').replace('MT', 'M')
    if normalize(variant.contig) != normalize(reference.contig):
        return None, None, 'outside_reference_locus'
    if not reference.start_complete or not reference.stop_complete:
        return None, None, 'incomplete_reference_orf'
    try:
        start, end = interbase_range_affected_by_variant_on_transcript(
            variant, reference.transcript)
    except NoExonicReferenceBases:
        return None, None, 'no_exonic_reference_bases'
    if start < reference.start + 3:
        return None, None, ('reference_start_affected' if end > reference.start
                            else 'outside_reference_orf')
    if start >= reference.end:
        return None, None, 'outside_reference_orf'
    key = ReferenceSequenceKey.from_variant_and_transcript(variant, reference.transcript, size)
    if key is None:
        return None, None, 'reference_sequence_or_allele_unresolved'
    prefix_start = start - len(key.sequence_before_variant_locus)
    contains_start = prefix_start <= reference.start
    offset = reference.start - prefix_start if contains_start else (reference.start - prefix_start) % 3
    first_amino = max(0, (prefix_start + offset - reference.start) // 3)
    before_amino = (start - reference.start) // 3
    context = ReferenceContext(
        strand=reference.strand,
        sequence_before_variant_locus=key.sequence_before_variant_locus,
        sequence_at_variant_locus=key.sequence_at_variant_locus,
        sequence_after_variant_locus=key.sequence_after_variant_locus,
        offset_to_first_complete_codon=offset,
        contains_start_codon=contains_start,
        overlaps_start_codon=prefix_start < reference.start + 3,
        contains_five_prime_utr=prefix_start < reference.start,
        amino_acids_before_variant=reference.protein_sequence[first_amino:before_amino],
        variant=variant, transcripts=[reference.transcript],
        translation_table=reference.translation_table)
    return context, (start, end), None


def _translate(sequence, table, starts_at_reference_start):
    amino, stop = translate_cdna(sequence, translation_table=table)
    # A source-supported non-AUG start is conditionally decoded as methionine.
    # Its sequence was checked against the reference by frame matching.
    if starts_at_reference_start and amino:
        amino = 'M' + amino[1:]
    return amino, stop


def _interval_support(sequence, start, end, reverse, evidence):
    """Exact single-read spans and fragment tiling are separate measurements."""
    if reverse:
        start, end = len(sequence) - end, len(sequence) - start
    covered = [set() for _ in range(end - start)]
    spanning = []
    for read in sequence.reads:
        shift = len(sequence.prefix) - len(read.prefix)
        left, right = max(start, shift), min(end, shift + len(read.sequence))
        if right <= left:
            continue
        if sequence.sequence[left:right] != read.sequence[left - shift:right - shift]:
            continue
        fragments = fragment_ids([read])
        for position in range(left, right):
            covered[position - start].update(fragments)
        if left == start and right == end:
            spanning.append(read)
    depths = [len(names) for names in covered]
    if reverse:
        depths.reverse()
    return dict(scope='observed_coding_interval_including_stop_when_present',
                fragment_depth=depths, minimum_fragment_depth=min(depths) if depths else 0,
                fully_tiled=bool(depths) and all(depths),
                single_read_span_support=evidence.support(spanning))


def _candidate(sequence, context, locus, reference, model_id, rna_id, evidence, creator):
    framed = match_variant_sequence_to_reference_context(
        sequence, context, creator.min_transcript_prefix_length,
        creator.max_transcript_mismatches, creator.count_mismatches_after_variant)
    if framed is None:
        return None
    cdna, offset = framed.cdna_sequence, framed.offset_to_first_complete_codon
    reference_prefix_start = locus[0] - len(framed.reference_cdna_sequence_before_variant)
    starts = reference_prefix_start <= reference.start
    if starts and cdna[offset:offset + 3] != reference.start_codon:
        return None
    protein, stop = _translate(cdna[offset:], reference.translation_table, starts)
    if creator.protein_sequence_length and len(protein) > creator.protein_sequence_length:
        protein, stop = protein[:creator.protein_sequence_length], False
    begin, end = framed.variant_cdna_interval_start, framed.variant_cdna_interval_end
    mutation_start, mutation_end, frameshift = find_mutant_amino_acid_interval(
        cdna, offset, begin, end, len(context.sequence_at_variant_locus), len(protein))
    if stop and mutation_start == len(protein):
        # A focal premature stop has no amino-acid residue to highlight.
        mutation_end = mutation_start
    contains = (0 < mutation_start <= len(protein) if mutation_start == mutation_end and stop else
                0 < mutation_start < len(protein) if mutation_start == mutation_end
                else mutation_start < len(protein))
    if not contains:
        return None
    counterfactual = cdna[:begin] + context.sequence_at_variant_locus + cdna[end:]
    counter_protein, _ = _translate(counterfactual[offset:], reference.translation_table, starts)
    if creator.protein_sequence_length:
        counter_protein = counter_protein[:creator.protein_sequence_length]
    reference_codon_start = max(reference.start, reference_prefix_start +
                                (reference.start - reference_prefix_start) % 3)
    reference_end = min(reference.end, locus[1] + len(cdna) - end)
    reference_window, _ = _translate(reference.transcript_sequence[reference_codon_start:reference_end],
                                     reference.translation_table, reference_codon_start == reference.start)
    if creator.protein_sequence_length:
        reference_window = reference_window[:creator.protein_sequence_length]
    oriented = reverse_complement_dna(sequence.sequence) if reference.strand == '-' else sequence.sequence
    oriented_prefix = len(sequence.suffix) if reference.strand == '-' else len(sequence.prefix)
    shift = oriented_prefix - begin
    if oriented[shift:shift + len(cdna)] != cdna:
        raise ValueError('Framed RNA does not map back to the original assembly')
    coding_end = offset + 3 * (len(protein) + int(stop))
    support = _interval_support(sequence, shift + offset, shift + coding_end,
                                reference.strand == '-', evidence)
    return dict(
        candidate_id=content_identifier('orf_candidate', [model_id, rna_id, cdna, offset, protein]),
        model_id=model_id, rna_candidate_id=rna_id, orf_id=reference.orf_id,
        transcript_id=reference.transcript_id, gene_name=reference.gene_name,
        nucleotide_sequence=cdna, coding_interval=[offset, coding_end],
        variant_cdna_interval=[begin, end], protein_sequence=protein,
        mutation_interval=[mutation_start, mutation_end], frameshift=frameshift,
        starts_at_reference_start=starts, ends_with_stop_codon=stop,
        sequence_extent='complete_coding_candidate' if starts and stop else 'partial_coding_candidate',
        reference_protein_window=reference_window,
        reference_protein_start_0based=(reference_codon_start - reference.start) // 3,
        reference_allele_counterfactual_protein=counter_protein,
        focal_amino_acid_change=protein != counter_protein,
        comparison_scope='focal allele replaced on the same observed RNA background',
        mismatches_before_variant=framed.num_mismatches_before_variant,
        mismatches_after_variant=framed.num_mismatches_after_variant,
        rna_support=evidence.support(sequence.reads), coding_support=support,
        initiation_observed=False, translation_observed=False,
        evidence_state='reference_orf_frame_on_mutation_supporting_rna')


def reconstruct_variant_orfs(variant, read_evidence, reference_orfs, *, sample_id, source,
                             protein_sequence_creator=None):
    """Reconstruct local RNA candidates in supplied ORFdata coding frames.

    Parameters
    ----------
    variant : varcode.Variant
        Literal small variant. Structural events use ``reconstruct_sv_rna``.
    read_evidence : ReadEvidence
        Original observations grouped by their focal allele.
    reference_orfs : iterable of orfdata.TranscriptORF
        Sequence-validated models, including ordinary CDS competitors. No
        catalogue download, canonical-CDS filter or strongest-frame choice.
    sample_id, source : str
        Independent RNA sample and library identities for evidence accounting.
    protein_sequence_creator : ProteinSequenceCreator, optional
        Assembly, matching and local protein-window budgets. Its protein-count
        cap does not discard models or RNA alternatives in this assessment.

    Returns
    -------
    dict
        JSON-ready RNA/model/protein candidates and scoped original witnesses.
        Missing candidates and unresolved models remain explicit. Reference
        frame assignment does not establish initiation, translation or HLA display.
    """
    require_literal_variant(variant)
    from . import __version__
    if (read_evidence.trimmed_base1_start, read_evidence.trimmed_ref,
            read_evidence.trimmed_alt) != trim_variant(variant):
        raise ValueError('Read evidence does not describe the supplied variant')
    if not sample_id or not source:
        raise ValueError('sample_id and source are required')
    creator = protein_sequence_creator or ProteinSequenceCreator(max_protein_sequences_per_variant=0)
    references = list(reference_orfs)
    variant_context = ORFReferenceContext.from_genome(variant.genome)
    transcripts = {}
    for reference in references:
        if reference.context.key != variant_context.key:
            raise ValueError('Variant and ORF reference contexts differ')
        transcript = reference.transcript
        prior = transcripts.get(transcript.id)
        if prior is not None and (prior.sequence != transcript.sequence or
                                  prior.exon_intervals != transcript.exon_intervals):
            raise ValueError('Ambiguous supplied transcript identity')
        transcripts[transcript.id] = transcript
    evidence = EvidenceSets([sample_id, source])
    reads = list(read_evidence.alt_reads)
    if reads:
        left, right = creator._variant_sequence_creator.flanking_sequence_lengths(len(reads[0].allele))
        classified = annotate_reads_with_transcript_compatibility(reads, transcripts.values(), left, right)
    else:
        classified = []
    compatible = [read for read in classified if read.compatible_transcript_ids]
    excluded = [read for read in classified if not read.compatible_transcript_ids]
    sequences = creator.variant_sequences_from_reads(variant, compatible)
    unresolved_sequences = creator.variant_sequences_from_reads(
        variant, [read.with_compatible_transcript_ids(None) for read in excluded])
    rna_rows, sequence_ids = [], {}
    for state, group in (('reference_path_compatible', sequences),
                         ('no_supplied_reference_path', unresolved_sequences)):
        for sequence in group:
            support = evidence.support(sequence.reads)
            placements = sorted(content_identifier('placement', [read.source_alignments,
                                                                  read.source_alignment_paths])
                                for read in sequence.reads)
            ident = content_identifier('rna_candidate', [state, sequence.sequence, len(sequence.prefix),
                                                         sorted(sequence.compatible_transcript_ids or ()),
                                                         support['evidence_set_id'], placements])
            sequence_ids[id(sequence)] = ident
            rna_rows.append(dict(rna_candidate_id=ident, state=state, sequence=sequence.sequence,
                                 variant_interval=[len(sequence.prefix), len(sequence.prefix) + len(sequence.alt)],
                                 compatible_transcript_ids=sorted(sequence.compatible_transcript_ids or ()),
                                 rna_support=support))
    models, candidates = [], []
    size = max(creator._cdna_sequence_length, creator.min_transcript_prefix_length)
    seen_models = set()
    for reference in references:
        model = reference.to_dict()
        model_id = content_identifier('orf_model', model)
        if model_id in seen_models:
            raise ValueError('Duplicate supplied ORF model')
        seen_models.add(model_id)
        context, locus, reason = _context(variant, reference, size)
        matching_sequences = [s for s in sequences if reference.transcript_id in
                              (s.compatible_transcript_ids or ())]
        accepted = []
        if context is not None:
            for sequence in matching_sequences:
                candidate = _candidate(sequence, context, locus, reference, model_id,
                                       sequence_ids[id(sequence)], evidence, creator)
                if candidate is not None:
                    accepted.append(candidate)
        if not accepted and reason is None:
            reason = ('no_alt_reads' if not reads else
                      'no_compatible_rna_reads' if not any(
                          reference.transcript_id in (r.compatible_transcript_ids or ()) for r in classified)
                      else 'no_rna_sequence_passed_coverage_and_assembly' if not matching_sequences
                      else 'no_retained_rna_frame_match')
        models.append(dict(model_id=model_id, reference=model,
                           state='rna_frame_candidates' if accepted else 'unresolved',
                           reason=None if accepted else reason,
                           compatible_rna_candidate_count=len(matching_sequences),
                           candidate_ids=[c['candidate_id'] for c in accepted]))
        candidates.extend(accepted)
    candidates = list({c['candidate_id']: c for c in candidates}.values())
    candidates.sort(key=lambda c: (-c['rna_support']['fragments'], -len(c['protein_sequence']), c['candidate_id']))
    return dict(schema=SCHEMA, producer=dict(name='isovar', version=__version__),
                sample_id=sample_id, source=source,
                variant=dict(reference_name=variant.reference_name, contig=variant.contig,
                             start=variant.start, ref=variant.ref, alt=variant.alt),
                assembly_state='assembled' if rna_rows else 'no_retained_rna' if reads else 'no_alt_reads',
                settings=creator.settings(), allele_support=evidence.support(reads),
                protein_candidate_limit_applied=None,
                reference_incompatible_support=evidence.support(excluded),
                rna_candidates=rna_rows, models=models, protein_candidates=candidates,
                evidence_sets=evidence.sets,
                limitations=['Local allele-centred assembly; no whole-transcript discovery claim.',
                             'Reference ORF evidence is not mutant initiation or translation evidence.'])
