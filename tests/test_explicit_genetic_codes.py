"""Source codes survive RNA reconstruction and Varcode's ORF interchange (#57)."""

from dataclasses import dataclass, replace
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from varcode import Variant

from isovar.genetic_code import translate_cdna
from isovar.protein_hypotheses import _check_translation, _translation, _EvidenceSets
from isovar.protein_sequence_creator import ProteinSequenceCreator
from isovar.reference_context import ReferenceContext
from isovar.reference_coding_sequence_key import ReferenceCodingSequenceKey
from isovar.variant_sequence import VariantSequence
from .test_protein_hypotheses import Gene, Transcript


@dataclass(frozen=True, order=True)
class CodingTranscript(Transcript):
    translation_table: int = 1


def context(table, transcripts=None):
    return ReferenceContext(
        strand="+", sequence_before_variant_locus="ATG", sequence_at_variant_locus="AAA",
        sequence_after_variant_locus="CTGTGAAGGTAA", offset_to_first_complete_codon=0,
        contains_start_codon=True, overlaps_start_codon=True, contains_five_prime_utr=False,
        amino_acids_before_variant="M", variant=Variant("MT", 103, "AAA", "GAA", "GRCh38"),
        transcripts=transcripts if transcripts is not None else [CodingTranscript(
            "tx", "TX", Gene("gene", "GENE"), translation_table=table)])


@pytest.mark.parametrize("table,protein", [(1, "MEL"), (2, "MELW"), (3, "METWR")])
def test_rna_translation_and_export_use_the_explicit_code(table, protein):
    reference = context(table)
    sequence = VariantSequence("ATG", "GAA", "CTGTGAAGGTAA", [])
    creator = ProteinSequenceCreator(min_transcript_prefix_length=0, protein_sequence_length=None)
    translated = creator.translation_from_variant_sequence_and_reference_context(sequence, reference)
    assert translated.amino_acids == protein
    assert translated.ends_with_stop_codon
    assert translated.reference_context.translation_table == table
    _check_translation(translated)
    exported = _translation(translated, "test-event", _EvidenceSets("synthetic-test"), None, None)
    assert exported["reference_context"]["translation_table"] == table
    assert json.loads(json.dumps(exported))["reference_context"]["translation_table"] == table


def test_equal_sequence_contexts_with_different_codes_do_not_group():
    fields = dict(strand="+", sequence_before_variant_locus="ATG", sequence_at_variant_locus="AAA",
                  sequence_after_variant_locus="TAA", offset_to_first_complete_codon=0,
                  contains_start_codon=True, overlaps_start_codon=True, contains_five_prime_utr=False,
                  amino_acids_before_variant="M")
    keys = [ReferenceCodingSequenceKey(**fields, translation_table=table) for table in (1, 2, 3)]
    assert len(set(keys)) == 3
    for key in keys:
        result = ReferenceContext.from_reference_coding_sequence_key(key, context(1).variant, ())
        assert result.translation_table == key.translation_table
    with pytest.raises(ValueError, match="conflicting genetic codes"):
        context(None, context(1).transcripts + context(3).transcripts)


@pytest.mark.parametrize("table", [0, 999, True, "2", 27])
def test_unsupported_or_context_dependent_codes_never_fall_back(table):
    with pytest.raises(ValueError, match="translation table|recoding"):
        translate_cdna("ATGTGA", translation_table=table)


def test_legacy_boolean_is_still_explicitly_table_two():
    assert translate_cdna("ATATGAAGG", mitochondrial=True) == ("MW", True)
    with pytest.raises(ValueError, match="conflicting"):
        translate_cdna("ATATGAAGG", mitochondrial=True, translation_table=3)
    assert translate_cdna("GTGCTGTGAAGATAG", translation_table=3, first_codon_is_start=True) == ("MTWR", True)


def test_reference_scope_not_a_mitochondrial_name_alone_selects_the_code():
    transcript = SimpleNamespace(contig="MT", genome=SimpleNamespace(species="saccharomyces_cerevisiae"))
    assert context(None, (transcript,)).translation_table == 3
    transcript.genome = SimpleNamespace()
    with pytest.raises(ValueError, match="No mitochondrial genetic-code assignment"):
        context(None, (transcript,))


EXAMPLES = json.loads((Path(__file__).parent / "data/references/mitochondrial.json").read_text())["examples"]


@pytest.fixture(params=EXAMPLES, ids=[e["name"] for e in EXAMPLES])
def native_orf(request, tmp_path):
    """Pinned native transcript slices, copied unchanged from ORFdata's qualification fixtures."""
    from orfdata import ORFModel, ORFRecord, ReferenceContext as ORFContext, EvidenceSummary, ORFEvidence
    from pyensembl import Genome
    from pyensembl.species import find_species_by_name

    example = request.param
    fields = dict(example["model"])
    fields["reference_context"] = ORFContext(**fields["reference_context"])
    record = ORFRecord.from_model(ORFModel(**fields), ORFEvidence(
        summaries=tuple(EvidenceSummary.from_dict(s) for s in example.get("summaries", ()))))
    rna = example["transcript"]
    attrs = (f'gene_id "{record.gene_id}"; gene_name "{record.gene_name}"; transcript_id "{rna["id"]}"; '
             f'gene_biotype "{rna["biotype"]}"; transcript_biotype "{rna["biotype"]}";')
    if rna["version"] is not None:
        attrs += f' transcript_version "{rna["version"]}";'
    rows = [("gene", rna["start"], rna["end"], ""), ("transcript", rna["start"], rna["end"], "")]
    rows += [("exon", a, b, f' exon_id "e{i}"; exon_number "{i + 1}";') for i, (a, b) in enumerate(rna["exons"])]
    gtf, fasta = tmp_path / "slice.gtf", tmp_path / "slice.fa"
    gtf.write_text("".join(f"{record.contig}\tfixture\t{kind}\t{a}\t{b}\t.\t{record.strand}\t.\t{attrs}{extra}\n"
                           for kind, a, b, extra in rows))
    versioned = rna["id"] + ("" if rna["version"] is None else "." + rna["version"])
    fasta.write_text(f">{versioned}\n{rna['sequence']}\n")
    scope = record.reference_context
    genome = Genome(scope.genome_build, scope.annotation_name, scope.annotation_version,
                    gtf_path_or_url=str(gtf), transcript_fasta_paths_or_urls=[str(fasta)],
                    cache_directory_path=str(tmp_path))
    genome.species = find_species_by_name(scope.species)
    genome.assembly_accession, genome.strain, genome.strain_taxon_id = (
        scope.assembly_accession, scope.strain, scope.strain_taxon_id)
    genome.index()
    try:
        yield example, record.bind(genome.transcript_by_id(rna["id"]))
    finally:
        genome.close()


def test_native_reference_and_modified_rna_products_keep_the_code(native_orf):
    from orfdata import ORFCatalog
    from isovar import annotate_sv_orf_evidence
    from varcode.effect_hypotheses import RealizedTranscriptProduct
    from varcode.orf import annotate_orf_rna_product, annotate_orf_variants

    example, bound = native_orf
    record, transcript = bound.record, bound.transcript
    assert translate_cdna(bound.coding_sequence, first_codon_is_start=True,
                          translation_table=record.translation_table) == (record.protein_sequence, True)
    # Synthetic synonymous change in a native context; this is not new RNA evidence.
    start, end = bound.cdna_interval
    sequence = str(transcript.sequence)
    code = bound.genetic_code.codons
    offset, alternate = next((i, alt) for i in range(start + 3, end - 3, 3)
                             for alt in sorted(code) if alt != sequence[i:i + 3]
                             and code[alt] == code[sequence[i:i + 3]])
    mutated = sequence[:offset] + alternate + sequence[offset + 3:]
    intervals = sorted(transcript.exon_intervals, reverse=transcript.strand == "-")
    positions = [p for a, b in intervals for p in (range(a, b + 1) if transcript.strand == "+" else range(b, a - 1, -1))]
    product = RealizedTranscriptProduct(transcript, mutated, None, evidence=dict(
        base_origins=[(record.contig, p, record.strand) for p in positions]))
    annotated = annotate_orf_rna_product(bound, product)
    assert annotated.product.protein_sequence == record.protein_sequence
    assert annotated.reference.model.translation_table == record.translation_table
    assert annotated.product.reference_evidence["mutant_translation_evidence"] == "not_assessed"
    # The same edit also passes Varcode's variant annotation path with the base code.
    from varcode.nucleotides import reverse_complement
    ref, alt = sequence[offset:offset + 3], alternate
    if record.strand == "-":
        ref, alt = reverse_complement(ref), reverse_complement(alt)
    variant = Variant(record.contig, min(positions[offset:offset + 3]), ref, alt, transcript.genome)
    # This fixture intentionally has transcript/exon features only; annotate
    # the supplied ORF directly, without a separate canonical CDS prediction.
    consequence = annotate_orf_variants([variant], bound)
    assert consequence.reference.model.translation_table == record.translation_table
    assert consequence.outcomes[0].mutant.protein_sequence == record.protein_sequence
    # Exact native projection, with synthetic RNA bases, through Isovar's consumer API.
    blocks, cursor = [], 0
    for a, b in intervals:
        blocks.append(dict(query_start=cursor, query_end=cursor + b - a + 1,
                           reference_start=a - 1, reference_end=b, contig=record.contig, strand=record.strand))
        cursor += b - a + 1
    result = dict(reference_name=record.genome_build, event_id="synthetic", paths=[dict(
        path_id="native-path", sequence=mutated, blocks=blocks, exploratory_orfs=dict(candidates=[dict(
            query_interval=[start, end], amino_acids=record.protein_sequence,
            full_interval_support={"fragments": 0}, initiation_observed=False, translation_observed=False)]))])
    catalog = ORFCatalog([record], dict(genome_build=record.genome_build,
                                      taxon_id=record.taxon_id, species=record.species))
    exported = annotate_sv_orf_evidence(result, [transcript], catalog=catalog)
    annotation = exported["paths"][0]["exploratory_orfs"]["candidates"][0]["reference_orf_annotations"][0]
    assert annotation["model"]["translation_table"] == record.translation_table
    assert annotation["product"]["protein_sequence"] == record.protein_sequence
    assert annotation["product"]["reference_evidence"]["mutant_translation_evidence"] == "not_assessed"
    assert tuple(tuple(b) for b in annotation["model"]["coding_blocks"]) == record.coding_blocks
    with pytest.raises(ValueError, match="mismatch|initiation codon"):
        replace(record, translation_table=1).bind(transcript)
