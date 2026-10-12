"""Synthetic RNA products exercise reference scope, not experimental evidence."""

import json
from dataclasses import replace

import pytest
from pyensembl import Genome
from pyensembl.species import find_species_by_name

from isovar import annotate_sv_orf_evidence
from orfdata import ORFCatalog, ORFModel, ORFRecord, ReferenceContext


@pytest.mark.parametrize("taxon,species", [
    (9606, "homo_sapiens"), (10090, "mus_musculus"),
    (9615, "canis_lupus_familiaris"), (4932, "saccharomyces_cerevisiae"),
])
def test_versioned_overlay_retains_each_species_context(tmp_path, taxon, species):
    attrs = 'gene_id "G"; gene_name "G"; transcript_id "T"; transcript_version "2"; gene_biotype "lncRNA"; transcript_biotype "lncRNA";'
    gtf = tmp_path / "slice.gtf"
    gtf.write_text("".join(f'1\ttest\t{kind}\t100\t108\t.\t+\t.\t{attrs}{extra}\n'
                           for kind, extra in (("gene", ""), ("transcript", ""),
                                               ("exon", ' exon_id "E"; exon_number "1";'))))
    fasta = tmp_path / "slice.fa"
    fasta.write_text(">T\nATGAAATAA\n")
    genome = Genome("shared_label", "synthetic", "1", gtf_path_or_url=str(gtf),
                    transcript_fasta_paths_or_urls=[str(fasta)], cache_directory_path=str(tmp_path))
    genome.species = find_species_by_name(species)
    genome.index()
    try:
        transcript = genome.transcript_by_id("T")
        record = ORFRecord.from_model(ORFModel(
            "shared_orf", "T.2", "G", "G", "1", "+", ((100, 108),), "MK", "ncRNA",
            "shared_label", taxon, species, source="synthetic:1"))
        metadata = dict(genome_build="shared_label", taxon_id=taxon, species=species)
        catalog = ORFCatalog([record], metadata)
        result = dict(reference_name="shared_label", event_id="synthetic-event", paths=[dict(
            path_id="path", sequence="ATGAAATAA",
            blocks=[dict(query_start=0, query_end=9, reference_start=99, reference_end=108,
                         contig="1", strand="+")],
            exploratory_orfs=dict(candidates=[dict(query_interval=[0, 9], amino_acids="MK",
                full_interval_support={"fragments": 2}, initiation_observed=False,
                translation_observed=False)]),
        )])

        def annotation(**kwargs):
            output = annotate_sv_orf_evidence(result, [transcript], catalog=catalog, **kwargs)
            return output["paths"][0]["exploratory_orfs"]["candidates"][0]["reference_orf_annotations"][0]

        observed = annotation()
        assert observed["status"] == "conditional_translation"
        assert observed["reference"]["reference_context"]["taxon_id"] == taxon
        assert observed["model"]["transcript_id"] == "T.2"
        assert observed["product"]["evidence"]["rna_support"] == {"fragments": 2}
        assert observed["product"]["reference_evidence"]["mutant_translation_evidence"] == "not_assessed"
        assert json.loads(json.dumps(observed))["reference"]["reference_key"] == ReferenceContext.from_genome(genome).key
        with pytest.raises(ValueError, match="Ambiguous supplied transcript"):
            annotate_sv_orf_evidence(result, [transcript, transcript], catalog=catalog)
        catalog = ORFCatalog([replace(record, transcript_id="T.3")], metadata)
        assert "version mismatch" in annotation()["reason"]
        catalog = ORFCatalog([record], metadata)
        explicit = ReferenceContext.from_genome(genome)
        del genome.species
        assert "explicit ReferenceContext" in annotation()["reason"]
        assert annotation(reference_context=explicit)["status"] == "conditional_translation"
    finally:
        genome.close()
