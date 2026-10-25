"""The ORFdata adapter keeps Isovar products and RNA provenance offline."""

from copy import deepcopy
import json

from orfdata import ProducerInfo, Proteome, ProteomeEntry, ReferenceContext
import pytest

from isovar import (
    ProteinSequenceCreator, export_protein_hypotheses, from_protein_hypotheses, run_isovar, to_proteome,
)

from .test_protein_hypotheses import BAM, VCF, synthetic_result


REFERENCE = ReferenceContext(9606, "human", "GRCh38", "Ensembl", "87")


def exported():
    return export_protein_hypotheses(
        [synthetic_result(settings={"max_protein_sequences_per_variant": 0})],
        sample_id="tumor-1", source="tumor-rna.bam")


def records_of(proteome, kind):
    return [record for record in proteome.records if record.kind == kind]


def test_native_and_offline_conversion_preserve_products_and_evidence(tmp_path):
    export = exported()
    offline = from_protein_hypotheses(json.loads(json.dumps(export)), reference=REFERENCE)
    native = to_proteome(
        [synthetic_result(settings={"max_protein_sequences_per_variant": 0})],
        sample_id="tumor-1", source="tumor-rna.bam", reference=REFERENCE)
    assert native.identity == offline.identity
    assert [entry.sequence for entry in native.entries] == ["MAQG", "MAQD", "AQG"]
    assert len(records_of(native, "translation")) == 4
    assert len(records_of(native, "evidence")) == len(export["evidence_sets"])
    assert {record.native_id: record.to_dict()["payload"]
            for record in records_of(native, "evidence")} == export["evidence_sets"]
    translations = {translation["translation_id"]: translation
                    for protein in export["events"][0]["protein_hypotheses"]
                    for translation in protein["translations"]}
    assert {record.payload["translation_id"]: record.to_dict()["payload"]
            for record in records_of(native, "translation")} == translations
    assert records_of(native, "event")[0].payload["filters"] == export["events"][0]["filters"]
    assert {record.payload["nucleotide_sequence"] for record in records_of(native, "translation")
            if record.payload["nucleotide_sequence"].startswith("CCATGGC")} >= {
                "CCATGGCTCAAGGCTAA", "CCATGGCCCAAGGCTAA"}
    assert native.entries[0].producer.start_complete is True
    assert native.entries[0].producer.stop_complete is True
    assert native.entries[2].producer.start_complete is False
    assert native.entries[2].producer.stop_complete is False
    assert native.entries[0].producer.gene_names == ("GENE1",)
    assert any(link.relation == "supports_rna_sequence" for link in native.links)
    assert all(record.support == "rna_sequence" for record in records_of(native, "evidence"))
    assert native.merge(offline).identity == native.identity

    path = tmp_path / "proteome.json.gz"
    native.save(path)
    loaded = Proteome.load(path)
    assert loaded.identity == native.identity
    fasta = tmp_path / "proteome.fasta"
    loaded.write_fasta(fasta)
    headers = [line for line in fasta.read_text().splitlines() if line.startswith(">")]
    assert len(headers) == 3
    assert all("isovar" in header.split()[0] and "GENE1" in header.split()[0] for header in headers)


def test_zero_product_unknown_and_capped_events_remain_distinct():
    export = exported()
    zero = deepcopy(export["events"][0])
    zero["protein_hypotheses"] = []
    zero["protein_hypotheses_complete"] = True
    unknown = deepcopy(zero)
    unknown["protein_hypotheses_complete"] = None
    capped = deepcopy(export["events"][0])
    capped["protein_hypotheses_complete"] = False
    capped["protein_sequence_limit"] = 3
    export["events"] = [zero, unknown, capped]
    proteome = from_protein_hypotheses(export, reference=REFERENCE)
    assert len(proteome.entries) == 3
    assert [record.status for record in records_of(proteome, "event")] == [
        "zero_product", "completeness_unknown", "capped"]
    assert records_of(proteome, "event")[-1].limits["max_protein_sequences_per_variant"] == 3
    assert [record.payload["event_id"] for record in records_of(proteome, "event")] == [
        zero["event_id"]] * 3


def test_reference_and_schema_must_be_explicit_and_consistent():
    export = exported()
    with pytest.raises(ValueError, match="Expected isovar.protein_hypotheses.v2"):
        from_protein_hypotheses({"schema": "isovar.sv_rna_orfs.v4"}, reference=REFERENCE)
    with pytest.raises(ValueError, match="conflicts"):
        from_protein_hypotheses(export, reference=ReferenceContext(
            9606, "human", "GRCh37", "Ensembl", "75"))
    with pytest.raises(ValueError, match="ReferenceContext"):
        from_protein_hypotheses(export, reference=None)
    broken = deepcopy(export)
    broken["events"][0]["protein_hypotheses"][0]["rna_support"]["evidence_set_id"] = "missing"
    with pytest.raises(ValueError, match="Unknown RNA evidence"):
        from_protein_hypotheses(broken, reference=REFERENCE)


def test_read_set_and_gene_names_remain_source_qualified_and_escaped():
    first = exported()
    first["sample_id"] = "tumor A"
    first["source"] = "reads / A.bam"
    first["evidence_scope"] = ["tumor A", "reads / A.bam"]
    first["events"][0]["protein_hypotheses"][0]["gene_names"] = ["GENE A/1", "GENE B"]
    other = deepcopy(first)
    other["source"] = "reads / B.bam"
    other["evidence_scope"] = ["tumor A", "reads / B.bam"]
    a = from_protein_hypotheses(first, reference=REFERENCE)
    b = from_protein_hypotheses(other, reference=REFERENCE)
    joint = a.merge(b)
    assert len(joint.entries) == 6
    assert a.entries[0].producer.gene_names == ("GENE A/1", "GENE B")
    assert a.entries[0].fasta_token() != b.entries[0].fasta_token()
    assert " " not in a.entries[0].fasta_token()
    assert len(joint.by_sequence()["MAQG"]) == 2


def test_unknown_gene_and_mixed_start_state_stay_explicit():
    export = exported()
    protein = export["events"][0]["protein_hypotheses"][0]
    protein["gene_ids"] = []
    protein["gene_names"] = []
    protein["n_terminus"] = "mixed"
    converted = from_protein_hypotheses(export, reference=REFERENCE)
    assert converted.entries[0].producer.gene_names == ()
    assert converted.entries[0].producer.start_complete is None
    assert "gene=unknown" in converted.entries[0].fasta_token()
    hypothesis = records_of(converted, "hypothesis")[0]
    assert hypothesis.payload["n_terminus"] == "mixed"


def test_reference_ncorf_and_isovar_sources_compose_without_sequence_dedup(tmp_path):
    isovar = from_protein_hypotheses(exported(), reference=REFERENCE)
    reference = Proteome(REFERENCE, [ProteomeEntry("ENSP1", "MAQG", "annotated", 9606,
                                                  "Ensembl:87")])
    ncorf = ProteomeEntry.from_producer("MAQG", taxon_id=9606, producer=ProducerInfo(
        "nc_orf", "catalogue", "orfdata", "ORF1", reference=REFERENCE))
    varcode = ProteomeEntry.from_producer("MAQG", taxon_id=9606, producer=ProducerInfo(
        "reference_cds", "dna_prediction", "varcode", "V1", sample_id="tumor-1",
        reference=REFERENCE))
    joint = reference.merge(Proteome.from_entries((ncorf, varcode), reference=REFERENCE), isovar)
    assert len(joint.by_sequence()["MAQG"]) == 4
    joint.write_fasta(tmp_path / "joint.fasta")
    assert (tmp_path / "joint.fasta").read_text().count(">") == len(joint.entries)
    grouped = joint.fasta_index(deduplicate=True)
    assert any(len(members) == 4 and {member.protein_source for member in members} == {
        "reference_cds", "nc_orf"} for members in grouped.values())
    path = tmp_path / "joint.proteome.json.gz"
    joint.save(path)
    assert len(Proteome.load(path).by_sequence()["MAQG"]) == 4


def test_small_real_b16_reads_convert_without_live_reference(tmp_path):
    results = run_isovar(
        VCF, BAM, protein_sequence_creator=ProteinSequenceCreator(max_protein_sequences_per_variant=0))
    reference = ReferenceContext(10090, "mouse", "GRCm38", "Ensembl", None)
    export = export_protein_hypotheses(results, sample_id="b16", source="b16.bam")
    native = to_proteome(results, sample_id="b16", source="b16.bam", reference=reference)
    offline = from_protein_hypotheses(json.loads(json.dumps(export)), reference=reference)
    assert native.identity == offline.identity
    assert len(native.entries) == sum(len(event["protein_hypotheses"]) for event in export["events"])
    assert len(records_of(native, "event")) == 4
    assert any(record.status == "zero_product" for record in records_of(native, "event"))
    path = tmp_path / "b16.proteome.json.gz"
    native.save(path)
    assert Proteome.load(path).identity == native.identity
