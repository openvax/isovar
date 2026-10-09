"""Partial context is uncertainty, and initiation priors never erase evidence."""

from dataclasses import replace
import gzip
import json

import pytest

from isovar import InitiationContext, UORFCatalog, UORFRecord, annotate_initiation_context, load_uorf_catalog
from isovar.cli.isovar_uorf_evidence import run
from examples.build_uorf_catalog import load_transcript_reference, model_initiation_context


def test_partial_and_ambiguous_positions_are_not_mismatches():
    context = annotate_initiation_context("ACCATG", 3)
    assert context.assessment["key_positions_observed"] == 1
    assert context.assessment["key_positions_matching"] == 1
    assert context.assessment["interpretation"] == "one_key_preference_present"
    assert context.assessment["kozak_positions"][-1]["matches"] is None
    unknown = annotate_initiation_context("NNNATGN", 3)
    assert unknown.assessment["key_positions_observed"] == 0
    assert unknown.assessment["key_positions_matching"] == 0
    mismatch = annotate_initiation_context("TTTTTTTTTATGCCCCCC", 9)
    assert mismatch.assessment["key_positions_observed"] == 2
    assert mismatch.assessment["key_positions_matching"] == 0
    assert "initiation_still_possible" in mismatch.assessment["interpretation"]


def test_atg_at_sequence_start_does_not_prove_a_zero_length_leader():
    context = annotate_initiation_context("AUGG", 0)
    assert (context.upstream, context.start_codon, context.downstream) == ("", "ATG", "G")
    assert context.assessment["is_atg"]
    assert context.assessment["key_positions_matching"] == 1
    assert context.assessment["cap_distance_nt"] is None
    assert "true_5prime_end_unverified" in context.assessment["sequence_boundary_note"]
    verified = annotate_initiation_context("ATGG", 0, five_prime_complete=True)
    assert verified.assessment["cap_distance_nt"] == 0
    assert "conventional_scanning_uncertain" in verified.assessment["sequence_boundary_note"]
    assert annotate_initiation_context("CTGG", 0).assessment["is_atg"] is False


@pytest.mark.parametrize("sequence,start", [("AT", 0), ("ATG", True), ("ATG", -1), ("ATG", 1), ("AT?", 0)])
def test_invalid_sequences_and_offsets_are_rejected(sequence, start):
    with pytest.raises(ValueError):
        annotate_initiation_context(sequence, start)


def test_native_tpst1_context_and_dual_role_snv():
    catalog = load_uorf_catalog()
    native = catalog.get("c7riboseqorf56")
    context = native.initiation_context
    assert context.assessment["display"] == "GAGGCCAGG[ATG]CCGTCC"
    assert context.source == "ensembl101" and context.source_transcript_id.startswith("ENST00000304842.")
    assert context.assessment["key_positions_matching"] == 1
    upstream = native.annotate_initiation_snv(66205480, "A", "C", genome_build="GRCh38")
    assert upstream["context_position"] == -3 and upstream["coding_annotation"] is None
    assert upstream["sequence_changes"] == ["kozak_preference_lost"]
    assert upstream["alternate_context"]["assessment"]["key_positions_matching"] == 0
    body = native.annotate_initiation_snv(66205486, "C", "G", genome_build="GRCh38")
    assert body["context_position"] == 4
    assert body["roles"] == ["initiation_context", "kozak_preference_position", "body"]
    assert body["coding_annotation"]["amino_acid_position"] == 2
    assert body["codon_consequence"]["reference_codon_translation"] == "P"
    assert body["codon_consequence"]["alternate_codon_translation"] == "A"
    assert body["sequence_changes"] == ["kozak_preference_gained"]
    assert body["translation_change"] == "not_measured"
    start = native.annotate_initiation_snv(66205483, "A", "T", genome_build="GRCh38")
    assert start["sequence_changes"] == ["ATG_lost"]
    assert start["alternate_context"]["start_codon"] == "TTG"
    assert native.has_evidence("ribosome") and native.has_evidence("hla")
    assert catalog.query(contig="7", start=66205480, end=66205480, genome_build="GRCh38") == ()
    assert catalog.query(contig="7", start=66205480, end=66205480, genome_build="GRCh38",
                         include_initiation_context=True) == (native,)


@pytest.mark.parametrize("position,ref,alt,build,match", [
    (66205480, "G", "C", "GRCh38", "REF mismatch"),
    (66205480, "A", "C", "GRCh37", "genome build"),
    (66210000, "A", "C", "GRCh38", "unavailable"),
    (66205480, "A", "AA", "GRCh38", "single-base"),
    (66205480, "A", "A", "GRCh38", "distinct")])
def test_snv_validation(position, ref, alt, build, match):
    with pytest.raises(ValueError, match=match):
        load_uorf_catalog().get("c7riboseqorf56").annotate_initiation_snv(position, ref, alt, genome_build=build)


def test_minus_strand_snv_and_spliced_context_do_not_include_introns():
    context = InitiationContext("GCCGCCACC", "ATG", "AAATAA",
                                genomic_positions=tuple(range(311, 299, -1)) + tuple(range(105, 99, -1)))
    record = UORFRecord("minus", "g", "G", "tx", "1", "-", ((100, 105), (300, 302)),
                        "MK", "uORF", (), (), initiation_context=context)
    catalog = UORFCatalog([record], {"genome_build": "GRCh38"})
    effect = record.annotate_initiation_snv(105, "T", "C", genome_build="GRCh38")
    assert effect["context_position"] == 4 and effect["variant"]["transcript_alt"] == "G"
    assert effect["codon_consequence"]["alternate_codon_translation"] == "E"
    assert effect["sequence_changes"] == ["kozak_preference_gained"]
    assert record.map_initiation_position(200, genome_build="GRCh38") is None
    assert catalog.query(contig="1", start=200, end=200, genome_build="GRCh38", include_initiation_context=True) == ()
    with pytest.raises(ValueError, match="strand"):
        UORFCatalog([replace(record, initiation_context=replace(context, genomic_positions=context.genomic_positions[::-1]))], {"genome_build": "GRCh38"})
    bad = tuple(range(312, 300, -1)) + tuple(range(105, 99, -1))
    with pytest.raises(ValueError, match="coding path"):
        UORFCatalog([replace(record, initiation_context=replace(context, genomic_positions=bad))], {"genome_build": "GRCh38"})


def test_reference_context_requires_entire_splice_path_and_exact_protein():
    model = dict(contig="1", strand="+", coding_blocks=((100, 101), (200, 206)), protein_sequence="MK")
    reference = dict(contig="1", strand="+", exons=((95, 101), (200, 206)),
                     sequence="GCCACATGAAATAA", transcript_id="tx.1")
    context = model_initiation_context(model, reference)
    assert context.start_codon == "ATG" and context.genomic_positions == tuple(range(95, 102)) + tuple(range(200, 207))
    assert context.upstream == "GCCAC" and context.assessment["status"] == "partial"
    assert model_initiation_context(model, None).unavailable_reason == "transcript_not_in_reference"
    assert model_initiation_context(dict(model, protein_sequence="MP"), reference).unavailable_reason == "ORF_protein_mismatch"
    assert model_initiation_context(dict(model, coding_blocks=((100, 102), (201, 206))), reference).unavailable_reason == "ORF_splice_path_mismatch"
    assert model_initiation_context(model, dict(reference, sequence="ATG")).unavailable_reason == "transcript_sequence_length_mismatch"


def test_reference_loading_preserves_version_and_checks_disagreement(tmp_path):
    gtf, fasta = tmp_path / "x.gtf.gz", tmp_path / "x.fa.gz"
    with gzip.open(gtf, "wt") as stream:
        stream.write('1\tX\texon\t1\t9\t.\t+\t.\ttranscript_id "tx"; transcript_version "2";\n')
    with gzip.open(fasta, "wt") as stream:
        stream.write(">tx.2 description\nATGAAATAA\n>other.1\nATG\n")
    reference = load_transcript_reference(gtf, fasta, {"tx"})
    assert reference["tx"]["transcript_id"] == "tx.2"
    assert reference["tx"]["sequence"] == "ATGAAATAA" and "other" not in reference
    with gzip.open(fasta, "wt") as stream:
        stream.write(">tx.3\nATGAAATAA\n")
    with pytest.raises(ValueError, match="versions"):
        load_transcript_reference(gtf, fasta, {"tx"})


def test_cli_context_lookup_and_ref_errors(capsys):
    run(["--position", "7:66205480", "--genome-build", "GRCh38", "--ref", "A", "--alt", "C"])
    output = json.loads(capsys.readouterr().out)
    assert output["position_annotations"] == [None]
    assert output["initiation_annotations"][0]["sequence_changes"] == ["kozak_preference_lost"]
    with pytest.raises(SystemExit) as exc:
        run(["--position", "7:66205480", "--genome-build", "GRCh38", "--ref", "T", "--alt", "C"])
    assert exc.value.code == 2 and "REF mismatch" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        run(["--ref", "A", "--alt", "C"])
    with pytest.raises(SystemExit):
        run(["--position", "7:66240380", "--genome-build", "GRCh38", "--ref", "A", "--alt", "C"])
    assert "No mapped initiation context" in capsys.readouterr().err


def test_old_catalogues_and_new_export_roundtrip(tmp_path):
    catalog = load_uorf_catalog()
    record = catalog.get("c7riboseqorf56").to_dict()
    data = dict(schema="isovar.uorf_evidence.v1", metadata=catalog.metadata, records=[record])
    path = tmp_path / "reference.json"
    path.write_text(json.dumps(data))
    assert load_uorf_catalog(path).records == (catalog.get("c7riboseqorf56"),)
    del record["initiation_context"]
    path.write_text(json.dumps(data))
    assert load_uorf_catalog(path).records[0].initiation_context is None
