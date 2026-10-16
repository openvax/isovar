"""Observed viral coding predictions use ORFdata and preserve sample evidence."""

import gzip
from hashlib import md5, sha256
import json
from pathlib import Path

import pysam
import pytest
from orfdata import ReferenceContext, SourceRNA, TranscriptORF, viral_orfs

from isovar import write_viral_proteins
from isovar.cli.commands import run
from isovar.dna import reverse_complement_dna
from isovar.viral_coding import viral_coding_predictions
from tests.test_viral_sequences import read, reconstruct


@pytest.fixture
def hpv_inputs(tmp_path):
    text = gzip.decompress((Path(__file__).parent / "data/viral/hpv16-GDC.fa.gz").read_bytes()).decode()
    dna = "".join(text.splitlines()[1:])
    assert len(dna) == 7906 and md5(dna.encode()).hexdigest() == "f2463721818c29998ef5ab917e6df74a"
    sequences = {"HPV16": dna, "host": "ACGT" * 250}
    fasta = tmp_path / "reference.fa"
    fasta.write_text("".join(">%s\n%s\n" % row for row in sequences.items()))
    header = pysam.AlignmentHeader.from_dict(dict(
        SQ=[dict(SN=name, LN=len(seq), M5=md5(seq.encode()).hexdigest()) for name, seq in sequences.items()],
        RG=[dict(ID="a", SM="sample")]))
    return sequences, fasta, header


def e7_prediction(row):
    return next(c for c in row["coding"] if any(r["protein_id"] == "NP_041326.1" for r in c["reference_orfs"]))


def test_complete_sample_e7_variant_uses_observed_rna_and_original_support(tmp_path, hpv_inputs):
    dna, _, _ = hpv_inputs
    observed = dna["HPV16"][538:882]
    assert observed[645 - 539] == "A"
    observed = observed[:645 - 539] + "C" + observed[645 - 539 + 1:]
    result = reconstruct(tmp_path, hpv_inputs,
                         [read(hpv_inputs, "observed-E7", "HPV16", 538, sequence=observed, length=344)],
                         contigs=["HPV16"])
    row, = result["sequences"]
    coding = e7_prediction(row)
    assert row["sequence"] == observed and row["rna_strand"] == "unknown"
    assert (coding["query_start"], coding["query_end"]) == (23, 320)
    assert coding["status"] == "complete_CDS"
    assert len(coding["amino_acids"]) == 98 and coding["amino_acids"][27] == "F"
    source, = coding["reference_orfs"]
    assert source["reference_amino_acids"][27] == "L"
    assert source["protein_matches_reference"] is False
    assert source["coding_path"] == "compatible_at_observed_positions"
    assert source["start_status"] == "observed_source_start" and source["stop_status"] == "source_stop"
    assert coding["all_bases_observed"] and len(coding["whole_cds_segments"]) == 1
    assert row["compatible_fragments"] == result["support"]["target_fragments"] == 1
    assert result["coding_reference"]["provider"] == "orfdata"
    coverage = result["coding_reference"]["contigs"][0]
    assert coverage["coding_models"] == 9
    assert coverage["reference"]["id"] == "GRCh38.d1.vd1:HPV16"
    assert coverage["reference"]["virus"]["type"] == "HPV16"
    assert coverage["reference"]["host_taxon_ids"] == [9606]
    assert "dna" not in coverage["reference"] and "orfs" not in coverage["reference"]
    json.dumps(result["coding_reference"])
    output = tmp_path / "proteins.fa"
    write_viral_proteins(result, output)
    with pysam.FastxFile(str(output)) as records:
        assert any(r.sequence == coding["amino_acids"] and "status=complete_CDS" in r.comment for r in records)


def test_partial_coding_path_keeps_missing_start_and_stop(tmp_path, hpv_inputs):
    result = reconstruct(tmp_path, hpv_inputs, [read(hpv_inputs, "partial", "HPV16", 582, length=60)],
                         contigs=["HPV16"])
    coding = e7_prediction(result["sequences"][0])
    assert coding["status"] == "partial_CDS"
    assert not coding["complete_start"] and not coding["observed_stop"]
    assert coding["amino_acids"].startswith("L")  # ordinary Leu, not an invented initiating Met
    assert len(coding["amino_acids"]) == 20
    assert coding["whole_cds_segments"] == []


@pytest.mark.parametrize("protein_id,blocks,length", [
    ("YP_009268707.1", ((83, 226), (409, 417)), 50),
    ("YP_009268708.1", ((865, 880), (3358, 3620)), 92),
])
def test_source_joined_products_preserve_phase_and_separate_other_coding_paths(
        tmp_path, hpv_inputs, protein_id, blocks, length):
    dna = hpv_inputs[0]["HPV16"]
    observed = "".join(dna[a - 1:b] for a, b in blocks)
    first, second = blocks
    cigar = "%dM%dN%dM" % (first[1] - first[0] + 1, second[0] - first[1] - 1,
                            second[1] - second[0] + 1)
    result = reconstruct(tmp_path, hpv_inputs,
                         [read(hpv_inputs, "joined", "HPV16", first[0] - 1,
                               sequence=observed, cigar=cigar)], contigs=["HPV16"])
    predictions = [c for c in result["sequences"][0]["coding"]
                   if any(r["protein_id"] == protein_id for r in c["reference_orfs"])]
    coding, = predictions
    reference = next(r for r in coding["reference_orfs"] if r["protein_id"] == protein_id)
    assert len(coding["amino_acids"]) == length
    assert coding["status"] == "complete_CDS" and reference["protein_matches_reference"] is True
    assert reference["coding_path"] == "compatible_at_observed_positions"
    assert any(r["coding_path"] == "altered" for r in coding["reference_orfs"])


def test_reverse_oriented_path_translates_without_changing_rna_coordinates(hpv_inputs):
    dna = hpv_inputs[0]["HPV16"]
    observed = dna[538:882]
    positions = [("HPV16", p, "-") for p in reversed(range(538, 882))]
    rows = viral_coding_predictions("reverse", reverse_complement_dna(observed), positions,
                                    viral_orfs(dna, contig="HPV16"))
    coding = e7_prediction(dict(coding=rows))
    assert coding["orientation"] == "-" and coding["status"] == "complete_CDS"
    assert (coding["query_start"], coding["query_end"]) == (24, 321)
    assert coding["reference_orfs"][0]["protein_matches_reference"] is True


def synthetic_model():
    coding = "ATGAAACCCGGGTTTTAA"
    context = ReferenceContext(999999, "synthetic_fixture", "fixture", "fixture", "1",
                               sequence_sha256=sha256(coding.encode()).hexdigest())
    rna = SourceRNA("fixture-CDS", coding, context, "synthetic coding fixture",
                    contig="fixture", strand="+", exons=((1, 18),))
    return TranscriptORF.at(rna, 0, orf_id="fixture-protein", source=rna.source)


def test_observed_insertion_changes_translation_and_never_fills_the_stop():
    model = synthetic_model()
    sequence = model.coding_sequence[:6] + "A" + model.coding_sequence[6:]
    positions = [("fixture", p, "+") for p in range(6)] + [None] + [
        ("fixture", p, "+") for p in range(6, 18)]
    coding, = viral_coding_predictions("inserted", sequence, positions, [model])
    assert coding["amino_acids"] != model.protein_sequence
    assert coding["coding_sequence"] == sequence[:18]
    assert not coding["observed_stop"] and coding["status"] == "partial_CDS"
    assert coding["reference_orfs"][0]["coding_path"] == "altered"


def test_lost_source_start_is_not_promoted_to_an_initiation_claim():
    model = synthetic_model()
    sequence = "GTC" + model.coding_sequence[3:]
    positions = [("fixture", p, "+") for p in range(len(sequence))]
    coding, = viral_coding_predictions("lost-start", sequence, positions, [model])
    assert not coding["complete_start"] and coding["amino_acids"].startswith("V")
    assert coding["reference_orfs"][0]["start_status"] == "source_start_altered"


def test_observed_early_stop_is_distinct_from_source_termination():
    model = synthetic_model()
    sequence = model.coding_sequence[:6] + "TAA" + model.coding_sequence[9:]
    coding, = viral_coding_predictions("early-stop", sequence,
        [("fixture", p, "+") for p in range(len(sequence))], [model])
    assert coding["amino_acids"] == "MK" and coding["observed_stop"]
    assert coding["reference_orfs"][0]["stop_status"] == "alternative_stop"
    assert coding["reference_orfs"][0]["protein_matches_reference"] is False


def test_same_name_unqualified_reference_retains_rna_with_explicit_coding_gap(tmp_path, hpv_inputs):
    dna, fasta, _ = hpv_inputs
    changed = "A" if dna["HPV16"][0] != "A" else "C"
    dna["HPV16"] = changed + dna["HPV16"][1:]
    fasta.write_text("".join(">%s\n%s\n" % row for row in dna.items()))
    header = pysam.AlignmentHeader.from_dict(dict(
        SQ=[dict(SN=name, LN=len(seq), M5=md5(seq.encode()).hexdigest()) for name, seq in dna.items()],
        RG=[dict(ID="a", SM="sample")]))
    inputs = dna, fasta, header
    result = reconstruct(tmp_path, inputs, [read(inputs, "still-RNA", "HPV16", 582, length=60)],
                         contigs=["HPV16"])
    assert result["sequences"] and not result["sequences"][0]["coding"]
    assert result["sequences"][0]["coding_status"] == "unavailable_reference"
    assert result["coding_reference"]["contigs"][0]["reason"] == "no_sequence_qualified_ORFdata_reference"


def test_command_exports_predicted_proteins_without_a_gtf(tmp_path, hpv_inputs):
    reconstruct(tmp_path, hpv_inputs, [read(hpv_inputs, "E7", "HPV16", 538, length=344)], contigs=["HPV16"])
    output = tmp_path / "out"
    run(["viral", "--bam", str(tmp_path / "rna.bam"), "--genome-fasta", str(hpv_inputs[1]),
         "--sample-id", "sample", "--contig", "HPV16", "--alignment-reference", "custom",
         "--reference-cache", str(tmp_path / "cache"), "--output-dir", str(output)])
    result = json.loads((output / "viral.json").read_text())
    assert e7_prediction(result["sequences"][0])["status"] == "complete_CDS"
    assert "status=complete_CDS" in (output / "proteins.fa").read_text()
    assert (output / "sequences.fa").is_file()
