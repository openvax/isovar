"""Original tumor RNA -> transcript matching -> ranked mutant protein windows."""

import copy
import gzip
import json
from pathlib import Path

import pysam
import pytest
from varcode import Variant

from isovar.protein_sequence_creator import ProteinSequenceCreator
from isovar.read_collector import ReadCollector

from .data.osteosarc import audit
from .osteosarc_protein_helpers import (
    REFERENCE, check_translation, expected_variant, load_references, reference_genome, translate,
)


DATA = Path(__file__).parent / "data" / "osteosarc"
SELECTION = json.loads((DATA / "selection.json").read_text())
MANIFEST = json.loads((DATA / "manifest.json").read_text())


@pytest.fixture(scope="module")
def protein_cases(tmp_path_factory):
    directory = tmp_path_factory.mktemp("osteosarc_proteins")
    genome = reference_genome(directory / "reference")
    references = load_references()
    cases = {}
    for name, metadata in MANIFEST["datasets"].items():
        sam = directory / (name + ".sam")
        sam.write_bytes(gzip.decompress((DATA / metadata["file"]).read_bytes()))
        bam_path = directory / (name + ".bam")
        pysam.sort("--no-PG", "-o", str(bam_path), str(sam))
        pysam.index(str(bam_path))
        for record in SELECTION["variants"]:
            variant = Variant(record["chrom"].removeprefix("chr"), int(record["pos"]),
                              record["ref"], record["alt"], ensembl=genome)
            expected = expected_variant(record, references[record["gene"]])
            with pysam.AlignmentFile(bam_path) as bam:
                evidence = ReadCollector().read_evidence_for_variant(variant, bam)
            cases[name, record["gene"]] = variant, evidence, expected
    return cases


@pytest.mark.parametrize("sample", MANIFEST["datasets"])
@pytest.mark.parametrize("gene", [r["gene"] for r in SELECTION["variants"]])
@pytest.mark.parametrize("coverage", [1, 2])
def test_real_rna_ranked_proteins_match_independent_expected_sequences(protein_cases, sample, gene, coverage):
    variant, evidence, expected = protein_cases[sample, gene]
    result = audit.protein_result(variant, evidence, expected, coverage)
    assert result["status"] == "ok", result
    if gene == "MAP2" or (sample == "bulk_star_t0" and gene == "PIP5K1A"):
        assert not evidence.alt_reads
        assert result["outcome"] == "no_alt_reads"
        assert result["top_proteins"] == []
    else:
        assert evidence.alt_reads
        assert result["outcome"] == "expected_top_protein", result
        assert result["expected_matching_translation_count"] > 0
        assert len(result["top_proteins"]) == 1
        top = result["top_proteins"][0]
        assert top["contains_mutation"]
        assert top["frameshift"] == (gene == "PIP5K1A")
        assert top["supporting_read_names"] >= coverage
        if gene == "EXOC4":
            assert top["amino_acids"][top["mutation_interval"][0]] == "I"
        if gene == "DYNC1H1":
            assert top["amino_acids"][top["mutation_interval"][0]] == "I"


@pytest.mark.parametrize("record", SELECTION["variants"], ids=lambda r: r["gene"])
def test_reference_edits_reproduce_reported_protein_changes(record):
    expected = expected_variant(record, load_references()[record["gene"]])
    protein, mutant = expected["protein"], expected["mutant_protein"]
    gene = record["gene"]
    if gene == "EXOC4":
        assert protein[33] == "S"
        assert mutant == protein[:33] + "I" + protein[34:]
    elif gene == "DYNC1H1":
        assert protein[313] == "V"
        assert mutant == protein[:313] + "I" + protein[314:]
    elif gene == "H1-2":
        assert expected["strand"] == "-"
        assert protein[196:201] == "AAKPK"
        assert mutant == protein[:196] + protein[201:]
    elif gene == "GTF3C5":
        # Equivalent deletions within the E repeat yield the same full protein;
        # genomic normalization need not use the rightmost HGVS protein junction.
        assert protein[502:506] == "EEEE"
        assert mutant == protein[:502] + protein[506:]
    else:
        first, mutant_length = (473, 494) if gene == "PIP5K1A" else (866, 906)
        assert protein[first] == ("G" if gene == "PIP5K1A" else "L")
        assert mutant[:first] == protein[:first]
        assert mutant[first:] != protein[first:first + len(mutant) - first]
        assert len(mutant) == mutant_length
    assert expected["mutant_has_stop"]


@pytest.mark.parametrize("attribute,value", [
    ("amino_acids", "WRONG"), ("frameshift", False),
    ("ends_with_stop_codon", True), ("mutation_start_idx", 0),
    ("mutation_end_idx", 0), ("contains_mutation", False),
])
def test_protein_oracle_detects_corrupted_translation(protein_cases, attribute, value):
    variant, evidence, expected = protein_cases["ont_t1", "PIP5K1A"]
    translation = ProteinSequenceCreator(protein_sequence_length=20).translate_variant_reads(
        variant, evidence.alt_reads)[0]
    assert check_translation(translation, expected, 20)["matches_expected"]
    corrupt = copy.copy(translation)
    setattr(corrupt, attribute, value)
    with pytest.raises(AssertionError):
        check_translation(corrupt, expected, 20)


def test_protein_oracle_detects_wrong_reading_frame(protein_cases):
    variant, evidence, expected = protein_cases["ont_t1", "PIP5K1A"]
    translation = ProteinSequenceCreator(protein_sequence_length=20).translate_variant_reads(
        variant, evidence.alt_reads)[0]
    assert check_translation(translation, expected, 20)["matches_expected"]
    corrupt = copy.copy(translation)
    corrupt.variant_orf = copy.copy(corrupt.variant_orf)
    corrupt.variant_orf.offset_to_first_complete_codon = (corrupt.offset_to_first_complete_codon + 1) % 3
    with pytest.raises(AssertionError):
        check_translation(corrupt, expected, 20)


def test_protein_stage_failure_preserves_read_counts(protein_cases, monkeypatch):
    variant, evidence, expected = protein_cases["ont_t1", "PIP5K1A"]
    monkeypatch.setattr(audit.ReadCollector, "read_evidence_for_variant", lambda *a: evidence)

    def fail(*args, **kwargs):
        raise RuntimeError("deliberate translation failure")

    monkeypatch.setattr(audit.ProteinSequenceCreator, "translate_variant_reads", fail)
    result = audit.isovar_result(None, variant, expected)
    assert result["status"] == "ok"
    assert result["counts"]["alt"] >= 3
    assert result["sequence_creation_default"]["status"] == "ok"
    for key in ("protein_default", "protein_coverage1"):
        assert result[key]["status"] == "error"
        assert result[key]["stage"] == "translation"
        assert "top_proteins" not in result[key]


def test_independent_translation_stops_and_ignores_partial_terminal_codon():
    assert translate("ATGTGGTAAATG") == ("MW", True)
    assert translate("ATGTGGA") == ("MW", False)


def test_protein_reference_drift_is_not_an_expected_sequence(tmp_path):
    import shutil

    shutil.copytree(REFERENCE, tmp_path / "reference")
    (tmp_path / "reference" / "reference.gtf.gz").write_bytes(b"corrupted reference")
    with pytest.raises(ValueError, match="Protein reference checksum mismatch"):
        load_references(tmp_path / "reference")


def test_subset_reference_does_not_poison_full_grch38_contig_validation(protein_cases):
    from .genomes_for_testing import grch38

    # A real subset query primes Varcode's contig cache. A subsequent query
    # on a chromosome absent from that subset must still use the full genome.
    variant, _, _ = protein_cases["ont_t1", "PIP5K1A"]
    assert variant.gene_names == ["PIP5K1A"]
    assert variant.reference_name == "GRCh38-osteosarc-six-transcript-subset"
    assert "TP53" in Variant("17", 7676589, "CTC", "", ensembl=grch38).gene_names


def test_full_region_snapshot_validates_proteins_without_calling_no_alt_a_failure():
    report = json.loads((DATA / "sample_audit.json").read_text())
    assert report["schema_version"] == 4
    supported = 0
    for row in report["rows"]:
        for mode in ("primary_only", "defaults"):
            result = row[mode]
            assert result["status"] == "ok"
            for key in ("protein_default", "protein_coverage1", "protein_support20", "protein_balanced95", "protein_balanced90"):
                protein = result[key]
                assert protein["status"] == "ok"
                assert protein["outcome"] == ("expected_top_protein" if result["counts"]["alt"] else "no_alt_reads")
                if key != "protein_support20":
                    assert protein["target_protein_length"] == 49
                    assert protein["peptide_length"] == 25
                    for top in protein["top_proteins"]:
                        assert top["fraction_of_best_candidate_support"] >= protein["min_support_fraction"]
                        assert top["mutation_containing_peptide_windows"] >= 0
                        assert min(top["retained_cdna_min_coverages"]) >= protein["min_coverage"]
            assert result["protein_context"]["status"] == "ok"
            if mode == "primary_only" and result["counts"]["alt"]:
                supported += 1
    assert supported == 9


@pytest.mark.parametrize("sample", MANIFEST["datasets"])
@pytest.mark.parametrize("gene", [r["gene"] for r in SELECTION["variants"]])
@pytest.mark.parametrize("fraction", [0.85, 0.9, 0.95, 1.0])
def test_real_rna_balanced_selection_respects_budget_and_independent_oracle(protein_cases, sample, gene, fraction):
    variant, evidence, expected = protein_cases[sample, gene]
    result = audit.protein_result(variant, evidence, expected, min_protein_sequence_support_fraction=fraction)
    assert result["status"] == "ok", result
    assert result["preference"] == "balanced"
    assert result["target_protein_length"] == 49
    assert result["candidate_context_lengths"] == [20, 25, 29, 33, 37, 41, 45, 49]
    for top in result["top_proteins"]:
        assert top["matches_expected"], result
        assert top["fraction_of_best_candidate_support"] >= fraction
        assert top["supporting_read_names"] <= len(evidence.alt_read_names)


def test_full_region_dync1h1_context_and_support_tradeoffs():
    report = json.loads((DATA / "sample_audit.json").read_text())
    expected = {
        "bulk_star_t0": {"protein_support20": (9, 121, 0), "protein_default": (47, 111, 23),
                         "protein_balanced95": (35, 115, 11), "protein_context": (47, 111, 23)},
        "ont_t1": {"protein_support20": (20, 766, 0), "protein_default": (37, 666, 13),
                   "protein_balanced95": (25, 744, 1), "protein_context": (49, 516, 25)},
    }
    for row in report["rows"]:
        if row["variant"]["gene"] != "DYNC1H1":
            continue
        for key, values in expected[row["sample"]].items():
            top, = row["primary_only"][key]["top_proteins"]
            assert top["matches_expected"]
            assert (len(top["amino_acids"]), top["supporting_read_names"],
                    top["mutation_containing_peptide_windows"]) == values


def test_full_region_context_first_can_select_weak_discordant_rna():
    report = json.loads((DATA / "sample_audit.json").read_text())
    row, = [r for r in report["rows"] if r["sample"] == "ont_t1" and r["variant"]["gene"] == "H1-2"]
    balanced, = row["primary_only"]["protein_default"]["top_proteins"]
    context, = row["primary_only"]["protein_context"]["top_proteins"]
    assert balanced["matches_expected"]
    assert balanced["fraction_of_best_candidate_support"] >= .85
    assert len(context["amino_acids"]) == 49
    assert context["supporting_read_names"] == 6
    assert not context["matches_expected"]


@pytest.mark.parametrize("sample", MANIFEST["datasets"])
@pytest.mark.parametrize("peptide", [15, 25, 30])
@pytest.mark.parametrize("floor", [2, 5, 10])
def test_real_dync1h1_adapts_to_peptide_size_and_coverage_floor(protein_cases, sample, peptide, floor):
    variant, evidence, expected = protein_cases[sample, "DYNC1H1"]
    result = audit.protein_result(variant, evidence, expected, coverage=floor,
                                 protein_context_peptide_length=peptide)
    assert result["status"] == "ok", result
    assert result["target_protein_length"] == 2 * peptide - 1
    assert result["peptide_length"] == peptide
    assert result["min_support_fraction"] == .85
    for top in result["top_proteins"]:
        assert top["matches_expected"], result
        assert top["fraction_of_best_candidate_support"] >= .85
        assert min(top["retained_cdna_min_coverages"]) >= floor
    if not result["top_proteins"]:
        # A selected fixture with too few reads must not bypass the floor.
        assert result["outcome"] in ("no_rna_candidate", "no_translation")


def test_full_region_adaptive_sweep_records_both_support_metrics():
    report = json.loads((DATA / "sample_audit.json").read_text())
    for row in report["rows"]:
        if row["variant"]["gene"] != "DYNC1H1":
            continue
        sweep = row["primary_only"]["protein_adaptive_sweep"]
        assert {(r["peptide_length"], r["min_coverage"]) for r in sweep} == {
            (p, c) for p in (15, 25, 30) for c in (2, 5, 10)}
        for result in sweep:
            assert result["outcome"] == "expected_top_protein"
            assert result["target_protein_length"] == 2 * result["peptide_length"] - 1
            top, = result["top_proteins"]
            assert min(top["retained_cdna_min_coverages"]) >= result["min_coverage"]
            assert top["fraction_of_best_candidate_support"] >= .85
        if row["sample"] == "ont_t1":
            p25 = {r["min_coverage"]: r["top_proteins"][0] for r in sweep if r["peptide_length"] == 25}
            assert len(p25[2]["amino_acids"]) == 37
            assert len(p25[5]["amino_acids"]) == 33
