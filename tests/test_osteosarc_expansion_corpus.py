"""Composed, network-free regressions for all 44 vaccine loci and native MT references."""

from collections import Counter
from copy import deepcopy
import json
from pathlib import Path
from unittest.mock import patch

import pysam
import pytest
from varcode import Variant

from tests.data.osteosarc.expansion.inventory import digest
from tests.data.osteosarc.expansion.references import apply_variant, load_reference, reference_genome, translate
from tests.data.osteosarc.expansion.runner import audit_mode, independent_counts, PRIMARY_EXCLUDE_FLAGS
from tests.real_rna_helpers import record_digest
from tests.test_shared_support_index import exhaustive_support


CORPUS = Path(__file__).parent / "data/osteosarc/expansion/corpus"
MANIFEST = json.loads((CORPUS / "manifest.json").read_text())
CASES = MANIFEST["cases"]


@pytest.fixture(scope="session")
def cohort_references(tmp_path_factory):
    cache = tmp_path_factory.mktemp("osteosarc-cohort-references")
    result = {}
    for name in sorted({case["reference"] for case in CASES}):
        directory = CORPUS / "references" / name
        manifest, models = load_reference(directory)
        result[name] = (manifest, models, reference_genome(directory, cache / name))
    return result


def canonical_result(value):
    result = deepcopy(value)
    for protein in result.get("proteins", []):
        # Translation contributions within one ranked protein are an
        # unordered group; the order of the ranked proteins remains tested.
        protein["checks"].sort(key=lambda c: json.dumps(c, sort_keys=True))
    return result


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["case_id"])
def test_original_records_and_primary_derivative_are_unchanged(case):
    for name, checksum in case["files"].items():
        assert digest(CORPUS / name) == checksum
    with pysam.AlignmentFile(CORPUS / case["bam"]) as bam:
        reads = list(bam)
    assert [record_digest(r) for r in reads] == [r["sam_sha256"] for r in case["selected_records"]]
    with pysam.AlignmentFile(CORPUS / case["primary_bam"]) as bam:
        primary = list(bam)
    assert [record_digest(r) for r in primary] == [record_digest(r) for r in reads if not r.flag & PRIMARY_EXCLUDE_FLAGS]


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["case_id"])
def test_independent_exact_cigar_and_quality_counts(case):
    record = case["variant"]
    with pysam.AlignmentFile(CORPUS / case["bam"]) as bam:
        canonical = record["chrom"].removeprefix("chr")
        aliases = {"M", "MT", "chrM", "chrMT"} if canonical in ("M", "MT") else {canonical, "chr" + canonical}
        names = aliases & set(bam.references)
        assert len(names) == 1
        result = independent_counts(bam, dict(record, chrom=next(iter(names))))
    assert result == case["independent_primary"]


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["case_id"])
@pytest.mark.parametrize("mode", ["defaults", "primary_only"])
def test_public_read_to_ranked_protein_pipeline(case, mode, cohort_references):
    manifest, models, genome = cohort_references[case["reference"]]
    record = case["variant"]
    contig = "MT" if record["chrom"] == "chrM" else record["chrom"].removeprefix("chr")
    variant = Variant(contig, record["pos"], record["ref"], record["alt"], ensembl=genome)
    expected = {tid: apply_variant(record, models[tid]) for tid in manifest["variant_transcripts"][record["variant_id"]]}
    filename = case["bam"] if mode == "defaults" else case["primary_bam"]
    actual = audit_mode(CORPUS / filename, variant, expected)
    assert actual["status"] == "ok"
    assert canonical_result(actual) == canonical_result(case[mode])
    for protein in actual["proteins"]:
        assert protein["validation_status"] == "ok"
        assert all(c["status"] == "ok" for c in protein["checks"])


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["case_id"])
def test_uncapped_ranking_capture_preserves_exact_public_default(case, cohort_references):
    manifest, models, genome = cohort_references[case["reference"]]
    record = case["variant"]
    contig = "MT" if record["chrom"] == "chrM" else record["chrom"].removeprefix("chr")
    variant = Variant(contig, record["pos"], record["ref"], record["alt"], ensembl=genome)
    expected = {tid: apply_variant(record, models[tid]) for tid in manifest["variant_transcripts"][record["variant_id"]]}
    actual = audit_mode(CORPUS / case["bam"], variant, expected, capture_all_ranked=True)
    uncapped = actual.pop("uncapped_ranked_proteins")
    assert actual.pop("returned_protein_limit") == 1
    assert actual.pop("uncapped_validation_status") == "ok"
    assert actual.pop("uncapped_all_match_expected") == (all(p["matches_expected"] for p in uncapped) if uncapped else None)
    assert canonical_result(actual) == canonical_result(case["defaults"])
    assert actual["proteins"] == uncapped[:1]
    if record["variant_id"] == "DYNC1H1-chr14-102030200":
        assert len(uncapped) > 1  # The test must really exercise hidden alternatives.


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["case_id"])
@pytest.mark.parametrize("mode", ["defaults", "primary_only"])
def test_all_original_read_rankings_match_exhaustive_support(case, mode, cohort_references):
    """Compare EVERY ranked window, not just the one returned by default."""
    manifest, models, genome = cohort_references[case["reference"]]
    record = case["variant"]
    contig = "MT" if record["chrom"] == "chrM" else record["chrom"].removeprefix("chr")
    variant = Variant(contig, record["pos"], record["ref"], record["alt"], ensembl=genome)
    expected = {tid: apply_variant(record, models[tid]) for tid in manifest["variant_transcripts"][record["variant_id"]]}
    path = CORPUS / case["bam" if mode == "defaults" else "primary_bam"]
    actual = audit_mode(path, variant, expected, capture_all_ranked=True)
    with patch("isovar.variant_sequence_helpers._variant_sequences_with_shared_read_support", exhaustive_support):
        baseline = audit_mode(path, variant, expected, capture_all_ranked=True)
    for result in (actual, baseline):
        assert result["status"] == "ok"
        for protein in result["uncapped_ranked_proteins"]:
            protein["checks"].sort(key=lambda c: json.dumps(c, sort_keys=True))
    assert canonical_result(actual) == canonical_result(baseline)


def test_every_original_reference_model_independently_translates(cohort_references):
    for manifest, models, _ in cohort_references.values():
        assert len(models) == manifest["transcript_count"]
        for model in models.values():
            assert translate(model["cdna"][model["cds_start"]:], model["genetic_code"], True) == (
                model["protein"], model["reference_has_stop"])


def test_corpus_scope_contains_every_vaccine_locus_and_explicit_negative_examples():
    cohort = CASES[:44]
    assert len({c["variant"]["variant_id"] for c in cohort}) == 44
    assert Counter("SNV" if len(c["variant"]["ref"]) == len(c["variant"]["alt"]) == 1 else "deletion"
                   for c in cohort) == {"SNV": 37, "deletion": 7}
    outcomes = Counter(c["defaults"]["outcome"] for c in cohort)
    assert outcomes["no_alt_reads"] == 2
    assert outcomes["no_rna_candidate"] == 4
    assert outcomes["different_top_protein"] == 1
    assert outcomes["expected_top_protein"] == 37


def test_real_mitochondrial_window_requires_table_two(cohort_references):
    case = next(c for c in CASES if c["variant"]["chrom"] == "chrM" and c["source_id"] == "f30f618fb76a0e49")
    assert case["defaults"]["counts"]["reads"] == {"ref": 10, "alt": 27, "other": 1}
    protein = case["defaults"]["proteins"][0]
    assert protein["amino_acids"] == "WDPQQMALLNANPSLTPLLGLLLATAGKSAQLGLHPWLPSAMEGPTPVS"
    assert protein["mutation_interval"] == [24, 25]
    assert protein["length"] == 49 and protein["mutation_containing_peptide_windows"] == 25
    manifest, models, _ = cohort_references["GRCh38"]
    tid, = manifest["variant_transcripts"][case["variant"]["variant_id"]]
    expected = apply_variant(case["variant"], models[tid])
    start = 3 * (protein["checks"][0]["protein_start_1based"] - 1) + expected["cds_start"]
    cdna = expected["mutant_cdna"][start:start + 147]
    assert translate(cdna, 2) == (protein["amino_acids"], False)
    assert translate(cdna, 1) == ("", True)
    assert {cdna[i:i + 3] for i in range(0, len(cdna), 3)} >= {"TGA", "ATA"}


def test_cegat_mitochondrial_reference_shift_is_explicit(cohort_references):
    manifest, models, _ = cohort_references["GRCh37-cegat"]
    case = next(c for c in CASES if c["reference"] == "GRCh37-cegat")
    assert case["variant"]["pos"] == 12995
    assert case["variant"]["original_identity"]["pos"] == 12994
    tid = "ENST00000361567"
    derivation = manifest["mitochondrial_annotation_derivation"]["sequence_validations"][tid]
    assert len(derivation["reference_substitutions"]) == 2
    assert models[tid]["cdna"] == cohort_references["GRCh37"][1][tid]["cdna"]
    assert case["defaults"]["outcome"] == "expected_top_protein"


def test_native_grch37_recovers_nr2f2_evidence_absent_from_grch38_corpus():
    native = next(c for c in CASES if c["variant"]["gene"] == "NR2F2" and c["reference"] == "GRCh37-cegat")
    grch38_case = next(c for c in CASES if c["variant"]["gene"] == "NR2F2" and c["reference"] == "GRCh38")
    assert native["variant"]["pos"] == 96875528
    assert native["defaults"]["counts"]["reads"]["alt"] == 7
    assert native["defaults"]["outcome"] == "different_top_protein"
    assert native["defaults"]["proteins"][0]["validation_status"] == "ok"
    assert grch38_case["defaults"]["outcome"] == "no_alt_reads"
