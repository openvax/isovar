"""Independent CIGAR semantics and original insertion/compound RNA probes."""

from collections import Counter
import json
from pathlib import Path

import pysam
import pytest
from varcode import Variant

from tests.data.osteosarc.expansion.inventory import digest
from tests.data.osteosarc.expansion.references import apply_variant, load_reference, minimal_edit, reference_genome
from tests.data.osteosarc.expansion.runner import audit_mode, PRIMARY_EXCLUDE_FLAGS
from tests.data.osteosarc.expansion.stress import exact_allele
from tests.real_rna_helpers import record_digest
from tests.test_osteosarc_expansion_corpus import canonical_result


DIRECTORY = Path(__file__).parent / "data/osteosarc/expansion/stress-corpus"
CASES = json.loads((DIRECTORY / "manifest.json").read_text())["cases"]


@pytest.fixture(scope="session")
def stress_reference(tmp_path_factory):
    manifest, models = load_reference(DIRECTORY / "reference")
    return manifest, models, reference_genome(DIRECTORY / "reference", tmp_path_factory.mktemp("stress-reference"))


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["case_id"])
def test_stress_original_records_and_independent_cigar(case):
    for name, checksum in case["files"].items():
        assert digest(DIRECTORY / name) == checksum
    with pysam.AlignmentFile(DIRECTORY / case["bam"]) as bam:
        reads = list(bam)
    assert [record_digest(r) for r in reads] == [r["sam_sha256"] for r in case["selected_records"]]
    _, ref, alt = minimal_edit(case["variant"])
    counts = Counter()
    for read in reads:
        if read.flag & PRIMARY_EXCLUDE_FLAGS or read.mapping_quality < 1:
            continue
        allele = exact_allele(read, case["variant"])
        counts["uncallable" if allele is None else "ref" if allele == ref else "alt" if allele == alt else "other"] += 1
    assert {g: counts[g] for g in ("ref", "alt", "other", "uncallable")} == case["independent_primary"]["counts"]


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["case_id"])
def test_stress_public_protein_pipeline(case, stress_reference):
    manifest, models, genome = stress_reference
    record = case["variant"]
    variant = Variant(record["chrom"].removeprefix("chr"), record["pos"], record["ref"], record["alt"], ensembl=genome)
    expected = {tid: apply_variant(record, models[tid]) for tid in manifest["variant_transcripts"][record["variant_id"]]}
    actual = audit_mode(DIRECTORY / case["bam"], variant, expected)
    assert canonical_result(actual) == canonical_result(case["defaults"])
    assert actual["status"] == "ok"
    assert all(p["validation_status"] == "ok" for p in actual["proteins"])


def test_compound_definition_explains_ntf3_rna_without_inventing_vaccine_membership():
    positives = [c for c in CASES if c["defaults"]["proteins"]]
    assert {c["source_id"] for c in positives} == {"7010dd389c50467d", "1b66c15da594a3ef"}
    for case in positives:
        assert case["variant"]["ref"] == "AG" and case["variant"]["alt"] == "GT"
        assert case["variant"]["vaccine_count"] is None
        assert case["defaults"]["outcome"] == "expected_top_protein"
    assert not any(c["independent_primary"]["counts"]["alt"] for c in CASES if c["variant"]["gene"] != "NTF3")


def test_pacbio_missing_qualities_explain_lost_compound_evidence():
    case = next(c for c in CASES if c["source_id"] == "0066232879babe83" and c["variant"]["gene"] == "NTF3")
    assert case["independent_primary"]["counts"] == dict(ref=1, alt=1, other=0, uncallable=0)
    with pysam.AlignmentFile(DIRECTORY / case["bam"]) as bam:
        assert all(r.query_qualities is None for r in bam)
    assert case["defaults"]["outcome"] == "no_callable_allele"
    assert case["defaults"]["counts"]["reads"] == dict(ref=0, alt=0, other=0)


@pytest.mark.parametrize("cigar,sequence,pos,ref,alt,expected", [
    ("6M", "AAGTCC", 103, "AG", "GT", "GT"),
    ("2M2I3M", "CAGGCCC", 102, "A", "AGG", "GG"),
    ("2M2D3M", "CACCC", 102, "AGG", "A", ""),
    ("2M2I2D3M", "CAGCCCC", 103, "AT", "GC", "GC"),
    ("2M2N3M", "CACCC", 102, "AGG", "A", None),
    ("2S5M", "TTCCGTA", 103, "AG", "GT", "GT"),
    ("2H5M", "CCGTA", 103, "AG", "GT", "GT"),
    ("2=2X2=", "AAGTCC", 103, "AG", "GT", "GT"),
    ("3M", "AGT", 100, "A", "G", None),
    ("3M", "AGT", 103, "T", "C", None),
])
def test_independent_exact_allele_cigar_anchors(cigar, sequence, pos, ref, alt, expected):
    header = pysam.AlignmentHeader.from_dict(dict(SQ=[dict(SN="chr1", LN=1000)]))
    read = pysam.AlignedSegment(header)
    read.reference_id, read.reference_start, read.cigarstring, read.query_sequence = 0, 100, cigar, sequence
    assert exact_allele(read, dict(chrom="chr1", pos=pos, ref=ref, alt=alt)) == expected
