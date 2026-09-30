"""Offline RNA evidence and protein reconstruction for the full Sid union panel."""

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import socket

from osteosarc import bundle_file
import pysam
import pytest
from varcode import Variant

from isovar import sid_data
from isovar.allele_read import AlleleRead
from isovar.read_collector import ReadCollector
from tests.data.osteosarc.expansion.references import (
    apply_variant, load_reference, minimal_edit, reference_genome, translate,
)
from tests.data.osteosarc.expansion.runner import audit_mode
from tests.osteosarc_union_helpers import canonical_audit, direct_evidence, target_record
from tests.real_rna_helpers import record_digest


DATA = Path(__file__).parent / "data/osteosarc/union"
MANIFEST = json.loads((DATA / "manifest.json").read_text())
CASES = [(n, c) for n, c in MANIFEST["cases"].items() if c["status"] == "tested"]


@pytest.fixture(scope="module", autouse=True)
def offline():
    def forbidden(*args, **kwargs):
        raise AssertionError("Sid union regressions must run without network access")
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(socket.socket, "connect", forbidden)
        patch.setattr(socket, "create_connection", forbidden)
        yield


@pytest.fixture(scope="module")
def panel(tmp_path_factory, offline):
    cache = tmp_path_factory.mktemp("sid-union")
    paths = {n: bundle_file(sid_data.bundle(), c["member"], format="bam", cache=cache) for n, c in CASES}
    reference, models = load_reference(DATA / "reference")
    genome = reference_genome(DATA / "reference", cache / "reference")
    return paths, reference, models, genome


def test_every_union_target_is_accounted_for():
    recipe = json.loads((sid_data.bundle() / "recipe.json").read_text())
    targets = {n: t for n, t in recipe["targets"].items()
               if t["kind"] == "small_variant" or (t["kind"] == "unresolved" and t.get("label") == "current")}
    assert targets == MANIFEST["targets"]
    assert set(MANIFEST["cases"]) == set(targets)
    assert Counter(t["kind"] for t in targets.values()) == {"small_variant": 187, "unresolved": 2}
    assert Counter(c["status"] for c in MANIFEST["cases"].values()) == {
        "tested": 184, "no_T2_member_for_reference": 3, "unresolved_allele": 2}
    members = {n: v for n, v in recipe["members"].items()
               if v["source"] == MANIFEST["source"] and v["target"] in targets}
    assert {c["member"] for _, c in CASES} == set(members)
    assert len(CASES) == len(members)
    for name, case in MANIFEST["cases"].items():
        if case["status"] == "tested":
            assert members[case["member"]]["target"] == name
            assert targets[name]["assembly"] == "GRCh38"
        elif case["status"] == "no_T2_member_for_reference":
            assert targets[name]["assembly"] == "GRCh37"
            assert not any(v["target"] == name for v in recipe["members"].values())
        else:
            assert targets[name]["kind"] == "unresolved" and targets[name]["reason"]
    assert MANIFEST["shared_manifest_sha256"] == sid_data.BUNDLE_MANIFEST_SHA256
    assert sha256((DATA / "reference/manifest.json").read_bytes()).hexdigest() == MANIFEST["reference_manifest_sha256"]


@pytest.mark.parametrize("name,case", CASES, ids=[n for n, _ in CASES])
def test_direct_alleles_match_independent_cigar_evidence(name, case, panel):
    paths, _, _, _ = panel
    record = target_record(name, MANIFEST["targets"][name])
    with pysam.AlignmentFile(paths[name]) as bam:
        reads = list(bam)
    assert len(reads) == case["records"]
    expected = direct_evidence(reads, record)
    assert expected == case["direct"]
    observations = {o["record"]: o for o in expected["observations"]}
    pos, ref, _ = minimal_edit(record)
    start, end = pos - 1, pos - 1 + len(ref)
    collector = ReadCollector(use_secondary_alignments=False, use_soft_clipped_bases=True,
                              merge_overlapping_fragments=False)
    for read in reads:
        observation = observations.get(record_digest(read))
        if observation is None:
            continue
        locus = collector.locus_read_from_pysam_aligned_segment(read, start, end)
        assert locus is not None, read.query_name
        allele = AlleleRead.from_locus_read(locus)
        assert allele is not None, read.query_name
        assert allele.allele == observation["allele"], read.query_name
        assert [locus.read_base0_start_inclusive, locus.read_base0_end_exclusive] == observation["query_interval"]


@pytest.mark.parametrize("name,case", [(n, c) for n, c in CASES
                                     if len(MANIFEST["targets"][n]["ref"]) == len(MANIFEST["targets"][n]["alt"]) == 1])
def test_snv_evidence_partition_matches_independent_records(name, case, panel):
    paths, _, _, genome = panel
    record = target_record(name, MANIFEST["targets"][name])
    variant = Variant(record["chrom"].removeprefix("chr"), record["pos"], record["ref"], record["alt"], ensembl=genome)
    with pysam.AlignmentFile(paths[name]) as bam:
        evidence = ReadCollector(use_secondary_alignments=False, merge_overlapping_fragments=False,
                                 read_filter=lambda read: not read.is_supplementary).read_evidence_for_variant(variant, bam)
    for group in ("ref", "alt", "other"):
        observed = Counter((r.name, r.allele) for r in getattr(evidence, group + "_reads"))
        expected = Counter((o["name"], o["allele"]) for o in case["direct"]["observations"]
                           if ("ref" if o["allele"] == record["ref"] else "alt" if o["allele"] == record["alt"] else "other") == group)
        assert observed == expected


@pytest.mark.parametrize("name,case", CASES, ids=[n for n, _ in CASES])
def test_public_rna_to_protein_pipeline(name, case, panel):
    paths, reference, models, genome = panel
    record = target_record(name, MANIFEST["targets"][name])
    variant = Variant(record["chrom"].removeprefix("chr"), record["pos"], record["ref"], record["alt"], ensembl=genome)
    expected = {tid: apply_variant(record, models[tid]) for tid in reference["variant_transcripts"][name]}
    actual = canonical_audit(audit_mode(paths[name], variant, expected, capture_all_ranked=True))
    assert actual["status"] == "ok", actual
    uncapped = actual.pop("uncapped_ranked_proteins")
    assert actual.pop("uncapped_validation_status") == "ok"
    assert actual.pop("uncapped_all_match_expected") == (all(p["matches_expected"] for p in uncapped) if uncapped else None)
    assert actual.pop("returned_protein_limit") == 1
    assert actual["proteins"] == uncapped[:1]
    assert actual == case["pipeline"]
    for protein in uncapped:
        assert protein["validation_status"] == "ok", protein
        assert protein["checks"]
        assert all(c["status"] == "ok" for c in protein["checks"])
    if not expected:
        assert not actual["proteins"]


def test_original_reference_models_translate_and_cover_all_runnable_targets(panel):
    _, reference, models, _ = panel
    assert set(reference["variant_transcripts"]) == {n for n, _ in CASES}
    assert sum(bool(tids) for tids in reference["variant_transcripts"].values()) == 171
    for model in models.values():
        assert translate(model["cdna"][model["cds_start"]:], model["genetic_code"], True) == (
            model["protein"], model["reference_has_stop"])
