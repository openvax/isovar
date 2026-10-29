"""Normal scope follows a cohort into its results, provenance, reports and resume identity."""

import gzip
from pathlib import Path

import pytest
from pypdf import PdfReader
from varcode import Variant

from isovar import assess_matched_normal
from isovar.cohort import run_cohort
from isovar.sv_rna_audit import file_sha256, read_json, write_json
from .test_normal_sources import source_file
from .testing_helpers import data_path


def test_normal_background_and_scope_are_preserved_in_cohort(tmp_path):
    normal_source = source_file(tmp_path, contig="chr14", start=92811640,
                                alleles=("GAGA", "G"), genotype=(0, 1))
    normal = assess_matched_normal([Variant("14", 92811641, "GAGA", "G")], [normal_source],
                                  sample_id="COLO829", reference_name="GRCh38")
    normal_path, calls = tmp_path / "normal.json", tmp_path / "calls.json"
    write_json(normal_path, normal)
    write_json(calls, [dict(contig="14", start=92811641, ref="GAGA", alt="G"),
                       dict(contig="not-a-contig", start=100, ref="A", alt="G")])
    bam = Path(data_path("data/colo829_whole_run/SRR8615617.all_sites.grch38.bam"))
    manifest = dict(schema="isovar.cohort_manifest.v1", samples=[dict(sample_id="COLO829")],
        files={name: dict(path=str(path), sha256=file_sha256(path)) for name, path in
               [("bam", bam), ("calls", calls), ("normal", normal_path)]},
        sources=[dict(source_id="original", sample_id="COLO829", library_id="colo-library", bam="bam",
                      scope="Regional RNA", provenance={"accession": "SRR8615617"})],
        small_variants=[dict(source_id="original", variants="calls", normal_evidence="normal", ensembl_release=87,
                             protein_settings=dict(protein_sequence_length=33))])
    manifest_path, output = tmp_path / "manifest.json", tmp_path / "out"
    write_json(manifest_path, manifest)
    reports = run_cohort(manifest_path, output)
    receipt = read_json(output / "sources/original/default/receipt.json")
    scope = receipt["matched_normal_scope"]
    assert scope["source_inventory"][0]["vcf_sample"] == "EXPLICIT.sample-column"
    assert scope["allele_states"] == {"alternate_observed": 1, "unavailable": 1}
    assert scope["receipt_sha256"] == file_sha256(normal_path)
    with gzip.open(output / "sources/original/default/variants.jsonl.gz", "rt") as stream:
        import json
        good, missing = [json.loads(line) for line in stream]
    assert good["nomination"]["origin_conflict"]
    assert good["nomination"]["matched_normal_evidence"]["state"] == "alternate_observed"
    assert missing["nomination"]["matched_normal_evidence"]["reason"] == "allele_not_assessed_in_normal_receipt"
    text = "\n".join(page.extract_text() for page in PdfReader(reports[0]).pages)
    assert "Matched-normal evidence scope" in text
    assert "public-source-1" in text and "library-1" in text
    # A changed policy/scope is a new scientific result, even if the RNA is unchanged.
    normal["source_inventory"][0]["policy"]["minimum_depth"] = 7
    write_json(normal_path, normal)
    manifest["files"]["normal"]["sha256"] = file_sha256(normal_path)
    write_json(manifest_path, manifest)
    with pytest.raises(ValueError, match="Resume inputs"):
        run_cohort(manifest_path, output, resume=True, render=False)
    with pytest.raises(ValueError, match="analysis identity differs"):
        run_cohort(manifest_path, output, render_only=True, render=False)


def test_cohort_rejects_unreconciled_background_inputs(tmp_path):
    from isovar.cohort import _normal_background
    with pytest.raises(ValueError, match="reconcile the normal inventory first"):
        _normal_background(dict(normal_evidence="normal", germline_variants="old"), {}, {}, [], [], None)
