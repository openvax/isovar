"""Normal source selection uses declared IDs and sample columns, not file naming."""

from dataclasses import asdict, replace
import json

import pysam
import pytest
from varcode import Variant

from isovar.normal_sources import NormalSource, assess_matched_normal
from isovar.sv_rna_audit import file_sha256


def source_file(tmp_path, *, contig="chr1", start=104, alleles=("A", "G", "T"), genotype=(0, 2)):
    path = tmp_path / "ArbitraryCase_and_not_a_Normal_suffix.vcf"
    header = pysam.VariantHeader()
    header.contigs.add(contig, length=1000)
    header.add_sample("EXPLICIT.sample-column")
    for name, kind in (("GT", "String"), ("GQ", "Integer"), ("DP", "Integer"), ("MIN_DP", "Integer")):
        header.formats.add(name, 1, kind, name)
    with pysam.VariantFile(path, "w", header=header) as writer:
        row = writer.new_record(contig=contig, start=start, stop=start + len(alleles[0]), alleles=alleles)
        row.samples["EXPLICIT.sample-column"].update(GT=genotype, GQ=40, DP=12)
        row.filter.add("PASS")
        writer.write(row)
    compressed = str(path) + ".gz"
    pysam.tabix_compress(str(path), compressed)
    pysam.tabix_index(compressed, preset="vcf")
    return NormalSource(source_id="public-source-1", reference_name="GRCh38", evidence_kind="variant_only",
        selected=True, selection_reason="Reviewed sample identity", specimen_id="blood-draw-1",
        library_id="library-1", vcf_sample="EXPLICIT.sample-column", caller="test-caller",
        processing_group="library-1-calling", path=compressed, sha256=file_sha256(compressed),
        index_path=compressed + ".tbi", index_sha256=file_sha256(compressed + ".tbi"))


def assess(variants, sources, **kwargs):
    return assess_matched_normal(variants, sources, sample_id="patient", reference_name="GRCh38", **kwargs)


def test_explicit_source_selection_and_genotype_indexed_background(tmp_path):
    source = source_file(tmp_path)
    result = assess([Variant("1", 105, "A", "G")], [source])
    assert result["alleles"][0]["state"] == "other_allele_observed"
    background, = result["background_variants"]
    assert background["alt"] == "T"
    assert background["origins"][0]["record"]["genotype"] == [0, 2]
    assert result["source_inventory"][0]["library_id"] == "library-1"


def test_missing_contig_is_unavailable_not_absent_normal_alt(tmp_path):
    source = source_file(tmp_path)
    result = assess([Variant("2", 105, "A", "G")], [source])
    row = result["alleles"][0]
    assert row["state"] == "unavailable"
    assert row["source_assessments"][0]["reason"] == "unresolved_chromosome"
    assert not row["origin_conflict"]


@pytest.mark.parametrize("contig,query", [("chrM", "MT"), ("HPV16", "HPV16"),
                                         ("chrUn_gl000220", "chrUn_gl000220")])
def test_mitochondrial_viral_and_decoy_names_use_shared_resolver(tmp_path, contig, query):
    source = source_file(tmp_path, contig=contig, genotype=(0, 1))
    result = assess([Variant(query, 105, "A", "G")], [source])
    assert result["alleles"][0]["state"] == "alternate_observed"


def test_subset_and_excluded_sources_are_not_complete_normal_scope(tmp_path):
    source = source_file(tmp_path)
    outside = replace(source, acquired_regions=(("chr1", 200, 300),))
    excluded = replace(source, source_id="another-product", selected=False, path=None,
                       selection_reason="Unresolved specimen assignment")
    result = assess([Variant("1", 105, "A", "G")], [outside, excluded])
    assert result["alleles"][0]["state"] == "unavailable"
    assert [r["reason"] for r in result["alleles"][0]["source_assessments"]] == ["outside_acquired_regions"]
    assert result["source_inventory"][1]["selection_reason"] == "Unresolved specimen assignment"


def test_partial_background_scope_remains_explicit(tmp_path):
    source = replace(source_file(tmp_path), acquired_regions=(("chr1", 104, 105),))
    result = assess([Variant("1", 105, "A", "G")], [source], background_flank=10)
    row = result["alleles"][0]["source_assessments"][0]
    assert row["state"] == "other_allele_observed"
    assert row["background_query_interval"] == [94, 115]
    assert row["background_interval_fully_acquired"] is False


def test_pinned_identity_and_reference_context_are_checked(tmp_path):
    source = source_file(tmp_path)
    with pytest.raises(ValueError, match="digest mismatch"):
        assess([Variant("1", 105, "A", "G")], [replace(source, sha256="0" * 64)])
    with pytest.raises(ValueError, match="index digest mismatch"):
        assess([Variant("1", 105, "A", "G")], [replace(source, index_sha256="0" * 64)])
    with pytest.raises(ValueError, match="incompatible reference"):
        assess([Variant("1", 105, "A", "G", genome="GRCh37")], [source])
    with pytest.raises(ValueError, match="Duplicate normal source"):
        assess([Variant("1", 105, "A", "G")], [source, source])


def test_cli_and_api_share_evidence_and_preserve_input_receipts(tmp_path, capsys):
    from isovar.cli.commands import run
    source = source_file(tmp_path)
    nominations = tmp_path / "nominations.json"
    nominations.write_text(json.dumps([dict(contig="1", start=105, ref="A", alt="G", origins=["tumor-call"])]))
    manifest = tmp_path / "normal.json"
    manifest.write_text(json.dumps(dict(schema="isovar.normal_manifest.v1", sample_id="patient", reference_name="GRCh38",
        snapshot="public-snapshot", variants=dict(path=nominations.name, sha256=file_sha256(nominations)),
        sources=[dict(asdict(source), path=source.path.split('/')[-1])], background_flank=0)))
    output = tmp_path / "normal-result.json"
    run(["normal-evidence", "--manifest", str(manifest), "--output", str(output)])
    cli = json.loads(output.read_text())
    api = json.loads(json.dumps(assess([Variant("1", 105, "A", "G")], [source])))
    assert cli["alleles"][0].pop("nomination")["origins"] == ["tumor-call"]
    assert cli["alleles"] == api["alleles"]
    assert cli["background_variants"] == api["background_variants"]
    prior = output.read_bytes()
    with pytest.raises(SystemExit):
        run(["normal-evidence", "--manifest", str(manifest), "--output", str(output)])
    assert output.read_bytes() == prior
    assert "already exists" in capsys.readouterr().err
