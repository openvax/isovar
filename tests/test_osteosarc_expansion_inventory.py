"""The expanded audit must not silently change its cohort or source identities."""

import json

import pytest

from .data.osteosarc.expansion.inventory import (
    alignment_scope, catalogue_sources, digest, fetch_snapshot,
    reconcile_bucket, vaccine_variants,
)


def index_row(identity="GENE-chr1-10", count=1, location="chr1:10"):
    return (f'<tr data-vaccines="{count}"><td><a href="/variant/{identity}">GENE</a></td>'
            f'<td>{location}</td><td>p.Ala1Val</td>' + "<td></td>" * 6 + "</tr>")


def vaf(alleles=("GENE-chr1-10\tchr1\t10\tC\tT",)):
    return "variant_id\tchrom\tpos\tref\talt\n" + "\n".join(alleles) + "\n"


def test_exact_variant_join_not_gene_join():
    html = index_row() + index_row("GENE-chr1-11", location="chr1:11") + index_row("OUT", count=0)
    variants = vaccine_variants(html, vaf(("GENE-chr1-10\tchr1\t10\tC\tT", "GENE-chr1-11\tchr1\t11\tG\tA")))
    assert [(v["pos"], v["ref"], v["alt"]) for v in variants] == [(10, "C", "T"), (11, "G", "A")]


@pytest.mark.parametrize("html,table", [
    (index_row(), vaf(())),
    (index_row() * 2, vaf()),
    (index_row(location="chr1:11"), vaf()),
    (index_row(), vaf(("GENE-chr1-10\tchr1\t10\tC\tT", "GENE-chr1-10\tchr1\t10\tC\tA"))),
    (index_row(), vaf(("GENE-chr1-10\tchr1\t10\tC\t<DEL>",))),
])
def test_invalid_cohort_fails_before_querying_reads(html, table):
    with pytest.raises(ValueError):
        vaccine_variants(html, table)


@pytest.mark.parametrize("key,expected", [
    ("rna-seq/tempus/x/WES/tumor.bam", "excluded_dna_alignment"),
    ("rna-seq/tempus/x/RNA/tumor.bam", "rna_candidate"),
    ("rna-seq/reprocessed/x.Aligned.toTranscriptome.out.bam", "excluded_transcript_coordinate_alignment"),
    ("ucsf/T2/x/vdj_t/consensus.bam", "excluded_receptor_or_targeted_alignment"),
    ("hudson_lab/PBMC/capture/count/sample_alignments.bam", "rna_candidate"),
    ("kamil/blood/output/TCRab/outs/all_contig.bam", "excluded_receptor_or_targeted_alignment"),
    ("vendor/cegat/P116686_2_S000048/P116686_3.bam", "rna_candidate"),
    ("vendor/cegat/P116686_2_S000048/P116686_2.bam", "excluded_dna_alignment"),
    ("ucsf/unknown.bam", "rna_candidate"),  # Must inspect its header, not assume genomic.
    ("unknown.bam", "unresolved_assay"),
    ("rna-seq/x.bam.bai", None),
])
def test_source_classification(key, expected):
    assert alignment_scope(key) == expected


def test_source_claims_preserve_conflicts_and_pool_attribution():
    source = {"name": "T0", "url": "rna-seq/BG009368.bam", "tissue": "tumor"}
    table = ("bam_file\tsample_label\ttimepoint\tsample_date\tassay_type\n"
             "BG009368.bam\tT1\tT1\t2024-06\tRNA\n")
    result = catalogue_sources(dict(baseUrl="https://example/", genome="hg38", categories=[
        dict(name="RNA", bams=[source]),
        dict(name="Blood scRNA", bams=[dict(source, url="blood/pooled.bam")]),
    ]), table)
    assert result[0]["name"] == "T0"
    assert result[0]["count_export_claims"][0]["timepoint"] == "T1"
    assert result[1]["attribution"] == "pooled_library_unresolved"


def test_bucket_union_preserves_catalogue_missing_products_and_unassigned():
    key = "rna-seq/catalogue-only.bam"
    value = dict(sources={}, variants=[], catalogue_alignments=[dict(key=key, name="catalogue")])
    listing = dict(generated_at="date", files=[
        ["kamil/tumor/outs/unassigned_alignments.bam", 10, 1],
        ["kamil/tumor/outs/unassigned_alignments.bam.bai", 2, 1],
    ])
    result = reconcile_bucket(value, listing)["alignments"]
    assert len(result) == 2
    missing = next(r for r in result if r["key"] == key)
    assert not missing["listed_in_bucket"]
    assert missing["bucket_bytes"] is None
    unassigned = next(r for r in result if r["key"] != key)
    assert unassigned["attribution"] == "unassigned_library"
    assert len(unassigned["listed_indexes"]) == 1


def test_cached_snapshot_is_checked_without_network(tmp_path, monkeypatch):
    path = tmp_path / "metadata"
    path.write_bytes(b"original source")
    receipt = {"url": "https://example/", "sha256": digest(path), "bytes": path.stat().st_size}
    path.with_name("metadata.receipt.json").write_text(json.dumps(receipt))
    monkeypatch.setattr("subprocess.run", lambda *a, **kw: pytest.fail("Unexpected network"))
    assert fetch_snapshot("https://example/", path) == receipt
    with pytest.raises(ValueError, match="drift"):
        fetch_snapshot("https://different/", path)
    path.write_bytes(b"tampered")
    with pytest.raises(ValueError, match="drift"):
        fetch_snapshot("https://example/", path)
