"""Competing interpretations of one locus are compared with exact RNA read sequences (#306)."""

import json

import pysam
import pytest
from varcode import Variant

from isovar import ReadCollector
from isovar.allele_interpretations import reconcile_allele_interpretations
from tests.test_sv_rna import HEADER, record, write_bam

# 0-based 1000-1039 on contig 1; offsets 20-21 (1-based 1021-1022) are the "AG" under study.
REFERENCE = "CATCGGATCCTAGGCTTACA" + "AG" + "TCCGATTGCAGTCAGGAT"
assert len(REFERENCE) == 40


def variant(position, ref, alt, contig="1"):
    return Variant(contig, position, ref, alt, normalize_contig_names=False)


def read(name, allele, flag=0, start=1000, group="a"):
    return record(name, "1", start, "40M", REFERENCE[:20] + allele + REFERENCE[22:], flag=flag, group=group)


INTERPRETATIONS = {
    "catalog A>G": [variant(1021, "A", "G")],
    "caller MNV": [variant(1021, "AG", "GT")],
    "caller SNV pair": [variant(1021, "A", "G"), variant(1022, "G", "T")],
    "G>T alone": [variant(1022, "G", "T")],
}


@pytest.fixture
def bam(tmp_path):
    reads = [read("c%d" % i, "GT") for i in range(3)] + [read("s0", "GG")]
    reads += [read("r%d" % i, "AG") for i in range(2)]
    # Overlapping mates of one fragment, both carrying the compound allele.
    reads += [read("pair", "GT", flag=99), read("pair", "GT", flag=147)]
    # A read ending inside the window cannot say which allele it carries.
    reads.append(record("partial", "1", 1000, "21M", REFERENCE[:21]))
    # One segment placed twice, reading different alleles: set aside as conflicting.
    reads += [read("conflict", "GT"), read("conflict", "GT", flag=256, start=999)]
    path = write_bam(tmp_path / "rna.bam", reads)
    with pysam.AlignmentFile(str(path)) as alignment_file:
        yield alignment_file


def reconcile(bam, interpretations=INTERPRETATIONS, **kwargs):
    return reconcile_allele_interpretations(interpretations, bam, sample_id="s", source="rna.bam", **kwargs)


def by_id(result):
    return {c["candidate_id"]: c for c in result["candidates"]}


def test_equivalent_representations_share_support_and_others_are_contradicted(bam):
    result = reconcile(bam)
    assert result["status"] == "assessed" and result["locus"]["window"] == [1020, 1022]
    assert result["locus"]["padding"] == [0, 0] and result["locus"]["reference_sequence"] == "AG"
    candidates = by_id(result)
    mnv, pair = candidates["caller MNV"], candidates["caller SNV pair"]
    assert mnv["rna_status"] == pair["rna_status"] == "supported"
    assert mnv["indistinguishable_from"] == ["caller SNV pair"] and pair["indistinguishable_from"] == ["caller MNV"]
    # Three single reads plus one merged mate pair; the pair is one fragment of two segments.
    assert (mnv["rna_support"]["fragments"], mnv["rna_support"]["segments"]) == (4, 5)
    assert mnv["rna_support"]["evidence_set_id"] == pair["rna_support"]["evidence_set_id"]
    # One fragment is below the two needed to call support, so not contradicted either.
    assert candidates["catalog A>G"]["rna_status"] == "insufficient_RNA"
    assert candidates["catalog A>G"]["rna_support"]["fragments"] == 1
    assert candidates["G>T alone"]["rna_status"] == "contradicted_by_informative_evidence"
    assert result["reference_allele"]["rna_support"]["fragments"] == 2
    assert result["informative"]["fragments"] == 7
    assert result["set_aside"] == dict(segments_not_spanning_window=1, conflicting_segments=1)
    assert json.loads(json.dumps(result)) == result


def test_observed_alleles_are_described_even_without_a_matching_candidate(bam):
    result = reconcile(bam, {"catalog A>G": [variant(1021, "A", "G")], "G>T alone": [variant(1022, "G", "T")]})
    alleles = {a["haplotype"]: a for a in result["observed_alleles"]}
    assert alleles["GT"]["matches"] == [] and alleles["GT"]["difference_from_reference"] == dict(
        start=1021, ref="AG", alt="GT")
    assert alleles["AG"]["matches"] == ["reference"] and alleles["AG"]["difference_from_reference"] is None
    assert alleles["GG"]["matches"] == ["catalog A>G"]
    ids = {a["rna_support"]["evidence_set_id"] for a in result["observed_alleles"]}
    assert ids <= set(result["evidence_sets"])


def test_min_fragments_sets_the_support_threshold(bam):
    result = reconcile(bam, min_fragments=1)
    assert by_id(result)["catalog A>G"]["rna_status"] == "supported"
    result = reconcile(bam, min_fragments=5)
    assert by_id(result)["caller MNV"]["rna_status"] == "insufficient_RNA"


def test_unknown_reference_bases_leave_the_locus_unassessed(bam):
    # Two SNVs ten bases apart: without an annotated exon the bases between are unknown.
    result = reconcile(bam, {"a": [variant(1005, "G", "A")], "b": [variant(1015, "T", "C")]})
    assert result["status"] == "unassessed" and "reference bases unknown" in result["reason"]
    assert {c["rna_status"] for c in result["candidates"]} == {"unassessed"}


def test_unusable_interpretations_are_rejected(bam):
    with pytest.raises(ValueError, match="No interpretations"):
        reconcile(bam, {})
    with pytest.raises(ValueError, match="overlap each other"):
        reconcile(bam, {"x": [variant(1021, "AG", "GT"), variant(1022, "G", "C")]})
    with pytest.raises(ValueError, match="disagree on the reference base"):
        reconcile(bam, {"x": [variant(1021, "A", "G")], "y": [variant(1021, "C", "G")]})
    with pytest.raises(ValueError, match="one contig"):
        reconcile(bam, {"x": [variant(1021, "A", "G")], "y": [variant(1021, "A", "G", contig="2")]})


def test_command_reads_json_interpretations(bam, tmp_path, capsys):
    from isovar.cli.commands import run as isovar_cli

    data = dict(reference_name="GRCh38", interpretations={
        "caller MNV": [dict(contig="1", start=1021, ref="AG", alt="GT")],
        "catalog A>G": [dict(contig="1", start=1021, ref="A", alt="G")]})
    path, output = tmp_path / "interpretations.json", tmp_path / "result.json"
    path.write_text(json.dumps(data))
    isovar_cli(["allele-interpretations", "--input", str(path), "--bam", bam.filename.decode(),
                "--sample-id", "s", "--output", str(output), "--log-level", "WARNING"])
    result = json.loads(output.read_text())
    assert by_id(result)["caller MNV"]["rna_status"] == "supported"
    assert result["source"] == bam.filename.decode()
    data["interpretations"]["bad"] = [dict(contig="1", start=1021, ref="A")]
    path.write_text(json.dumps(data))
    with pytest.raises(SystemExit) as exit_info:
        isovar_cli(["allele-interpretations", "--input", str(path), "--bam", bam.filename.decode(),
                    "--sample-id", "s", "--output", str(output)])
    assert exit_info.value.code == 2 and "each variant needs exactly" in capsys.readouterr().err


# Original osteosarc RNA with the calls that different sources made at each locus.

def real_case(gene):
    import tempfile
    from pathlib import Path

    from examples import osteosarc_assembly_figures as figures
    from tests.data.osteosarc.expansion.references import reference_genome

    manifest = json.loads((figures.CORPUS / "manifest.json").read_text())
    case = next(c for c in manifest["cases"] if c["variant"]["gene"] == gene
                and c.get("reference", "GRCh38") == "GRCh38")
    genome = reference_genome(figures.CORPUS / "references" / "GRCh38", Path(tempfile.mkdtemp()))
    return figures.CORPUS / case["primary_bam"], genome


def run_real(gene, calls):
    path, genome = real_case(gene)
    interpretations = {name: [Variant(*v, ensembl=genome) for v in vs] for name, vs in calls.items()}
    with pysam.AlignmentFile(str(path)) as bam:
        return reconcile_allele_interpretations(
            interpretations, bam, sample_id="sid", source=path.name,
            read_collector=ReadCollector(use_secondary_alignments=False))


def test_ntf3_rna_carries_the_compound_allele_not_the_catalogued_snv():
    result = run_real("NTF3", {
        "catalog A>G": [("12", 5494381, "A", "G")],
        "mutect2 AG>GT": [("12", 5494381, "AG", "GT")],
        "strelka A>G + G>T": [("12", 5494381, "A", "G"), ("12", 5494382, "G", "T")]})
    candidates = by_id(result)
    assert result["locus"]["padding"] == [10, 10]
    assert candidates["catalog A>G"]["rna_status"] == "contradicted_by_informative_evidence"
    assert candidates["mutect2 AG>GT"]["rna_status"] == "supported"
    assert candidates["mutect2 AG>GT"]["indistinguishable_from"] == ["strelka A>G + G>T"]
    assert candidates["mutect2 AG>GT"]["rna_support"]["fragments"] == 9
    assert result["reference_allele"]["rna_support"]["fragments"] == 5


def test_glis3_rna_carries_the_catalogued_deletion_with_a_co_somatic_snv():
    result = run_real("GLIS3", {
        "catalog deletion": [("9", 3856149, "CTGATGTGG", "C")],
        "catalog deletion + strelka A>T": [("9", 3856149, "CTGATGTGG", "C"), ("9", 3856160, "A", "T")],
        "strelka A>T": [("9", 3856160, "A", "T")],
        "mutect2 deletion": [("9", 3856153, "TGTGGTGA", "T")]})
    status = {key: c["rna_status"] for key, c in by_id(result).items()}
    assert status == {"catalog deletion": "contradicted_by_informative_evidence",
                      "catalog deletion + strelka A>T": "supported",
                      "strelka A>T": "contradicted_by_informative_evidence",
                      "mutect2 deletion": "contradicted_by_informative_evidence"}
    compound = next(a for a in result["observed_alleles"] if a["matches"] == ["catalog deletion + strelka A>T"])
    assert compound["rna_support"]["fragments"] == 10
    # The same haplotype, written as one replacement of the reference.
    assert compound["difference_from_reference"] == dict(start=3856152, ref="ATGTGGTGA", alt="T")


def test_reference_allele_must_agree_with_the_annotation_and_inputs_are_strict(bam):
    from isovar.allele_interpretations import interpretations_from_dict

    with pytest.raises(ValueError, match="disagrees with the annotation"):
        run_real("NTF3", {"wrong reference": [("12", 5494381, "C", "G")]})
    with pytest.raises(ValueError, match="at least one variant"):
        reconcile(bam, {"empty": []})
    with pytest.raises(ValueError, match="exactly the keys"):
        interpretations_from_dict(dict(interpretations={}))
