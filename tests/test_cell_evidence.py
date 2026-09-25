"""Cell barcodes and UMIs behind small-variant allele reads (#226)."""

import json
from types import SimpleNamespace

import pysam
import pytest
from varcode import Variant

from isovar import ReadCollector, cell_umi_allele_evidence
from isovar.cell_evidence import CellUmiAlleles
from tests.test_sv_rna import record, write_bam

REFERENCE = "CATCGGATCCTAGGCTTACA" + "AG" + "TCCGATTGCAGTCAGGAT"
VARIANT = Variant("1", 1021, "A", "G", normalize_contig_names=False)
# Read groups a and b share one declared library; c declares none.
HEADER = pysam.AlignmentHeader.from_dict(dict(
    SQ=[dict(SN="1", LN=20000), dict(SN="2", LN=20000)],
    RG=[dict(ID="a", SM="tumor", LB="lib1"), dict(ID="b", SM="tumor", LB="lib1"), dict(ID="c", SM="tumor")]))


def read(name, allele, group="a", cell=None, umi=None, flag=0):
    segment = record(name, "1", 1000, "40M", REFERENCE[:20] + allele + REFERENCE[21:], flag=flag, group=group)
    for tag, value in (("CB", cell), ("UB", umi)):
        if value is not None:
            segment.set_tag(tag, value)
    return segment


def evidence(tmp_path, reads, proteins=()):
    path = write_bam(tmp_path / "rna.bam", reads, header=HEADER)
    with pysam.AlignmentFile(str(path)) as bam:
        read_evidence = ReadCollector().read_evidence_for_variant(VARIANT, bam)
        labels = CellUmiAlleles(bam, sample_id="s", source="rna.bam")
        protein_reads = [[r for r in read_evidence.alt_reads if r.name in names] for names in proteins]
        return labels.evidence(VARIANT, read_evidence,
                               [SimpleNamespace(supporting_reads=reads) for reads in protein_reads])


ALT = [read("a1", "G", "a", "C1", "U1"), read("a2", "G", "b", "C1", "U1"),
       read("a3", "G", "a", "C2", "U2"), read("a4", "G", "a", "C2", "U3")]
REF = [read("r1", "A", "a", "C1", "U9"), read("r2", "A", "a", "C3", "U10")]


def test_labels_and_cells_are_counted_per_allele_within_the_declared_library(tmp_path):
    result = evidence(tmp_path, ALT + REF, proteins=[{"a1", "a3"}])
    alt, ref = result["alleles"]["alt"], result["alleles"]["ref"]
    assert result["policy"] == "isovar.cell_umi_labels.v1"
    # a1 and a2 share library lib1 and the C1/U1 label, so they count once.
    assert (alt["segments"], alt["observed_labels"], alt["observed_cells"]) == (4, 3, 2)
    assert (alt["complete_label_count"], alt["complete_cell_count"]) == (3, 2)
    assert (ref["observed_labels"], ref["observed_cells"]) == (2, 2)
    assert result["alleles"]["other"]["segments"] == 0
    # Cell C1 has a reference and an alternate label.
    assert result["cells_with_ref_and_alt"] == 1
    protein, = result["protein_hypotheses"]
    assert (protein["segments"], protein["observed_cells"], protein["complete_cell_count"]) == (2, 2, 2)
    assert alt["evidence_set_id"].startswith("evidence_") and "segment_ids" not in alt
    assert json.loads(json.dumps(result)) == result


def test_missing_umis_and_unknown_libraries_make_counts_lower_bounds(tmp_path):
    reads = ALT + [read("a5", "G", "a", "C4"), read("a6", "G", "c", "C5", "U5")]
    alt = evidence(tmp_path, reads)["alleles"]["alt"]
    assert alt["unresolved_segments"] == 1 and alt["unknown_library_segments"] == 1
    assert alt["status_counts"]["missing_umi"] == 1
    assert (alt["observed_labels"], alt["observed_cells"]) == (4, 3)
    assert alt["complete_label_count"] is None and alt["complete_cell_count"] is None


def test_mates_with_conflicting_labels_are_unresolved(tmp_path):
    reads = ALT + [read("pair", "G", "a", "C6", "U6", flag=99), read("pair", "G", "a", "C7", "U6", flag=147)]
    alt = evidence(tmp_path, reads)["alleles"]["alt"]
    assert alt["status_counts"]["conflicting_template_labels"] == 2
    assert alt["complete_label_count"] is None


def test_command_adds_cell_columns_and_requires_a_sample(tmp_path, capsys):
    import pandas as pd
    from isovar.cli.commands import run as isovar_cli

    path = write_bam(tmp_path / "rna.bam", ALT + REF, header=HEADER)
    vcf = tmp_path / "variants.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                   "1\t1021\t.\tA\tG\t.\tPASS\t.\n")
    output = tmp_path / "counts.csv"
    isovar_cli(["allele-counts", "--vcf", str(vcf), "--bam", str(path), "--cell-umi-labels",
                "--sample-id", "s", "--output", str(output), "--log-level", "WARNING"])
    row = pd.read_csv(output).iloc[0]
    assert (row.num_alt_cells, row.num_alt_cell_umi_labels, row.num_ref_cells) == (2, 3, 2)
    assert row.num_cells_with_ref_and_alt == 1 and row.num_alt_unlabeled_segments == 0
    with pytest.raises(SystemExit) as exit_info:
        isovar_cli(["allele-counts", "--vcf", str(vcf), "--bam", str(path), "--cell-umi-labels",
                    "--output", str(output)])
    assert exit_info.value.code == 2 and "--sample-id" in capsys.readouterr().err


def test_single_cell_ont_rna_reports_cells_behind_the_alt_allele_and_protein():
    import tempfile
    from pathlib import Path

    from examples import osteosarc_assembly_figures as figures
    from isovar import export_protein_hypotheses, run_isovar
    from tests.data.osteosarc.expansion.references import load_reference, reference_genome

    manifest = json.loads((figures.CORPUS / "manifest.json").read_text())
    case = next(c for c in manifest["cases"] if c["case_id"] == "26-NME1-chr17-51154443-53f498a544883d51")
    reference_dir = figures.CORPUS / "references" / "GRCh38"
    reference_manifest, _ = load_reference(reference_dir)
    genome = reference_genome(reference_dir, Path(tempfile.mkdtemp()))
    r = case["variant"]
    variant = Variant("17", r["pos"], r["ref"], r["alt"], ensembl=genome)
    collector = ReadCollector(use_secondary_alignments=False)
    with pysam.AlignmentFile(str(figures.CORPUS / case["primary_bam"])) as bam:
        results = run_isovar([variant], bam, read_collector=collector,
                             transcript_id_whitelist=set(reference_manifest["variant_transcripts"][r["variant_id"]]))
        cells, = cell_umi_allele_evidence(results, bam, sample_id="sid", source="nme1", read_collector=collector)
        export = export_protein_hypotheses(results, sample_id="sid", source="nme1",
                                           cell_umi_alignment_file=bam, read_collector=collector)
    alt = cells["alleles"]["alt"]
    assert (alt["segments"], alt["observed_labels"], alt["observed_cells"]) == (48, 48, 44)
    # This selected fixture declares no read groups, so the library is unknown and counts are lower bounds.
    assert alt["unknown_library_segments"] == 48 and alt["complete_cell_count"] is None
    assert cells["cells_with_ref_and_alt"] == 5
    event, = export["events"]
    assert event["cell_umi_evidence"]["alleles"]["alt"] == alt
    top = event["protein_hypotheses"][0]["rna_support"]["cell_umi_support"]
    assert top == cells["protein_hypotheses"][0] and top["observed_cells"] == 38
