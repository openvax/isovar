"""Cell barcodes and UMIs behind small-variant allele reads (#226)."""

import json
from types import SimpleNamespace

import pysam
import pytest
from varcode import Variant

from isovar import ReadCollector, cell_umi_allele_evidence, export_protein_hypotheses, run_isovar
from isovar.cell_evidence import CellUmiAlleles
from isovar.rna_evidence import CELL_UMI_FIELDS
from tests.test_sv_rna import record, write_bam
from tests.testing_helpers import complete_umis

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
    assert (alt["reads"], alt["umis"], alt["cells"]) == (4, 3, 2)
    assert (complete_umis(alt), (alt["cells"] if alt["cells_complete"] else None)) == (3, 2)
    assert (ref["umis"], ref["cells"]) == (2, 2)
    assert result["alleles"]["other"]["reads"] == 0
    # Cell C1 has a reference and an alternate label.
    assert result["cells_with_ref_and_alt"] == 1
    protein, = result["protein_hypotheses"]
    assert (protein["reads"], protein["cells"], (protein["cells"] if protein["cells_complete"] else None)) == (2, 2, 2)
    assert set(alt) == {"reads", "fragments", *CELL_UMI_FIELDS}
    assert json.loads(json.dumps(result)) == result


def test_missing_umis_and_unknown_libraries_leave_complete_counts_null(tmp_path):
    reads = ALT + [read("a5", "G", "a", "C4"), read("a6", "G", "c", "C5", "U5")]
    alt = evidence(tmp_path, reads)["alleles"]["alt"]
    assert alt["unlabeled_reads"] == 1 and alt["unknown_library_reads"] == 1
    assert alt["label_statuses"]["missing_umi"] == 1
    # A barcode without a UMI still identifies its cell (C4); C5 is scoped to read group c.
    assert (alt["umis"], alt["cells"]) == (4, 4)
    assert complete_umis(alt) is None and (alt["cells"] if alt["cells_complete"] else None) is None


def test_one_label_in_two_read_groups_without_a_library_counts_twice(tmp_path):
    # Documented: without LB, labels stay local to their read group.
    alt = evidence(tmp_path, [read("x1", "G", "c", "C1", "U1"), read("x2", "G", "d", "C1", "U1")])["alleles"]["alt"]
    assert (alt["umis"], alt["cells"], alt["unknown_library_reads"]) == (2, 2, 2)


def test_reads_ending_at_an_insertion_are_labelled(tmp_path):
    insertion = Variant("1", 1020, "A", "ATT", normalize_contig_names=False)
    reads = [record("i%d" % i, "1", 1000, "20M2I20M", REFERENCE[:20] + "TT" + REFERENCE[20:]) for i in range(3)]
    reads.append(record("end", "1", 1000, "20M2I", REFERENCE[:20] + "TT"))
    for i, segment in enumerate(reads):
        segment.set_tag("CB", "C%d" % i)
        segment.set_tag("UB", "U%d" % i)
    path = write_bam(tmp_path / "rna.bam", reads, header=HEADER)
    with pysam.AlignmentFile(str(path)) as bam:
        read_evidence = ReadCollector().read_evidence_for_variant(insertion, bam)
        alt = CellUmiAlleles(bam, sample_id="s", source="rna.bam").evidence(
            insertion, read_evidence)["alleles"]["alt"]
    assert (alt["reads"], alt["umis"], alt["unlabeled_reads"]) == (4, 4, 0)


def test_reads_the_collector_would_reject_are_unresolved_and_warned(tmp_path, caplog):
    duplicates = [read("d%d" % i, "G", "a", "C%d" % i, "U%d" % i, flag=1024) for i in range(3)]
    path = write_bam(tmp_path / "rna.bam", ALT[:1] + duplicates, header=HEADER)
    with pysam.AlignmentFile(str(path)) as bam, caplog.at_level("WARNING", logger="isovar.cell_evidence"):
        read_evidence = ReadCollector(use_duplicate_reads=True).read_evidence_for_variant(VARIANT, bam)
        alt = CellUmiAlleles(bam, sample_id="s", source="rna.bam").evidence(
            VARIANT, read_evidence)["alleles"]["alt"]
    assert alt["label_statuses"]["metadata_unavailable"] == 3
    assert (alt["unlabeled_reads"], complete_umis(alt)) == (3, None)
    assert "no eligible record" in caplog.text


def test_mates_with_conflicting_labels_are_unresolved(tmp_path):
    reads = ALT + [read("pair", "G", "a", "C6", "U6", flag=99), read("pair", "G", "a", "C7", "U6", flag=147)]
    alt = evidence(tmp_path, reads)["alleles"]["alt"]
    assert alt["label_statuses"]["conflicting_template_labels"] == 2
    assert complete_umis(alt) is None


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
    assert (row.num_alt_reads, row.num_alt_fragments, row.num_alt_umis, row.num_alt_cells) == (4, 4, 3, 2)
    assert row.num_ref_cells == 2 and row.num_cells_with_ref_and_alt == 1
    assert row.num_alt_unlabeled_reads == row.num_alt_unknown_library_reads == 0
    assert row.alt_umis_complete and row.alt_cells_complete
    # No read has the other allele, so its counts are zero but not marked exact.
    assert row.num_other_umis == 0 and not row.other_umis_complete
    with pytest.raises(SystemExit) as exit_info:
        isovar_cli(["allele-counts", "--vcf", str(vcf), "--bam", str(path), "--cell-umi-labels",
                    "--output", str(output)])
    assert exit_info.value.code == 2 and "--sample-id" in capsys.readouterr().err
    with pytest.raises(SystemExit) as exit_info:
        isovar_cli(["allele-counts", "--vcf", str(vcf), "--bam", str(path), "--sample-id", "s",
                    "--output", str(output)])
    assert exit_info.value.code == 2 and "only apply with --cell-umi-labels" in capsys.readouterr().err
    output = tmp_path / "hypotheses.json"
    isovar_cli(["protein-hypotheses", "--vcf", str(vcf), "--bam", str(path), "--cell-umi-labels",
                "--sample-id", "s", "--output", str(output), "--log-level", "WARNING"])
    event, = json.loads(output.read_text())["events"]
    assert event["allele_support"]["alt"]["cells"] == 2 and event["allele_support"]["cells_with_ref_and_alt"] == 1


def test_allele_counts_still_streams_without_the_flag(tmp_path, monkeypatch):
    import types
    from isovar.cli import rna_args
    from isovar.cli.isovar_allele_counts import parser

    seen = []
    monkeypatch.setattr(rna_args, "allele_counts_dataframe", lambda pairs: seen.append(pairs) or "table")
    path = write_bam(tmp_path / "rna.bam", ALT, header=HEADER)
    args = parser.parse_args(["--variant", "1", "1021", "A", "G", "--genome", "GRCh38", "--bam", str(path)])
    assert rna_args.allele_counts_dataframe_from_args(args) == "table"
    assert isinstance(seen[0], types.GeneratorType)


def test_untagged_reads_are_warned_about(tmp_path, caplog):
    path = write_bam(tmp_path / "rna.bam", [read("u%d" % i, "G") for i in range(2)], header=HEADER)
    from isovar.cell_evidence import warn_if_unlabelled
    with pysam.AlignmentFile(str(path)) as bam, caplog.at_level("WARNING", logger="isovar.cell_evidence"):
        read_evidence = ReadCollector().read_evidence_for_variant(VARIANT, bam)
        warn_if_unlabelled(CellUmiAlleles(bam, sample_id="s", source="rna.bam").evidence(
            VARIANT, read_evidence)["alleles"].values())
    assert "No read carried a usable CB" in caplog.text


def test_protein_export_warns_about_reads_without_an_eligible_record(tmp_path, caplog):
    # Reads collected with duplicates allowed, then labelled with the default collector.
    duplicates = [read("d%d" % i, "G", "a", "C%d" % i, "U%d" % i, flag=1024) for i in range(3)]
    path = write_bam(tmp_path / "rna.bam", duplicates, header=HEADER)
    with pysam.AlignmentFile(str(path)) as bam, caplog.at_level("WARNING", logger="isovar.cell_evidence"):
        results = run_isovar([VARIANT], bam, read_collector=ReadCollector(use_duplicate_reads=True))
        event, = export_protein_hypotheses(results, sample_id="s", source="rna.bam",
                                           cell_umi_alignment_file=bam)["events"]
    assert event["allele_support"]["alt"]["label_statuses"] == {"metadata_unavailable": 3}
    assert "3 read(s) at 1:1021-1021 had no eligible record" in caplog.text
    # The reads do carry barcodes, so the not-single-cell warning would mislead.
    assert "No read carried a usable CB" not in caplog.text


def test_single_cell_ont_rna_reports_cells_behind_the_alt_allele_and_protein(tmp_path):
    from examples import osteosarc_assembly_figures as figures
    from isovar import export_protein_hypotheses, run_isovar
    from tests.data.osteosarc.expansion.references import load_reference, reference_genome

    manifest = json.loads((figures.CORPUS / "manifest.json").read_text())
    case = next(c for c in manifest["cases"] if c["case_id"] == "26-NME1-chr17-51154443-53f498a544883d51")
    reference_dir = figures.CORPUS / "references" / "GRCh38"
    reference_manifest, _ = load_reference(reference_dir)
    genome = reference_genome(reference_dir, tmp_path)
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
    assert (alt["reads"], alt["umis"], alt["cells"]) == (48, 48, 44)
    # This selected fixture declares no read groups, so the library is unknown and counts are lower bounds.
    assert alt["unknown_library_reads"] == 48 and (alt["cells"] if alt["cells_complete"] else None) is None
    assert cells["cells_with_ref_and_alt"] == 5
    event, = export["events"]
    # The export's supports are the same records, plus a scope and evidence set ID.
    assert {k: event["allele_support"]["alt"][k] for k in alt} == alt
    top = event["protein_hypotheses"][0]["rna_support"]
    assert {k: top[k] for k in cells["protein_hypotheses"][0]} == cells["protein_hypotheses"][0]
    assert top["cells"] == 38 and top["evidence_set_id"] in export["evidence_sets"]
