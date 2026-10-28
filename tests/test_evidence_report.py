"""Exon-footprint, denominator and command checks for evidence reports."""

import csv
from pathlib import Path

import pytest

from isovar.evidence_report import (
    changed_reference_interval, mutant_ninemers, overlapping_catalog_orfs, overlaps_exon,
    reference_ninemer_matches, self_novelty, summarize,
)


class TinyExons:
    def __init__(self, intervals):
        self.intervals = intervals

    def exon_ids_at_locus(self, contig, first, end=None):
        assert contig == "1"
        last = first if end is None else end
        return ["exon"] if any(a <= last and b >= first
                               for a, b in self.intervals) else []


def test_exon_funnel_uses_changed_bases_and_both_insertion_flanks():
    genome = TinyExons([(103, 104), (200, 200)])
    assert changed_reference_interval(100, "ATTTTTA", "A") == (101, 106, False)
    assert overlaps_exon(genome, "chr1", 100, "ATTTTTA", "A")
    assert not overlaps_exon(genome, "1", 100, "AT", "A")
    assert changed_reference_interval(199, "A", "AT") == (199, 200, True)
    assert overlaps_exon(genome, "1", 199, "A", "AT")
    assert changed_reference_interval(200, "AT", "T") == (200, 200, False)
    assert overlaps_exon(genome, "1", 200, "AT", "T")


def test_ncorf_overlap_accepts_chr_prefixed_vcf_contigs():
    index = {"1": ([99], [110], [(99, 110, "orf-1")])}
    assert overlapping_catalog_orfs(index, "chr1", 100, 101) == ["orf-1"]


def test_nested_counts_and_top_context_aa_denominators():
    calls = [
        dict(call_id="s", **{"class": "SNV"}, ref="A", alt="G", rna_status="alternate_observed",
             alt_reads="2", protein_count="1"),
        dict(call_id="d", **{"class": "indel"}, ref="AT", alt="A",
             rna_status="alternate_observed", alt_reads="3", protein_count="1"),
        dict(call_id="i", **{"class": "indel"}, ref="A", alt="AT",
             rna_status="alternate_observed", alt_reads="1", protein_count="0"),
        dict(call_id="z", **{"class": "SNV"}, ref="A", alt="C",
             rna_status="no_alignment_overlap", alt_reads="0", protein_count="0"),
    ]
    annotations = {
        "s": dict(exon_overlap=True, reference_coding_effect=True,
                  top_effect_type="Substitution"),
        "d": dict(exon_overlap=True, top_effect_type="FrameShift"),
        "i": dict(exon_overlap=False, top_effect_type="Intronic"),
        "z": dict(exon_overlap=False, top_effect_type="Intergenic"),
    }
    proteins = [
        dict(call_id="s", rank="1", protein_sequence="ABCDEFGHIJK",
             mutation_start_idx="5", mutation_end_idx="6", comparison_status="amino_acid_divergent",
             focal_comparison_status="amino_acid_divergent", changed_residue_count=1),
        dict(call_id="d", rank="1", protein_sequence="ABCDEFGHIJK",
             mutation_start_idx="5", mutation_end_idx="5", comparison_status="amino_acid_divergent",
             focal_comparison_status="amino_acid_divergent", changed_residue_count=1),
    ]
    result = summarize(calls, annotations, proteins)
    assert [result["counts"][key] for key in (
        "total", "exonic", "exonic_alt", "exonic_product", "frameshift",
        "frameshift_alt")] == [4, 2, 2, 2, 1, 1]
    assert result["aa_counts"] == [1, 1, None]
    assert result["marked_aa_counts"] == [1, 0, None]
    assert result["reads_by_class"] == {"SNV": [2, 0], "deletion": [3], "insertion": [1]}
    assert result["exonic_reads_by_class"] == {"SNV": [2], "deletion": [3]}
    pytest.importorskip("reportlab")
    from isovar.evidence_pdf import (
        EvidenceFunnel, ProteinComparison, altered_aa_figures, read_support_figures,
    )

    assert [n for _, n in EvidenceFunnel(result["counts"]).stages] == [4, 2, 2, 2, 2]
    assert [n for _, n in EvidenceFunnel(stages=(("SV", 3), ("exonic", 2))).stages] == [3, 2]
    with pytest.raises(ValueError, match="nested counts"):
        EvidenceFunnel(stages=(("SV", 2), ("exonic", 3)))
    with pytest.raises(ValueError, match="no amino-acid divergence"):
        ProteinComparison("ABCDE", "ABCDE", 2, 3, point_change=True)
    read_figures = read_support_figures(result["reads_by_class"])
    assert read_figures[0].values["0"] == (1, 0, 0)
    assert "0" not in read_figures[1].bins
    exonic_figures = read_support_figures(result["exonic_reads_by_class"],
                                          scope="exon-involved variant")
    assert exonic_figures[0].values["0"] == (0, 0)
    assert "0" not in exonic_figures[1].bins
    aa_figures = altered_aa_figures(result["aa_counts"])
    assert aa_figures[0].values["1"] == 2
    assert "1" not in aa_figures[1].bins
    with pytest.raises(ValueError, match="Exon overlap"):
        summarize(calls, {**annotations, "s": {}}, proteins)


def test_mutation_overlapping_reference_ninemers_include_deletion_junction(tmp_path):
    protein = dict(protein_sequence="ABCDEFGHIJKL", mutation_start_idx="5",
                   mutation_end_idx="5")
    kmers = mutant_ninemers(protein)
    assert kmers and all("F" in item for item in kmers)
    fasta = tmp_path / "reference.fasta"
    fasta.write_text(">self\nABCDEFGHI\n>other\nZZZZZZZZZ\n")
    matches = reference_ninemer_matches([protein], fasta)
    assert "ABCDEFGHI" in matches
    novel, total = self_novelty(protein, matches)
    assert novel == total - 1


def test_unresolved_chromosomes_keep_unmeasured_counts_and_do_not_enter_depth_bins(tmp_path):
    import json
    import pysam
    from pyensembl import EnsemblRelease
    from varcode import Variant
    from isovar import ReadCollector, ProteinSequenceCreator
    from isovar.cli.isovar_evidence_report import collect, write_per_call

    path = tmp_path / "empty.bam"
    with pysam.AlignmentFile(path, "wb", header={"SQ": [
            {"SN": "chr1", "LN": 1000}, {"SN": "CHR1", "LN": 1000}]}):
        pass
    pysam.index(str(path))
    genome = EnsemblRelease(87)
    variants = [Variant(name, 100 + i, "A", "G", genome, convert_ucsc_contig_names=True)
                for i, name in enumerate(("chr1", "1", "2"))]
    with pysam.AlignmentFile(path) as bam:
        calls, annotations, proteins = collect(variants, bam, ReadCollector(), ProteinSequenceCreator(),
                                               sample_id="fixture", source="empty.bam")
    assert [c["rna_status"] for c in calls] == ["no_alignment_overlap", "unresolved_chromosome",
                                               "unresolved_chromosome"]
    assert [c["contig_resolution"]["status"] for c in calls[1:]] == ["ambiguous", "missing"]
    assert all(c["alt_reads"] is None and c["protein_count"] is None for c in calls[1:])
    result = summarize(calls, annotations, proteins)
    assert result["counts"]["total"] == 3 and result["counts"]["unresolved_chromosomes"] == 2
    assert result["reads_by_class"] == {"SNV": [0]}
    output = tmp_path / "calls.tsv"
    write_per_call(output, calls, annotations, result,
                   dict(sample_id="fixture", dna_source="calls", rna_source="empty.bam"), {})
    with output.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert rows[0]["alt_reads"] == "0" and rows[1]["alt_reads"] == ""
    assert json.loads(rows[1]["contig_resolution"])["candidates"] == ["CHR1", "chr1"]


def test_vcf_bam_command_writes_per_call_and_visual_pdf(tmp_path):
    pytest.importorskip("reportlab")
    pypdf = pytest.importorskip("pypdf")
    from isovar.cli.isovar_evidence_report import run

    data = Path(__file__).parent / "data/somatic-variant-with-2-supporting-rna-reads"
    output = tmp_path / "report.pdf"
    run(["--vcf", str(data / "single_variant.vcf"),
         "--bam", str(data / "rna.14.105849746.G.A.2-alt.sorted.bam"),
         "--no-merge-overlapping-fragments", "--sample-id", "fixture",
         "--output", str(output)])
    rows = list(csv.DictReader(output.with_suffix(".tsv").open(), delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["sample_id"] == "fixture"
    assert rows[0]["alt_reads"] == "2"
    assert rows[0]["exon_overlap"] in ("True", "False")
    pdf = pypdf.PdfReader(output)
    content = "\n".join(page.extract_text() for page in pdf.pages)
    assert "Variant-to-RNA funnel" in content
    assert "ALT reads per variant: positive support only" in content
    assert "ALT reads per exon-involved variant: all records" in content
    assert "Verified residue differences: excluding one-residue calls" in content
    import json
    receipt = json.loads(output.with_suffix('.json').read_text())
    assert receipt['schema'] == 'isovar.evidence_report.v1'
    native = receipt['calls'][0]['protein_hypotheses']
    assert native['sample_id'] == 'fixture'
    assert native['events'][0]['allele_support']['alt']['reads'] == 2
    assert 'evidence_sets' in native


def test_report_rejects_output_collisions_before_loading_inputs(tmp_path, capsys):
    from isovar.cli.isovar_evidence_report import run
    source = tmp_path / "source.vcf"
    source.write_text("preserved input")
    with pytest.raises(SystemExit):
        run(["--vcf", str(source), "--bam", "unopened.bam", "--sample-id", "s",
             "--output", str(tmp_path / "report.pdf"), "--evidence-json", str(source)])
    assert source.read_text() == "preserved input"
    assert "must not overwrite input files" in capsys.readouterr().err
