# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Every command's CSV identifies variants the same way and names counts as counts."""

import re

import pandas as pd
import pytest

from varcode import Variant

from isovar import AlleleRead, IsovarResult, VariantSequence
from isovar.cli import commands
from isovar.dataframe_helpers import isovar_results_to_dataframe, variant_sequences_generator_to_dataframe

from .testing_helpers import data_path

VCF = data_path("data/b16.f10/b16.vcf")
BAM = data_path("data/b16.f10/b16.combined.sorted.bam")
VARIANT_COLUMNS = ["variant", "chr", "pos", "ref", "alt"]


def run(tmp_path, command, *options, vcf=VCF):
    output = tmp_path / (command + ".csv")
    args = [command, "--vcf", vcf, "--output", str(output), "--log-level", "ERROR", *options]
    if command != "reference-contexts":
        args += ["--bam", BAM]
    commands.run(args)
    return output


@pytest.mark.parametrize("command", [
    "run", "protein-sequences", "translations", "variant-sequences",
    "reference-contexts", "allele-counts", "allele-reads", "variant-reads",
])
def test_every_command_uses_the_same_variant_key_with_input_contig_names(tmp_path, command):
    df = pd.read_csv(run(tmp_path, command))
    assert list(df.columns[:5]) == VARIANT_COLUMNS
    # The VCF names contigs "chr13" etc.; keep them rather than a normalized "13".
    assert set(df["chr"]) <= {"chr13", "chr4", "chr9", "chrX"} and df["chr"].str.startswith("chr").all()
    assert df["variant"].str.startswith("chr").all()


def test_count_columns_are_named_as_counts(tmp_path):
    proteins = pd.read_csv(run(tmp_path, "protein-sequences"))
    assert {
        "gene_names", "transcript_ids", "transcript_names",
        "num_translations", "num_supporting_reads", "num_supporting_fragments",
    } <= set(proteins.columns)
    assert not {"translations", "gene", "known_variants"} & set(proteins.columns)

    sequences = pd.read_csv(run(tmp_path, "variant-sequences"))
    assert {
        "gene_names", "transcript_ids", "num_supporting_reads", "num_supporting_fragments",
    } <= set(sequences.columns)
    assert not {"reads", "num_reads", "num_fragments", "compatible_transcript_ids", "gene"} & set(sequences.columns)
    # This command does not classify transcripts; when present, IDs are listed, not counted.
    assert sequences["transcript_ids"].isna().all()
    reads = [AlleleRead("AC", "G", "TT", "read%d" % i, compatible_transcript_ids={"T2", "T1"}) for i in range(2)]
    sequence = VariantSequence(prefix="AC", alt="G", suffix="TT", reads=reads)
    frame = variant_sequences_generator_to_dataframe([(Variant("1", 3, "A", "G"), [sequence])])
    assert frame["transcript_ids"].tolist() == ["T1;T2"]
    assert frame["num_supporting_reads"].tolist() == frame["num_supporting_fragments"].tolist() == [2]

    translations = pd.read_csv(run(tmp_path, "translations"))
    assert {
        "gene_names", "in_frame_cdna_sequence", "transcript_ids", "transcript_names",
        "num_supporting_reads", "num_supporting_fragments",
    } <= set(translations.columns)
    assert not {
        "variant_orf", "reference_context", "reference_transcript_names",
        "untrimmed_variant_sequence_read_count",
    } & set(translations.columns)

    contexts = pd.read_csv(run(tmp_path, "reference-contexts"))
    assert {"gene_names", "transcript_ids", "transcript_names"} <= set(contexts.columns)
    assert not {"transcripts", "gene"} & set(contexts.columns)

    counts = pd.read_csv(run(tmp_path, "allele-counts"))
    assert {
        "gene_names", "num_ref_reads", "num_alt_reads", "num_other_reads",
        "num_ref_fragments", "num_alt_fragments", "num_other_fragments",
    } <= set(counts.columns)


def test_run_record_columns_have_logical_non_alphabetical_order():
    cols = IsovarResult.record_columns()
    props = IsovarResult.record_property_names()
    all_prefixed = {
        key for key in dir(IsovarResult)
        if key.startswith(("num_", "fraction_", "ratio_", "has_"))
    }
    assert set(IsovarResult.RECORD_PROPERTY_NAMES) == all_prefixed
    assert set(props) == all_prefixed
    # Locus counts precede fractions, which precede ratios, which precede boolean flags.
    assert cols.index("num_total_reads") < cols.index("num_alt_reads")
    assert cols.index("num_other_fragments") < cols.index("fraction_alt_reads")
    assert cols.index("fraction_other_fragments") < cols.index("ratio_alt_to_other_reads")
    assert cols.index("ratio_other_to_ref_fragments") < cols.index("has_mutant_protein_sequence_from_rna")


def test_default_output_filenames_tsv_delimiter_and_stdout(tmp_path, monkeypatch, capsys):
    expected_defaults = {
        "run": "isovar-results.csv",
        "protein-sequences": "isovar-protein-sequences.csv",
        "translations": "isovar-translations.csv",
        "variant-sequences": "isovar-variant-sequences.csv",
        "variant-reads": "isovar-variant-reads.csv",
        "allele-reads": "isovar-allele-reads.csv",
        "allele-counts": "isovar-allele-counts.csv",
        "reference-contexts": "isovar-reference-contexts.csv",
    }
    for name, expected in expected_defaults.items():
        parser = commands.load_subcommand_parser(name)
        assert parser.get_default("output") == expected, name

    tsv = tmp_path / "counts.tsv"
    commands.run(["allele-counts", "--vcf", VCF, "--bam", BAM, "--output", str(tsv), "--log-level", "ERROR"])
    raw = tsv.read_bytes()
    assert b"\r\n" not in raw
    header = raw.decode().splitlines()[0]
    assert "\t" in header and "," not in header
    assert list(pd.read_csv(tsv, sep="\t").columns[:5]) == VARIANT_COLUMNS

    tsv_gz = tmp_path / "counts.tsv.gz"
    commands.run(["allele-counts", "--vcf", VCF, "--bam", BAM, "--output", str(tsv_gz), "--log-level", "ERROR"])
    assert list(pd.read_csv(tsv_gz, sep="\t").columns[:5]) == VARIANT_COLUMNS

    monkeypatch.chdir(tmp_path)
    commands.run(["allele-counts", "--vcf", VCF, "--bam", BAM, "--output", "-", "--log-level", "INFO"])
    captured = capsys.readouterr()
    assert not (tmp_path / "-").exists()
    assert captured.out.splitlines()[0].startswith("variant,chr,pos,ref,alt,")
    assert "Wrote 4 rows to stdout" in captured.err


def test_allele_reads_reconciles_competing_secondary_and_mate_labels():
    from isovar.dataframe_helpers import allele_counts_dataframe, read_evidence_to_allele_reads_dataframe
    from isovar.read_collector import ReadCollector
    from .mock_objects import MockAlignmentFile, make_pysam_read

    variant = Variant("1", 5, "A", "C", normalize_contig_names=False)
    # Read "conflict" has primary=ALT and secondary=REF -> counted as "other".
    primary_alt = make_pysam_read(seq="TTTTCTTTT", cigar="9M", name="conflict", reference_start=0)
    primary_alt.flag = 0
    secondary_ref = make_pysam_read(seq="TTTTATTTT", cigar="9M", name="conflict", reference_start=0)
    secondary_ref.flag = 256
    # Paired-end overlapping mates carrying ALT -> merged as R1+R2.
    mate1 = make_pysam_read(seq="TTTTCT", cigar="6M", name="pair1", reference_start=0)
    mate1.flag = 65
    mate2 = make_pysam_read(seq="TTCTTTT", cigar="7M", name="pair1", reference_start=2)
    mate2.flag = 129
    # Unmerged R1 carrying REF.
    r1_ref = make_pysam_read(seq="TTTTATTTT", cigar="9M", name="pair2", reference_start=0)
    r1_ref.flag = 65

    bam = MockAlignmentFile(["1"], [primary_alt, secondary_ref, mate1, mate2, r1_ref])
    evidence = ReadCollector().read_evidence_for_variant(variant, bam)
    counts = allele_counts_dataframe([(variant, evidence)]).iloc[0]
    reads_df = read_evidence_to_allele_reads_dataframe([(variant, evidence)])

    assert (counts.num_ref_reads, counts.num_alt_reads, counts.num_other_reads) == (1, 2, 1)
    assert (counts.num_ref_fragments, counts.num_alt_fragments, counts.num_other_fragments) == (1, 1, 1)

    conflict_rows = reads_df[reads_df["name"] == "conflict"]
    assert len(conflict_rows) == 2
    assert set(conflict_rows["allele_class"]) == {"other"}
    assert set(conflict_rows["is_secondary"]) == {False, True}
    assert set(conflict_rows["mate"]) == {"unpaired"}

    pair1_row = reads_df[reads_df["name"] == "pair1"].iloc[0]
    assert pair1_row["allele_class"] == "alt"
    assert pair1_row["mate"] == "R1+R2"
    assert bool(pair1_row["is_secondary"]) is False

    pair2_row = reads_df[reads_df["name"] == "pair2"].iloc[0]
    assert pair2_row["allele_class"] == "ref"
    assert pair2_row["mate"] == "R1"
    assert bool(pair2_row["is_secondary"]) is False


def test_variant_sequences_cli_uses_length_ladder_by_default(monkeypatch):
    from isovar.cli import isovar_variant_sequences as cli_vs
    from isovar.read_evidence import ReadEvidence

    variant = Variant("1", 50, "A", "G", normalize_contig_names=False)
    # Two noisy long reads agree over a 65-bp core (32 + 1 + 32) around the variant,
    # but disagree in the outer bases (positions 33..45) required by a fixed 90-nt window.
    core_left = "ACGT" * 8
    core_right = "TGCA" * 8
    reads = [
        AlleleRead("A" * 20 + core_left, "G", core_right + "A" * 20, "ont1"),
        AlleleRead("C" * 20 + core_left, "G", core_right + "C" * 20, "ont2"),
    ]
    evidence = ReadEvidence(50, "A", "G", ref_reads=[], alt_reads=reads, other_reads=[])
    monkeypatch.setattr(cli_vs, "read_evidence_generator_from_args", lambda args: [(variant, evidence)])

    ladder_args = cli_vs.parser.parse_args(["--bam", "unused.bam"])
    ladder_df = cli_vs.variant_sequences_dataframe_from_args(ladder_args)
    assert len(ladder_df) >= 1
    assert any(seq in (core_left + "G" + core_right) and len(seq) == 62 for seq in ladder_df["sequence"])

    fixed_args = cli_vs.parser.parse_args(["--bam", "unused.bam", "--variant-sequence-length", "90"])
    fixed_df = cli_vs.variant_sequences_dataframe_from_args(fixed_args)
    assert fixed_df.empty


def test_floats_are_written_with_bounded_precision(tmp_path):
    text = run(tmp_path, "run").read_text()
    assert not re.search(r"\d\.\d{7,}", text)


def test_empty_results_keep_their_header(tmp_path):
    header_only = tmp_path / "header-only.vcf"
    header_only.write_text("".join(line for line in open(VCF) if line.startswith("#")))
    df = pd.read_csv(run(tmp_path, "run", vcf=str(header_only)))
    assert df.empty and list(df.columns) == IsovarResult.record_columns()
    assert list(isovar_results_to_dataframe([]).columns) == IsovarResult.record_columns()
