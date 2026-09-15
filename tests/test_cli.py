# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import tempfile
from os import remove
from os.path import getsize, exists
import pandas as pd

from .testing_helpers import data_path

from isovar.cli.isovar_translations import run as isovar_translations
from isovar.cli.isovar_allele_counts import run as isovar_allele_counts
from isovar.cli.isovar_allele_reads import run as isovar_allele_reads
from isovar.cli.isovar_protein_sequences import run as isovar_protein_sequences
from isovar.cli.isovar_reference_contexts import run as isovar_reference_contexts
from isovar.cli.isovar_variant_reads import run as isovar_variant_reads
from isovar.cli.isovar_variant_sequences import run as isovar_variant_sequences
from isovar.cli.isovar_main import run as isovar_main
from isovar.cli.rna_args import (
    allele_counts_dataframe_from_args,
    make_rna_reads_arg_parser,
    variants_reads_dataframe_from_args,
)

vcf_args = [
    "--vcf",
    data_path("data/b16.f10/b16.vcf")
]

args_with_bam = vcf_args + [
    "--bam",
    data_path("data/b16.f10/b16.combined.sorted.bam")
]


def run_cli_fn(fn, include_bam_in_args=True, return_dataframe=False):
    with tempfile.NamedTemporaryFile(delete=False) as f:
        output_path = f.name
    assert not exists(output_path) == 0
    output_args = ["--output", output_path]
    if include_bam_in_args:
        args = args_with_bam + output_args
    else:
        args = vcf_args + output_args
    fn(args)
    assert getsize(output_path) > 0
    if return_dataframe:
        df = pd.read_csv(output_path)
    remove(output_path)
    if return_dataframe:
        return df


def test_cli_allele_counts():
    df = run_cli_fn(isovar_allele_counts, return_dataframe=True)
    args = make_rna_reads_arg_parser().parse_args(args_with_bam)
    expected_df = allele_counts_dataframe_from_args(args)
    assert list(df.columns) == list(expected_df.columns)
    assert "ref_reads" not in df.columns
    assert len(df) == 4
    assert df[["num_ref_reads", "num_alt_reads", "num_other_reads"]].to_dict(orient="records") == \
        expected_df[["num_ref_reads", "num_alt_reads", "num_other_reads"]].to_dict(orient="records")
    assert df[["num_ref_fragments", "num_alt_fragments", "num_other_fragments"]].to_dict(orient="records") == \
        expected_df[["num_ref_fragments", "num_alt_fragments", "num_other_fragments"]].to_dict(orient="records")


def test_cli_allele_reads():
    df = run_cli_fn(isovar_allele_reads, return_dataframe=True)
    assert set(["prefix", "allele", "suffix", "name", "sequence", "gene"]).issubset(df.columns)
    assert "ref_reads" not in df.columns
    assert len(df) == 294


def test_cli_reference_contexts():
    run_cli_fn(isovar_reference_contexts, include_bam_in_args=False)


def test_cli_protein_sequences():
    run_cli_fn(isovar_protein_sequences)


def test_cli_translations():
    run_cli_fn(isovar_translations)


def test_cli_variant_reads():
    df = run_cli_fn(isovar_variant_reads, return_dataframe=True)
    assert set(["prefix", "allele", "suffix", "name", "sequence", "gene"]).issubset(df.columns)
    assert "ref_reads" not in df.columns
    assert len(df) == 42


def test_cli_variant_sequences():
    run_cli_fn(isovar_variant_sequences)

def test_cli_main():
    run_cli_fn(isovar_main)


def test_variant_reads_dataframe_helper():
    args = make_rna_reads_arg_parser().parse_args(args_with_bam)
    df = variants_reads_dataframe_from_args(args)
    assert set(["prefix", "allele", "suffix", "name", "sequence", "gene"]).issubset(df.columns)
    assert len(df) == 42


def _cli_error_message(fn, args, capsys):
    import pytest
    with pytest.raises(SystemExit) as exit_info:
        fn(args)
    assert exit_info.value.code == 2
    return capsys.readouterr().err.strip().splitlines()[-1]


def test_cli_missing_vcf_is_a_one_line_error(capsys):
    message = _cli_error_message(
        isovar_main,
        ["--vcf", "missing.vcf", "--bam", args_with_bam[-1]],
        capsys)
    assert message.endswith("error: --vcf file not found: missing.vcf")


def test_cli_missing_bam_is_a_one_line_error(capsys):
    message = _cli_error_message(
        isovar_allele_counts, vcf_args + ["--bam", "missing.bam"], capsys)
    assert message.endswith("error: --bam file not found: missing.bam")


def test_cli_unindexed_bam_is_a_one_line_error(capsys):
    unindexed_bam = data_path("data/primary.chr1.unsorted.bam")
    message = _cli_error_message(
        isovar_variant_reads, vcf_args + ["--bam", unindexed_bam], capsys)
    assert "must be a coordinate-sorted, indexed BAM or CRAM" in message


def test_cli_sam_passed_as_bam_is_a_one_line_error(capsys):
    sam_path = data_path("data/b16.f10/b16.combined.sam")
    message = _cli_error_message(
        isovar_allele_reads, vcf_args + ["--bam", sam_path], capsys)
    assert "must be a coordinate-sorted, indexed BAM or CRAM" in message


def test_cli_no_variants_is_a_one_line_error(capsys):
    message = _cli_error_message(
        isovar_reference_contexts, [], capsys)
    assert message.endswith(
        "error: no variants given; use --vcf, --maf, --variant or --json-variants")


def test_cli_variant_without_genome_is_a_one_line_error(capsys):
    message = _cli_error_message(
        isovar_reference_contexts, ["--variant", "9", "82927102", "G", "T"], capsys)
    assert message.endswith("error: --genome is required when using --variant")


def test_cli_unknown_genome_is_a_one_line_error(capsys):
    message = _cli_error_message(
        isovar_reference_contexts, vcf_args + ["--genome", "not-a-genome"], capsys)
    assert "error: --genome not-a-genome:" in message


def test_cli_missing_output_directory_fails_before_running(capsys, tmp_path):
    output_path = str(tmp_path / "missing-dir" / "out.csv")
    message = _cli_error_message(
        isovar_main, args_with_bam + ["--output", output_path], capsys)
    assert message.endswith(
        "error: --output directory does not exist: %s" % (tmp_path / "missing-dir"))
