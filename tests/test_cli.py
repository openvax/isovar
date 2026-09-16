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
import pytest

from .testing_helpers import data_path

from isovar.cli.isovar_translations import run as isovar_translations
from isovar.cli.isovar_allele_counts import run as isovar_allele_counts
from isovar.cli.isovar_allele_reads import run as isovar_allele_reads
from isovar.cli.isovar_protein_sequences import run as isovar_protein_sequences
from isovar.cli.isovar_reference_contexts import run as isovar_reference_contexts
from isovar.cli.isovar_variant_reads import run as isovar_variant_reads
from isovar.cli.isovar_variant_sequences import run as isovar_variant_sequences
from isovar.cli.isovar_main import run as isovar_main
from isovar.cli.main_args import make_isovar_arg_parser
from isovar.cli.rna_args import (
    allele_counts_dataframe_from_args,
    make_rna_reads_arg_parser,
    variants_reads_dataframe_from_args,
)

# the modules themselves, for the parsers they build at import time
from isovar.cli import isovar_allele_counts as isovar_allele_counts_module
from isovar.cli import isovar_allele_reads as isovar_allele_reads_module
from isovar.cli import isovar_protein_sequences as isovar_protein_sequences_module
from isovar.cli import isovar_reference_contexts as isovar_reference_contexts_module
from isovar.cli import isovar_translations as isovar_translations_module
from isovar.cli import isovar_variant_reads as isovar_variant_reads_module
from isovar.cli import isovar_variant_sequences as isovar_variant_sequences_module

vcf_args = [
    "--vcf",
    data_path("data/b16.f10/b16.vcf")
]

args_with_bam = vcf_args + [
    "--bam",
    data_path("data/b16.f10/b16.combined.sorted.bam")
]


def run_cli_fn(fn, include_bam_in_args=True, return_dataframe=False, extra_args=()):
    with tempfile.NamedTemporaryFile(delete=False) as f:
        output_path = f.name
    assert not exists(output_path) == 0
    output_args = ["--output", output_path]
    if include_bam_in_args:
        args = args_with_bam + output_args
    else:
        args = vcf_args + output_args
    fn(args + list(extra_args))
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


@pytest.mark.parametrize("extra_args,expected_count", [
    ([], 274), (["--no-merge-overlapping-fragments"], 294),
])
def test_cli_allele_reads(extra_args, expected_count):
    df = run_cli_fn(isovar_allele_reads, return_dataframe=True, extra_args=extra_args)
    assert set(["prefix", "allele", "suffix", "name", "sequence", "gene"]).issubset(df.columns)
    assert "ref_reads" not in df.columns
    assert len(df) == expected_count


def test_cli_reference_contexts():
    run_cli_fn(isovar_reference_contexts, include_bam_in_args=False)


def test_cli_protein_sequences():
    run_cli_fn(isovar_protein_sequences)


def test_cli_translations():
    run_cli_fn(isovar_translations)


@pytest.mark.parametrize("extra_args,expected_count", [
    ([], 36), (["--no-merge-overlapping-fragments"], 42),
])
def test_cli_variant_reads(extra_args, expected_count):
    df = run_cli_fn(isovar_variant_reads, return_dataframe=True, extra_args=extra_args)
    assert set(["prefix", "allele", "suffix", "name", "sequence", "gene"]).issubset(df.columns)
    assert "ref_reads" not in df.columns
    assert len(df) == expected_count


def test_cli_variant_sequences():
    run_cli_fn(isovar_variant_sequences)

def test_cli_main():
    run_cli_fn(isovar_main)


@pytest.mark.parametrize("extra_args,expected_count", [
    ([], 36), (["--no-merge-overlapping-fragments"], 42),
])
def test_variant_reads_dataframe_helper(extra_args, expected_count):
    args = make_rna_reads_arg_parser().parse_args(args_with_bam + extra_args)
    df = variants_reads_dataframe_from_args(args)
    assert set(["prefix", "allele", "suffix", "name", "sequence", "gene"]).issubset(df.columns)
    assert len(df) == expected_count


def test_cli_help_describes_every_command_and_option(capsys):
    # the eight commands build their parsers at import time, except isovar-main
    # which builds its own inside run(), so check that one separately
    described_parsers = [
        module.parser
        for module in (
            isovar_allele_counts_module,
            isovar_allele_reads_module,
            isovar_protein_sequences_module,
            isovar_reference_contexts_module,
            isovar_translations_module,
            isovar_variant_reads_module,
            isovar_variant_sequences_module,
        )
    ]
    for parser in described_parsers:
        assert parser.description, parser.prog

    # argparse only exposes its actions through the private _actions list
    for parser in described_parsers + [make_isovar_arg_parser()]:
        for action in parser._actions:
            assert action.help, (parser.prog, action.option_strings)

    with pytest.raises(SystemExit) as exit_info:
        isovar_main(["--help"])
    assert exit_info.value.code == 0
    # argparse wraps the description to the terminal width, so collapse
    # whitespace before looking for the phrase
    help_text = " ".join(capsys.readouterr().out.split())
    assert "Collect RNA evidence for each variant" in help_text
