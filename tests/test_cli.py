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
    run_cli_fn(isovar_allele_counts)


def test_cli_allele_reads():
    df = run_cli_fn(isovar_allele_reads, return_dataframe=True)
    assert set(["prefix", "allele", "suffix", "name", "sequence", "gene"]).issubset(df.columns)
    assert "ref_reads" not in df.columns
    assert len(df) == 293


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
