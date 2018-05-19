import tempfile
from os import remove
from os.path import getsize, exists

from testing_helpers import data_path

from isovar.cli.isovar_translations import run as isovar_translations
from isovar.cli.isovar_allele_counts import run as isovar_allele_counts
from isovar.cli.isovar_allele_reads import run as isovar_allele_reads
from isovar.cli.isovar_protein_sequences import run as isovar_protein_sequences
from isovar.cli.isovar_reference_contexts import run as isovar_reference_contexts
from isovar.cli.isovar_variant_reads import run as isovar_variant_reads
from isovar.cli.isovar_variant_sequences import run as isovar_variant_sequences


vcf_args = [
    "--vcf",
    data_path("data/b16.f10/b16.vcf")
]

args_with_bam = vcf_args + [
    "--bam",
    data_path("data/b16.f10/b16.combined.sorted.bam")
]


def run_cli_fn(fn, include_bam_in_args=True):
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
    remove(output_path)


def test_cli_allele_counts():
    run_cli_fn(isovar_allele_counts)


def test_cli_allele_reads():
    run_cli_fn(isovar_allele_reads)


def test_cli_reference_contexts():
    run_cli_fn(isovar_reference_contexts, include_bam_in_args=False)


def test_cli_protein_sequences():
    run_cli_fn(isovar_protein_sequences)


def test_cli_translations():
    run_cli_fn(isovar_translations)


def test_cli_variant_reads():
    run_cli_fn(isovar_variant_reads)


def test_cli_variant_sequences():
    run_cli_fn(isovar_variant_sequences)
