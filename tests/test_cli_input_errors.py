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

"""Bad command input is a one-line usage error (exit 2), never a traceback."""

import json

import pytest

from isovar.cli import commands
from isovar.fusion import fusion_from_dict
from isovar.sv_rna import sv_rna_input_from_dict

from .testing_helpers import data_path

VCF = data_path("data/b16.f10/b16.vcf")
BAM = data_path("data/b16.f10/b16.combined.sorted.bam")


def usage_error(capsys, args):
    with pytest.raises(SystemExit) as error:
        commands.run(args)
    assert error.value.code == 2
    message = capsys.readouterr().err.strip().splitlines()[-1]
    assert "Traceback" not in message
    return message


@pytest.mark.parametrize("args,expected", [
    (["run", "--vcf", "nope.vcf", "--bam", BAM], "No such file"),
    (["run", "--vcf", VCF, "--bam", "nope.bam"], "Cannot open alignment file nope.bam"),
    (["run", "--bam", BAM], "No variants loaded"),
    (["reference-contexts", "--json-variants", "nope.json"], "No such file"),
    (["allele-counts", "--variant", "chr9", "82927102", "G", "T", "--bam", BAM], "--genome"),
    (["run", "--vcf", VCF, "--bam", BAM, "--genome", "foo"], "specify the reference with --genome"),
    (["run", "--vcf", VCF, "--bam", BAM, "--genome", "GRCh38:75"], "release 75 of homo_sapiens provides GRCh37"),
    (["run", "--vcf", VCF, "--bam", data_path("data/primary.chr1.unsorted.bam")], "has no index"),
    (["run", "--vcf", VCF, "--bam", data_path("data/b16.f10/b16.combined.sam")], "is not BAM or CRAM"),
    (["run", "--vcf", VCF, "--bam", VCF], "does not contain alignment data"),
    (["run", "--vcf", VCF, "--bam", BAM, "--trim-adapters"], "read_end_profile"),
    (["variant-sequences", "--vcf", VCF, "--bam", BAM, "--variant-sequence-length", "-5"], "positive integer"),
    (["protein-sequences", "--vcf", VCF, "--bam", BAM, "--protein-sequence-length", "0"], "positive integer"),
    (["run", "--vcf", VCF, "--bam", BAM, "--min-alt-rna-fraction", "2.0"], "between 0 and 1"),
    (["run", "--vcf", VCF, "--bam", BAM, "--min-protein-sequence-support-fraction", "2"], "between 0 and 1"),
    (["reference-contexts", "--vcf", VCF, "--reference-context-size", "-7"], "positive integer"),
    (["fusion", "--input", "x.json", "--output", "y.json", "--dpi", "10"], "at least 72 dpi"),
])
def test_bad_inputs_are_usage_errors(capsys, tmp_path, args, expected):
    if "--output" not in args:
        args = args + ["--output", str(tmp_path / "out.csv")]
    assert expected in usage_error(capsys, args)


def test_genome_can_choose_an_ensembl_release():
    # Previously "GRCh38:93" silently gave the most recent release (#122).
    from isovar.cli.isovar_allele_counts import parser
    from isovar.cli.validation import variant_collection_from_args

    args = parser.parse_args(["--variant", "chr9", "82927102", "G", "T", "--bam", BAM,
                              "--genome", "GRCh38:93"])
    assert {v.genome.release for v in variant_collection_from_args(args)} == {93}


def test_output_directory_is_checked_before_reading_any_input(capsys, tmp_path):
    message = usage_error(capsys, ["run", "--vcf", "nope.vcf", "--bam", "nope.bam",
                                   "--output", str(tmp_path / "missing" / "out.csv")])
    assert "Output directory does not exist" in message


def test_unknown_output_column_is_a_usage_error(capsys, tmp_path):
    message = usage_error(capsys, ["allele-counts", "--vcf", VCF, "--bam", BAM, "--output",
                                   str(tmp_path / "out.csv"), "--output-columns", "chr", "nope"])
    assert "Unknown --output-columns nope" in message and "num_alt_reads" in message


@pytest.mark.parametrize("decode,data,expected", [
    (fusion_from_dict, {}, "missing 'fusion'"),
    (fusion_from_dict, {"fusion": {"donor": {"unexpected": 1}}}, "Invalid fusion input"),
    (sv_rna_input_from_dict, {}, "missing 'event_id'"),
])
def test_input_decoders_report_malformed_json_as_value_errors(decode, data, expected):
    with pytest.raises(ValueError, match=expected):
        decode(data)


@pytest.mark.parametrize("command", ["fusion", "sv-rna"])
def test_json_commands_report_malformed_input_as_usage_errors(capsys, tmp_path, command):
    event = tmp_path / "input.json"
    event.write_text(json.dumps({}))
    args = [command, "--input", str(event), "--output", str(tmp_path / "out.json")]
    if command == "sv-rna":
        args += ["--bam", BAM]
    assert "missing" in usage_error(capsys, args)
