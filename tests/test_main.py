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

from isovar import ProteinSequenceCreator, run_isovar, isovar_results_to_dataframe
from isovar.cli.variant_sequences_args import make_variant_sequences_arg_parser
from isovar.variant_sequence_creator import VariantSequenceCreator
from isovar.cli.main_args import make_isovar_arg_parser, run_isovar_from_parsed_args
from pandas.testing import assert_frame_equal
from .common import eq_
from .testing_helpers import data_path


def test_isovar_main_to_dataframe():
    results = run_isovar(
        variants=data_path("data/b16.f10/b16.vcf"),
        alignment_file=data_path("data/b16.f10/b16.combined.sorted.bam"))
    df = isovar_results_to_dataframe(results)
    print(df)
    eq_(len(df), 4)
    # B16 test data has 2/4 variants with enough coverage
    # to translate protein sequences
    eq_(df["passes_all_filters"].sum(), 2)
    cli_args = make_isovar_arg_parser().parse_args([
        "--vcf", data_path("data/b16.f10/b16.vcf"),
        "--bam", data_path("data/b16.f10/b16.combined.sorted.bam"),
    ])
    cli_df = isovar_results_to_dataframe(run_isovar_from_parsed_args(cli_args))
    assert_frame_equal(cli_df, df)


def test_python_and_cli_enable_variant_sequence_assembly_by_default():
    protein_creator = ProteinSequenceCreator()
    variant_creator = VariantSequenceCreator()
    cli_args = make_variant_sequences_arg_parser().parse_args(["--bam", "unused.bam"])

    assert protein_creator.variant_sequence_assembly is True
    assert variant_creator.variant_sequence_assembly is True
    assert cli_args.variant_sequence_assembly is True


def test_cli_can_disable_variant_sequence_assembly():
    cli_args = make_variant_sequences_arg_parser().parse_args([
        "--bam", "unused.bam",
        "--disable-variant-sequence-assembly",
    ])

    assert cli_args.variant_sequence_assembly is False
