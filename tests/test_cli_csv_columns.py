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
    assert {"num_translations", "num_supporting_reads", "num_supporting_fragments"} <= set(proteins.columns)
    assert "translations" not in proteins.columns
    sequences = pd.read_csv(run(tmp_path, "variant-sequences"))
    assert {"num_reads", "num_fragments"} <= set(sequences.columns) and "reads" not in sequences.columns
    # This command does not classify transcripts; when present, IDs are listed, not counted.
    assert sequences["compatible_transcript_ids"].isna().all()
    reads = [AlleleRead("AC", "G", "TT", "read%d" % i, compatible_transcript_ids={"T2", "T1"}) for i in range(2)]
    sequence = VariantSequence(prefix="AC", alt="G", suffix="TT", reads=reads)
    frame = variant_sequences_generator_to_dataframe([(Variant("1", 3, "A", "G"), [sequence])])
    assert frame["compatible_transcript_ids"].tolist() == ["T1;T2"]
    assert frame["num_reads"].tolist() == frame["num_fragments"].tolist() == [2]
    translations = pd.read_csv(run(tmp_path, "translations"))
    assert {"in_frame_cdna_sequence", "reference_transcript_names"} <= set(translations.columns)
    assert not {"variant_orf", "reference_context"} & set(translations.columns)


def test_floats_are_written_with_bounded_precision(tmp_path):
    text = run(tmp_path, "run").read_text()
    assert not re.search(r"\d\.\d{7,}", text)


def test_empty_results_keep_their_header(tmp_path):
    header_only = tmp_path / "header-only.vcf"
    header_only.write_text("".join(line for line in open(VCF) if line.startswith("#")))
    df = pd.read_csv(run(tmp_path, "run", vcf=str(header_only)))
    assert df.empty and list(df.columns) == IsovarResult.record_columns()
    assert list(isovar_results_to_dataframe([]).columns) == IsovarResult.record_columns()
