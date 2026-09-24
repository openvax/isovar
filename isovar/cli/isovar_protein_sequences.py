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

"""
Translate variant cDNA assembled from RNA reads into candidate mutant protein
sequences, combining translations with identical amino acids.
"""

from .commands import run_dataframe_command
from .output_args import add_output_args
from .protein_sequence_args import (
    make_protein_sequences_arg_parser,
    protein_sequences_dataframe_from_args,
)

parser = add_output_args(
    make_protein_sequences_arg_parser(description=__doc__),
    filename="isovar-protein-sequences-result.csv",
    description="CSV of candidate protein sequences")


def run(args=None, *, prog=None):
    run_dataframe_command(parser, protein_sequences_dataframe_from_args, args, prog)
