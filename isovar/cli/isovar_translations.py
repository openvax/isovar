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
Export every translation of each assembled variant cDNA sequence in the reading
frame of each compatible reference transcript, before candidates are grouped
and ranked into protein sequences.
"""

from ..dataframe_helpers import translations_generator_to_dataframe
from ..protein_sequence_creator import ProteinSequenceCreator
from .commands import run_dataframe_command
from .output_args import add_output_args
from .rna_args import read_evidence_generator_from_args
from .translation_args import make_translation_arg_parser, protein_sequence_creator_kwargs_from_args

parser = add_output_args(
    make_translation_arg_parser(description=__doc__),
    filename="isovar-translate-variants-results.csv",
    description="CSV of translations")


def translations_generator_from_args(args):
    """
    Given parsed commandline arguments, returns a generator whose elements
    are (varcode.Variant, [Translation])
    """
    read_evidence_generator = read_evidence_generator_from_args(args)
    protein_sequence_creator = ProteinSequenceCreator(**protein_sequence_creator_kwargs_from_args(args))
    return protein_sequence_creator.translate_variants(read_evidence_generator)


def translations_dataframe_from_args(args):
    """
    Collects Translation objects based on commandline arguments and
    converts them into a DataFrame.
    """
    return translations_generator_to_dataframe(translations_generator_from_args(args))


def run(args=None, *, prog=None):
    run_dataframe_command(parser, translations_dataframe_from_args, args, prog)
