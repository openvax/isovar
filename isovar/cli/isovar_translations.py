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
Translate each non-synonymous coding variants into possible mutant protein
sequences using an RNAseq BAM from the same tissuie.
"""

from __future__ import print_function, division, absolute_import
import sys

from ..logging import get_logger
from ..protein_sequence_creator import ProteinSequenceCreator
from ..dataframe_helpers import translations_generator_to_dataframe

from .translation_args import make_translation_arg_parser
from .rna_args import read_evidence_generator_from_args
from .output_args import add_output_args, write_dataframe

logger = get_logger(__name__)

parser = make_translation_arg_parser()
parser = add_output_args(
    parser,
    filename="isovar-translate-variants-results.csv",
    description="Name of CSV file which contains predicted sequences")


def translations_generator_from_args(args):
    """
    Given parsed commandline arguments, returns a generator whose elements
    are (varcode.Variant, [Translation])
    """
    read_evidence_generator = read_evidence_generator_from_args(args)
    protein_sequence_creator = ProteinSequenceCreator(
        protein_sequence_length=args.protein_sequence_length,
        min_variant_sequence_coverage=args.min_variant_sequence_coverage,
        variant_sequence_assembly=args.variant_sequence_assembly,
        min_transcript_prefix_length=args.min_transcript_prefix_length,
        max_transcript_mismatches=args.max_reference_transcript_mismatches,
        count_mismatches_after_variant=args.count_mismatches_after_variant)
    return protein_sequence_creator.translate_variants(read_evidence_generator)


def translations_dataframe_from_args(args):
    """
    Collects Translation objects based on commandline arguments and
    converts them into a DataFrame.
    """
    translations_generator = translations_generator_from_args(args)
    return translations_generator_to_dataframe(translations_generator)


def run(args=None):
    if args is None:
        args = sys.argv[1:]
    args = parser.parse_args(args)
    logger.info(args)
    df = translations_dataframe_from_args(args)
    logger.info(df)
    write_dataframe(df, args)
