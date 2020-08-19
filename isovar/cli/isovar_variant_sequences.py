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

from __future__ import print_function, division, absolute_import

import sys

from ..logging import get_logger

from ..variant_sequence_creator import VariantSequenceCreator
from ..dataframe_helpers import variant_sequences_generator_to_dataframe


from .rna_args import read_evidence_generator_from_args
from .variant_sequences_args import make_variant_sequences_arg_parser
from .output_args import add_output_args, write_dataframe

logger = get_logger(__name__)

parser = make_variant_sequences_arg_parser(add_sequence_length_arg=True)
parser = add_output_args(
    parser,
    filename="isovar-variant-sequences-results.csv",
    description="Name of CSV file which contains predicted sequences")


def variant_sequences_generator_from_args(args):
    """
    Use parsed commandline arguments to load variants and RNA reads and
    generate a sequence of (Variant, list of VariantSequence) pairs.
    """
    read_evidence_generator = read_evidence_generator_from_args(args)
    variant_sequence_creator = VariantSequenceCreator(
        min_variant_sequence_coverage=args.min_variant_sequence_coverage,
        preferred_sequence_length=args.variant_sequence_length,
        variant_sequence_assembly=args.variant_sequence_assembly)
    return variant_sequence_creator.sequences_from_read_evidence_generator(read_evidence_generator)

def variant_sequences_dataframe_from_args(args):
    variant_sequences_generator = variant_sequences_generator_from_args(args)
    return variant_sequences_generator_to_dataframe(variant_sequences_generator)


def run(args=None):
    if args is None:
        args = sys.argv[1:]
    args = parser.parse_args(args)
    logger.info(args)
    df = variant_sequences_dataframe_from_args(args)
    logger.info(df)
    write_dataframe(df, args)
