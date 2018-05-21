# Copyright (c) 2016-2018. Mount Sinai School of Medicine
#
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

from ..variant_sequences import (
    reads_generator_to_sequences_generator,
    variant_sequences_generator_to_dataframe
)
from .rna_args import allele_reads_generator_from_args
from .variant_sequences_args import make_variant_sequences_arg_parser
from .output_args import add_output_args, write_dataframe

logger = get_logger(__name__)

parser = make_variant_sequences_arg_parser(add_sequence_length_arg=True)
parser = add_output_args(
    parser,
    filename="isovar-variant-sequences-results.csv",
    description="Name of CSV file which contains predicted sequences")


def variant_sequences_generator_from_args(args):
    allele_reads_generator = allele_reads_generator_from_args(args)
    return reads_generator_to_sequences_generator(
        allele_reads_generator,
        min_alt_rna_reads=args.min_alt_rna_reads,
        min_variant_sequence_coverage=args.min_variant_sequence_coverage,
        preferred_sequence_length=args.variant_sequence_length,
        variant_sequence_assembly=args.variant_sequence_assembly)


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
