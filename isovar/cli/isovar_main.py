# Copyright (c) 2019. Mount Sinai School of Medicine
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

"""
Primary Isovar command, used to collect information about variants,
the RNA reads which overlap and protein sequences which can be constructed
from reads that support the variant.
"""
from __future__ import print_function, division, absolute_import

import sys

from varcode.cli import variant_collection_from_args

from ..logging import get_logger
from ..main import run_isovar
from ..dataframe_helpers import isovar_results_to_dataframe

from .protein_sequence_args import (
    make_protein_sequences_arg_parser,
    protein_sequence_creator_from_args
)
from .output_args import add_output_args, write_dataframe
from .filter_args import add_filter_args
from .rna_args import read_collector_from_args, alignment_file_from_args

logger = get_logger(__name__)

parser = make_protein_sequences_arg_parser()
parser = add_output_args(
    parser,
    filename="isovar-result.csv")
add_filter_args(parser)

def run(args=None):
    if args is None:
        args = sys.argv[1:]
    args = parser.parse_args(args)
    logger.info(args)
    variants = variant_collection_from_args(args)
    read_collector = read_collector_from_args(args)
    alignment_file = alignment_file_from_args(args)
    protein_sequences_creator = protein_sequence_creator_from_args(args)
    isovar_results = run_isovar(
        variants=variants,
        alignment_file=alignment_file,
        read_collector=read_collector,
        protein_sequences_creator=protein_sequences_creator,
        filter_thresholds={})
    df = isovar_results_to_dataframe(isovar_results)
    logger.info(df)
    write_dataframe(df, args)

