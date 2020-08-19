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

from varcode.cli.variant_args import (
    make_variants_parser,
    variant_collection_from_args
)

from ..reference_context_helpers import reference_contexts_generator
from ..dataframe_helpers import variants_to_reference_contexts_dataframe


def add_reference_context_args(parser):
    """
    Extends an ArgumentParser instance with the following commandline arguments:
        --context-size
    """
    reference_context_group = parser.add_argument_group("Reference Transcripts")
    reference_context_group.add_argument(
        "--reference-context-size",
        type=int,
        default=30,
        help=(
            "Number of nucleotides used to match assembled sequence to "
            "reference transcript to establish reading frame."))
    return reference_context_group


def make_reference_context_arg_parser(**kwargs):
    """
    Parameters
    ----------
    add_context_size_arg : bool
        If True then add a `--context-size` argument, which is otherwise
        inferred from cDNA sequence length.
    **kwargs : dict
        Parameters passed directly to argparse.ArgumentParser.

    Returns an argparse.ArgumentParser instance.
    """
    parser = make_variants_parser(**kwargs)
    add_reference_context_args(parser)
    return parser


def reference_contexts_dataframe_from_args(args):
    """
    Generate a DataFrame for variants and their associated reference contexts
    loaded based on parsed commandline arguments.
    """
    variants = variant_collection_from_args(args)
    reference_context_gen = reference_contexts_generator(
        variants=variants,
        context_size=args.reference_context_size)
    return variants_to_reference_contexts_dataframe(reference_context_gen)
