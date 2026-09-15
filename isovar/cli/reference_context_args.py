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

from argparse import ArgumentTypeError

from varcode.cli.variant_args import (
    make_variants_parser,
    variant_collection_from_args
)

from ..reference_context_helpers import reference_contexts_generator
from ..dataframe_helpers import variants_to_reference_contexts_dataframe
from ..default_parameters import REFERENCE_CONTEXT_SIZE


def _positive_context_size(value):
    value = int(value)
    if value < 1:
        raise ArgumentTypeError("reference context size must be a positive integer")
    return value


def add_reference_context_args(parser):
    """
    Extends an ArgumentParser instance with the following commandline arguments:
        --reference-context-size
    """
    reference_context_group = parser.add_argument_group("Reference Transcripts")
    reference_context_group.add_argument(
        "--reference-context-size",
        type=_positive_context_size,
        default=REFERENCE_CONTEXT_SIZE,
        help=(
            "Maximum reference nucleotides on each side of the variant "
            "in the exported contexts (default %(default)s)."))
    return reference_context_group


def make_reference_context_arg_parser(**kwargs):
    """
    Parameters
    ----------
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
