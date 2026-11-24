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

from ..default_parameters import REFERENCE_CONTEXT_SIZE
from .help import HelpArgumentParser
from .validation import positive_int, variant_collection_from_args

# Lazily resolved on first use; kept as module attributes for test monkeypatching.
reference_contexts_generator = None
variants_to_reference_contexts_dataframe = None


def add_reference_context_args(parser):
    """
    Extends an ArgumentParser instance with the following commandline arguments:
        --reference-context-size
    """
    reference_context_group = parser.add_argument_group("Reference Transcripts")
    reference_context_group.add_argument(
        "--reference-context-size",
        type=positive_int,
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
        Parameters passed directly to HelpArgumentParser.

    Returns an HelpArgumentParser instance.
    """
    from varcode.cli.variant_args import add_variant_args

    parser = HelpArgumentParser(**kwargs)
    add_variant_args(parser)
    add_reference_context_args(parser)
    return parser


def reference_contexts_dataframe_from_args(args):
    """
    Generate a DataFrame for variants and their associated reference contexts
    loaded based on parsed commandline arguments.
    """
    generator_fn = reference_contexts_generator
    if generator_fn is None:
        from ..reference_context_helpers import reference_contexts_generator as generator_fn
    to_dataframe_fn = variants_to_reference_contexts_dataframe
    if to_dataframe_fn is None:
        from ..dataframe_helpers import variants_to_reference_contexts_dataframe as to_dataframe_fn

    variants = variant_collection_from_args(args)
    reference_context_gen = generator_fn(
        variants=variants,
        context_size=args.reference_context_size)
    return to_dataframe_fn(reference_context_gen)
