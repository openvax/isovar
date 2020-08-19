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
Common command-line arguments for all Isovar commands which translate
cDNA into protein sequences.
"""

from __future__ import print_function, division, absolute_import

from ..default_parameters import (
    MIN_TRANSCRIPT_PREFIX_LENGTH,
    MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
    COUNT_MISMATCHES_AFTER_VARIANT,
    PROTEIN_SEQUENCE_LENGTH,
)
from .variant_sequences_args import make_variant_sequences_arg_parser
from .reference_context_args import add_reference_context_args

def add_translation_args(parser):
    translation_group = parser.add_argument_group(
        "Translation from cDNA to protein sequence")

    translation_group.add_argument(
        "--protein-sequence-length",
        default=PROTEIN_SEQUENCE_LENGTH,
        type=int)

    translation_group.add_argument(
        "--max-reference-transcript-mismatches",
        type=int,
        default=MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
        help=(
            "Maximum number of mismatches between variant sequence"
            " reference sequence before a candidate reading frame is ignored."))

    translation_group.add_argument(
        "--count-mismatches-after-variant",
        action="store_true",
        default=COUNT_MISMATCHES_AFTER_VARIANT,
        help="If true, mismatches after the variant locus will count toward the "
             "--max-reference-transcript-mismatches filter.")

    translation_group.add_argument(
        "--min-transcript-prefix-length",
        type=int,
        default=MIN_TRANSCRIPT_PREFIX_LENGTH,
        help=(
            "Number of nucleotides before the variant we try to match against "
            "a reference transcript. Values greater than zero exclude variants "
            "near the start codon of transcripts without 5' UTRs."))

    return translation_group


def make_translation_arg_parser(**kwargs):
    """
    Parameters
    ----------
    **kwargs : dict
        Passed directly to argparse.ArgumentParser

    Creates argparse.ArgumentParser instance with all of the options
    needed to translate each distinct cDNA sequence determined from
    variants & RNAseq.

    See `args.variant_sequences` for commandline parameters which aren't added
    in this module.
    """
    parser = make_variant_sequences_arg_parser(**kwargs)
    add_reference_context_args(parser)
    add_translation_args(parser)
    return parser
