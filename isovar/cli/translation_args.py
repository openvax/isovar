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

from ..default_parameters import (
    MIN_TRANSCRIPT_PREFIX_LENGTH,
    MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
    COUNT_MISMATCHES_AFTER_VARIANT,
    PROTEIN_SEQUENCE_PREFERENCE,
    PROTEIN_CONTEXT_PEPTIDE_LENGTH,
    MIN_PROTEIN_SEQUENCE_SUPPORT_FRACTION,
)
from .variant_sequences_args import make_variant_sequences_arg_parser
from .reference_context_args import add_reference_context_args

def add_translation_args(parser):
    translation_group = parser.add_argument_group(
        "Translation from cDNA to protein sequence")

    translation_group.add_argument(
        "--protein-sequence-length",
        default=None,
        type=int,
        help="Requested translated context length; default is 2 * protein-context-peptide-length - 1 (49).")

    translation_group.add_argument(
        "--protein-sequence-preference", choices=("balanced", "support", "context"),
        default=PROTEIN_SEQUENCE_PREFERENCE,
        help="Balance vaccine context within an RNA support budget, prioritize support, or prioritize context.")
    translation_group.add_argument(
        "--protein-context-peptide-length", type=int, default=None,
        help="Peptide size for counting mutation-containing windows; default 25, or vaxrank's vaccine peptide size.")
    translation_group.add_argument(
        "--min-protein-sequence-support-fraction", type=float,
        default=MIN_PROTEIN_SEQUENCE_SUPPORT_FRACTION,
        help="Balanced preference: minimum fraction of best candidate's read-name support (default 0.9, not confidence).")

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


def protein_sequence_creator_kwargs_from_args(args):
    """One adapter for translation/protein CLIs and vaxrank's namespace."""
    peptide_length = getattr(args, "protein_context_peptide_length", None)
    if peptide_length is None:
        peptide_length = getattr(args, "vaccine_peptide_length", PROTEIN_CONTEXT_PEPTIDE_LENGTH)
    if peptide_length is None:
        peptide_length = PROTEIN_CONTEXT_PEPTIDE_LENGTH
    return dict(
        protein_sequence_length=args.protein_sequence_length,
        min_variant_sequence_coverage=args.min_variant_sequence_coverage,
        min_transcript_prefix_length=args.min_transcript_prefix_length,
        max_transcript_mismatches=args.max_reference_transcript_mismatches,
        count_mismatches_after_variant=getattr(args, "count_mismatches_after_variant", COUNT_MISMATCHES_AFTER_VARIANT),
        variant_sequence_assembly=args.variant_sequence_assembly,
        protein_sequence_preference=getattr(args, "protein_sequence_preference", PROTEIN_SEQUENCE_PREFERENCE),
        protein_context_peptide_length=peptide_length,
        min_protein_sequence_support_fraction=getattr(
            args, "min_protein_sequence_support_fraction", MIN_PROTEIN_SEQUENCE_SUPPORT_FRACTION))


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
