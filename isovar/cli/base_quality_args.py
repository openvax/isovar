"""Explicit, shared reported-quality assessment for RNA evidence exports."""

from ..base_quality import BaseQualityPolicy
from .validation import CommandInputError, non_negative_int


def add_base_quality_args(parser):
    parser.add_argument(
        "--min-base-quality", type=non_negative_int,
        help="Also report focal-allele support meeting this reported Phred threshold (0-254). "
             "Missing quality is unassessed; reconstruction and raw counts are unchanged. "
             "Empty alleles assess flanking anchors only. Protein/cell assessments are in JSON.")
    return parser


def base_quality_policy_from_args(args):
    threshold = getattr(args, "min_base_quality", None)
    try:
        return None if threshold is None else BaseQualityPolicy(threshold)
    except ValueError as error:
        raise CommandInputError(str(error)) from error
