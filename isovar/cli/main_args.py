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
Create parser and run Isovar from parsed args
"""

import json

from ..main import run_isovar
from ..default_parameters import PHASING_ERROR_RATE, MAX_P_VALUE_FOR_PHASING
from ..phasing import _validate_locus_error_rates, _validate_phasing_rates

from .protein_sequence_args import (
    make_protein_sequences_arg_parser,
    protein_sequence_creator_from_args
)
from .filter_args import add_filter_args, filter_threshold_dict_from_args
from .rna_args import read_collector_from_args, alignment_file_from_args
from .validation import germline_variants_from_args, variant_collection_from_args

def make_isovar_arg_parser(**kwargs):
    """
    Create argument parser with all options required to run Isovar

    Parameters
    ----------
    **kwargs : dict
        Passed directly to argparse.ArgumentParser
    """
    parser = make_protein_sequences_arg_parser(**kwargs)
    add_filter_args(parser)
    parser.add_argument_group("Matched germline").add_argument(
        "--germline-vcf",
        help="Variants from a matched normal. Assembled cDNA edits that are one of "
             "these are reported as known germline; edits that are other input "
             "variants are co-somatic; the rest stay unexplained.")
    group = parser.add_argument_group("Read phasing")
    group.add_argument("--phasing-error-rate", type=float, default=PHASING_ERROR_RATE,
                       help="Wrong-allele probability floor; use a calibrated value for the reads.")
    group.add_argument("--max-p-value-for-phasing", type=float, default=MAX_P_VALUE_FOR_PHASING,
                       help="Maximum one-sided p-value for allele support beyond read errors.")
    group.add_argument("--phasing-error-rates", metavar="JSON",
                       help="Per-locus calibration: JSON rows with contig, start, ref, alt, "
                            "ref_to_alt and alt_to_ref. Coordinates must match the input variants.")
    return parser


def _phasing_error_rates_from_args(args, variants):
    path = getattr(args, "phasing_error_rates", None)
    if path is None:
        return None
    by_key = {(v.contig, v.start, v.ref, v.alt): v for v in variants}
    with open(path) as handle:
        rows = json.load(handle)
    if not isinstance(rows, list):
        raise ValueError("Phasing calibration must be a JSON list of locus records")
    rates = {}
    for row in rows:
        try:
            key = (row["contig"], row["start"], row["ref"], row["alt"])
            pair = (row["ref_to_alt"], row["alt_to_ref"])
            variant = by_key[key]
        except (KeyError, TypeError) as error:
            raise ValueError("Phasing calibration row has missing fields or does not match an input variant") from error
        if variant in rates:
            raise ValueError("Duplicate locus in phasing calibration: %s" % (key,))
        rates[variant] = pair
    return _validate_locus_error_rates(rates)

def phasing_parameters_from_args(args, variants):
    """Validated phasing keywords for Isovar and downstream CLI callers."""
    rate = getattr(args, "phasing_error_rate", PHASING_ERROR_RATE)
    alpha = getattr(args, "max_p_value_for_phasing", MAX_P_VALUE_FOR_PHASING)
    _validate_phasing_rates(rate, alpha)
    return dict(phasing_error_rate=rate, max_p_value_for_phasing=alpha,
                phasing_error_rates=_phasing_error_rates_from_args(args, variants))


def run_isovar_from_parsed_args(args):
    """
    Extract parameters from parsed arguments and use them to run Isovar
    """
    read_collector = read_collector_from_args(args)
    protein_sequence_creator = protein_sequence_creator_from_args(args)
    filter_thresholds = filter_threshold_dict_from_args(args)
    variants = variant_collection_from_args(args)
    germline_variants = germline_variants_from_args(args, variants)
    alignment_file = alignment_file_from_args(args)
    return run_isovar(
        variants=variants,
        alignment_file=alignment_file,
        read_collector=read_collector,
        protein_sequence_creator=protein_sequence_creator,
        filter_thresholds=filter_thresholds,
        germline_variants=germline_variants,
        **phasing_parameters_from_args(args, variants))
