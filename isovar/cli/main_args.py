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

from varcode import Variant

from ..default_parameters import MAX_P_VALUE_FOR_PHASING, PHASING_ERROR_RATE
from ..main import run_isovar
from ..phasing import validate_locus_error_rates, validate_phasing_rates
from .filter_args import add_filter_args, filter_threshold_dict_from_args
from .protein_sequence_args import (
    make_protein_sequences_arg_parser,
    protein_sequence_creator_from_args,
)
from .rna_args import alignment_file_from_args, read_collector_from_args
from .validation import (
    CommandInputError,
    germline_variants_from_args,
    open_fraction,
    variant_collection_from_args,
)


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
    group.add_argument("--phasing-error-rate", type=open_fraction, default=PHASING_ERROR_RATE,
                       help="Wrong-allele probability floor; use a calibrated value for the reads "
                            "(default: %(default)s)")
    group.add_argument("--max-p-value-for-phasing", type=open_fraction, default=MAX_P_VALUE_FOR_PHASING,
                       help="Maximum one-sided p-value for allele support beyond read errors "
                            "(default: %(default)s)")
    group.add_argument("--phasing-error-rates", metavar="JSON",
                       help="Per-locus calibration, as written by isovar calibrate-phasing: JSON rows "
                            "with contig, start, ref, alt, ref_to_alt and alt_to_ref, one per input "
                            "variant they apply to (VCF-style or normalized alleles).")
    return parser


def _locus(variant):
    return (variant.contig, variant.start, variant.ref, variant.alt)


def _phasing_error_rates_from_args(args, variants):
    path = getattr(args, "phasing_error_rates", None)
    if path is None:
        return None
    # Rows may spell a locus as VCF does ("chr7", anchored alleles) or as
    # calibrate-phasing and Varcode do; normalize each the way the inputs were.
    genome = next(iter(variants)).genome if len(variants) else None
    by_key = {_locus(v): v for v in variants}
    try:
        with open(path) as handle:
            rows = json.load(handle)
    except (OSError, ValueError) as error:
        raise CommandInputError("Cannot read --phasing-error-rates %s: %s" % (path, error)) from None
    if not isinstance(rows, list):
        raise CommandInputError("--phasing-error-rates must be a JSON list of locus records")
    rates = {}
    for number, row in enumerate(rows, 1):
        try:
            locus = _locus(Variant(str(row["contig"]), int(row["start"]), row["ref"], row["alt"],
                                   genome, convert_ucsc_contig_names=True))
            pair = (row["ref_to_alt"], row["alt_to_ref"])
        except (KeyError, TypeError, ValueError) as error:
            raise CommandInputError(
                "--phasing-error-rates row %d needs contig, start, ref, alt, ref_to_alt and "
                "alt_to_ref (%s)" % (number, error)) from None
        variant = by_key.get(locus)
        if variant is None:
            raise CommandInputError("--phasing-error-rates row %d (%s %s %s>%s) matches no input variant" % (
                number, row["contig"], row["start"], row["ref"], row["alt"]))
        if variant in rates:
            raise CommandInputError("--phasing-error-rates row %d repeats %s" % (number, variant))
        rates[variant] = pair
    try:
        return validate_locus_error_rates(rates)
    except ValueError as error:
        raise CommandInputError("--phasing-error-rates: %s" % error) from None


def phasing_parameters_from_args(args, variants):
    """Validated phasing keywords for Isovar and downstream CLI callers."""
    rate = getattr(args, "phasing_error_rate", PHASING_ERROR_RATE)
    alpha = getattr(args, "max_p_value_for_phasing", MAX_P_VALUE_FOR_PHASING)
    validate_phasing_rates(rate, alpha)
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
