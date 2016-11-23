# Copyright (c) 2016. Mount Sinai School of Medicine
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

from ..default_parameters import (
    MIN_ALT_RNA_READS,
    MIN_VARIANT_SEQUENCE_COVERAGE,
    VARIANT_SEQUENCE_LENGTH,
)
from ..variant_sequences import (
    reads_generator_to_sequences_generator,
    variant_sequences_generator_to_dataframe
)
from .rna_reads import (
    allele_reads_generator_from_args,
    make_rna_reads_arg_parser
)
def add_variant_sequence_args(
        parser,
        add_sequence_length_arg=False):
    rna_sequence_group = parser.add_argument_group(
        "Determine coding sequence from RNA")

    rna_sequence_group.add_argument(
        "--min-alt-rna-reads",
        type=int,
        default=MIN_ALT_RNA_READS,
        help="Minimum number of reads supporting variant allele (default %(default)s)")

    rna_sequence_group.add_argument(
        "--min-variant-sequence-coverage",
        type=int,
        default=MIN_VARIANT_SEQUENCE_COVERAGE,
        help="Minimum number of reads supporting a variant sequence (default %(default)s)")

    rna_sequence_group.add_argument(
        "--disable-variant-sequence-assembly",
        dest="variant_sequence_assembly",
        default=True,
        action="store_false",
        help="Disable assemble variant cDNA sequence from overlapping reads")
    # when cDNA sequence length can be inferred from a protein length then
    # we may want to omit this arg
    if add_sequence_length_arg:
        rna_sequence_group.add_argument(
            "--variant-sequence-length",
            default=VARIANT_SEQUENCE_LENGTH,
            type=int)
    return parser

def make_variant_sequences_arg_parser(add_sequence_length_arg=False, **kwargs):
    """
    Parameters
    ----------
    add_sequence_length_arg : bool
        If True then add the `--cdna-sequence-length` argument. This may be
        omitted if the cDNA sequence length is inferred from a protein length.

    **kwargs : dict
        Passed directly to argparse.ArgumentParser

    Creates argparse.ArgumentParser instance with all of the options
    needed to translate each distinct cDNA sequence determined from
    variants & RNAseq.

    See `args.variant_sequences` for commandline parameters which aren't added
    in this module.
    """
    parser = make_rna_reads_arg_parser(**kwargs)
    add_variant_sequence_args(
        parser=parser,
        add_sequence_length_arg=add_sequence_length_arg)
    return parser

def variant_sequences_generator_from_args(args):
    allele_reads_generator = allele_reads_generator_from_args(args)
    return reads_generator_to_sequences_generator(
        allele_reads_generator,
        min_alt_rna_reads=args.min_alt_rna_reads,
        min_variant_sequence_coverage=args.min_variant_sequence_coverage,
        preferred_sequence_length=args.cdna_sequence_length,
        variant_sequence_assembly=args.variant_sequence_assembly)

def variant_sequences_dataframe_from_args(args):
    variant_sequences_generator = variant_sequences_generator_from_args(args)
    return variant_sequences_generator_to_dataframe(variant_sequences_generator)
