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
Common command-line arguments for filtering Isovar results
"""

from __future__ import print_function, division, absolute_import

from collections import OrderedDict

from ..default_parameters import (
    MIN_NUM_RNA_ALT_FRAGMENTS,
    MIN_NUM_RNA_ALT_READS,
    MIN_FRACTION_RNA_ALT_FRAGMENTS,
    MIN_RATIO_RNA_ALT_TO_OTHER_FRAGMENTS
)

def add_filter_args(parser):
    """
    Extends an ArgumentParser instance with commandline arguments related
    to filtering variants and/or their associated protein sequences.
    """
    filter_group = parser.add_argument_group("Filtering")
    filter_group.add_argument(
        "--min-alt-rna-reads",
        type=int,
        default=MIN_NUM_RNA_ALT_READS,
        help="Minimum number of reads supporting variant allele (default %(default)s)")

    filter_group.add_argument(
        "--min-alt-rna-fragments",
        type=int,
        default=MIN_NUM_RNA_ALT_FRAGMENTS,
        help=(
            "Minimum number of fragments supporting variant allele (default %(default)s). "
            "Note that this option is the same as --min-alt-rna-reads for single-end "
            "sequencing."))

    filter_group.add_argument(
        "--min-alt-rna-fraction",
        type=float,
        default=MIN_FRACTION_RNA_ALT_FRAGMENTS,
        help=(
            "Minimum ratio of fragments supporting variant allele to total RNA fragments "
            "(default %(default)s)."))

    filter_group.add_argument(
        "--min-ratio-alt-to-other-fragments",
        type=float,
        default=MIN_RATIO_RNA_ALT_TO_OTHER_FRAGMENTS,
        help=(
            "At loci where alleles other than the ref and a single alt are supported, "
            "this parameter controls how many more times fragments supporting "
            "the variant allele are required relative to other non-reference "
            "alleles (default %(default)s)."))
    return filter_group


def filter_threshold_dict_from_args(args):
    """
    Convert names of filters from external CLI options to more
    consistent naming scheme of {min|max}_{Isovar property}.

    Returns OrderedDict
    """
    d = OrderedDict()
    d["min_ratio_alt_to_other_fragments"] = args.min_ratio_alt_to_other_fragments
    d["min_fraction_alt_fragments"] = args.min_alt_rna_fraction
    d["min_num_alt_fragments"] = args.min_alt_rna_fragments
    d["min_num_alt_reads"] = args.min_alt_rna_reads
    return d
