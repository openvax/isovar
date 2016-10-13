#!/usr/bin/env python

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

"""
Prints number of reads supporting ref, alt, and other alleles at variant loci.
"""

from __future__ import division, absolute_import
import logging
import logging.config
import pkg_resources

from isovar.cli.rna_reads import (
    make_rna_reads_arg_parser,
    allele_reads_generator_from_args
)
from isovar.allele_counts import allele_counts_dataframe


logging.config.fileConfig(pkg_resources.resource_filename('isovar.cli', 'logging.conf'))
logger = logging.getLogger(__name__)

parser = make_rna_reads_arg_parser()
parser.add_argument(
    "--output",
    default="isovar-allele-counts-result.csv",
    help="Name of CSV file which contains read sequences")

if __name__ == "__main__":
    args = parser.parse_args()
    logger.info(args)
    variants_and_allele_reads_generator = allele_reads_generator_from_args(args)
    allele_counts_df = allele_counts_dataframe(variants_and_allele_reads_generator)
    logger.info(allele_counts_df)
    allele_counts_df.to_csv(args.output)
