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
Prints all supporting read sequences containing a variant, along with their name
"""

from __future__ import print_function, division, absolute_import
import argparse

import varcode
from pysam import AlignmentFile

from isovar.variant_read import variant_reads_dataframe
from isovar.args import (
    extend_parser_with_somatic_vcf_args,
)

parser = argparse.ArgumentParser()

extend_parser_with_somatic_vcf_args(parser)

parser.add_argument(
    "--output",
    default="isovar-variant-reads-result.csv",
    help="Name of CSV file which contains read sequences")

if __name__ == "__main__":
    args = parser.parse_args()

    print(args)

    variants = varcode.load_vcf(
        args.vcf,
        genome=args.genome)

    samfile = AlignmentFile(args.bam)

    df = variant_reads_dataframe(variants, samfile)

    print(df)

    df.to_csv(args.output)
