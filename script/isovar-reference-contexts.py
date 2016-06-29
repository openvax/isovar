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

from __future__ import print_function, division, absolute_import

import argparse

from isovar.default_parameters import CDNA_CONTEXT_SIZE
from isovar.args import (
    add_somatic_vcf_args,
    reference_contexts_dataframe_from_args
)

parser = argparse.ArgumentParser()
parser = add_somatic_vcf_args(parser)

parser.add_argument(
    "--context-size",
    default=CDNA_CONTEXT_SIZE,
    type=int)

parser.add_argument(
    "--output",
    default="isovar-reference-contexts-result.csv",
    help="Name of output CSV")

if __name__ == "__main__":
    args = parser.parse_args()
    reference_contexts_df = reference_contexts_dataframe_from_args(args)
    print(reference_contexts_df)
    reference_contexts_df.to_csv(args.output)
