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

from varcode import load_vcf_fast

def add_somatic_vcf_args(parser):
    """
    Extends an ArgumentParser instance with the following commandline arguments:
        --vcf
        --genome
    """
    variant_group = parser.add_argument_group("Variants")
    variant_group.add_argument(
        "--vcf",
        required=True,
        help="Path to VCF file containing somatic variants")

    variant_group.add_argument(
        "--genome",
        default=None,
        required=False,
        help="Name of reference genome for VCF of somatic variants")
    return parser

def variants_from_args(args):
    return load_vcf_fast(args.vcf, genome=args.genome)
