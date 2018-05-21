# Copyright (c) 2018. Mount Sinai School of Medicine
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
Common helper functions for writing CSV output files, shared by all
the CLI commands
"""

from __future__ import print_function, division, absolute_import


def add_output_args(
        parser,
        filename="output.csv",
        description="Output CSV file"):
    output_group = parser.add_argument_group("Output")
    output_group.add_argument(
        "--output",
        default=filename,
        help=description)
    output_group.add_argument(
        "--output-columns",
        default=None,
        nargs="+",
        help="Subset of columns to write")
    return parser


def write_dataframe(df, args):
    assert len(args.output) > 0
    if args.output_columns is not None and len(args.output_columns) > 0:
        valid_columns = set(df.columns)
        for col in args.output_columns:
            if col not in valid_columns:
                raise ValueError("Column not found '%s', valid options: %s" % (
                    col, list(df.columns)))
        df = df[args.output_columns]
    df.to_csv(args.output, index=False)
