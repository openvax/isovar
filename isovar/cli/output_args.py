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

from .validation import CommandInputError

LOG_LEVELS = ("DEBUG", "INFO", "WARNING", "ERROR")


def add_log_level_arg(parser):
    parser.add_argument(
        "--log-level",
        type=str.upper,
        choices=LOG_LEVELS,
        default="INFO",
        help="Verbosity of progress messages written to stderr (default: %(default)s)")
    return parser


def add_output_args(parser, filename, description="Output CSV file"):
    add_log_level_arg(parser)
    output_group = parser.add_argument_group("Output")
    output_group.add_argument(
        "--output",
        default=filename,
        help=description + " (default: %(default)s)")
    output_group.add_argument(
        "--output-columns",
        default=None,
        nargs="+",
        help="Subset of columns to write")
    return parser


def write_dataframe(df, args):
    """
    Write a DataFrame to location specified in commandline arguments,
    optionally filtered by specific columns
    """
    if args.output_columns is not None and len(args.output_columns) > 0:
        missing = [col for col in args.output_columns if col not in set(df.columns)]
        if missing:
            raise CommandInputError("Unknown --output-columns %s; valid columns: %s" % (
                ", ".join(missing), ", ".join(df.columns)))
        df = df[args.output_columns]
    df.to_csv(args.output, index=False, float_format="%.6g")
