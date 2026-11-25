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
Common helper functions for writing CSV/TSV output tables, shared by all
table-export CLI commands.
"""

from pathlib import Path
import sys

from .validation import CommandInputError

LOG_LEVELS = ("DEBUG", "INFO", "WARNING", "ERROR")
COMPRESSED_SUFFIXES = frozenset({".gz", ".bz2", ".xz", ".zst"})
TSV_SUFFIXES = frozenset({".tsv", ".tab"})


def table_separator_for_path(path):
    """Choose delimiter from an output table path.

    Paths ending in ``.tsv`` or ``.tab`` (optionally followed by a compression
    suffix such as ``.gz``) use tab delimiters; all other paths (including
    ``"-"`` for stdout) use commas.

    Parameters
    ----------
    path : str or os.PathLike
        Destination file path or ``"-"`` for stdout.

    Returns
    -------
    str
        ``"\\t"`` for TSV/TAB paths, ``","`` otherwise.
    """
    text = str(path)
    if text == "-":
        return ","
    suffixes = [suffix.lower() for suffix in Path(text).suffixes]
    while suffixes and suffixes[-1] in COMPRESSED_SUFFIXES:
        suffixes.pop()
    if suffixes and suffixes[-1] in TSV_SUFFIXES:
        return "\t"
    return ","


def add_log_level_arg(parser):
    parser.add_argument(
        "--log-level",
        type=str.upper,
        choices=LOG_LEVELS,
        default="INFO",
        help="Verbosity of progress messages written to stderr (default: %(default)s)")
    return parser


def add_output_args(parser, filename, description="Output table file"):
    add_log_level_arg(parser)
    output_group = parser.add_argument_group("Output")
    output_group.add_argument(
        "--output",
        default=filename,
        help=description + " (.tsv/.tab writes TSV, '-' writes stdout; default: %(default)s)")
    output_group.add_argument(
        "--output-columns",
        default=None,
        nargs="+",
        help="Subset of columns to write")
    return parser


def write_dataframe(df, args):
    """
    Write a DataFrame to location specified in commandline arguments,
    optionally filtered by specific columns.

    Infers tab vs. comma delimiter from ``args.output`` (``.tsv`` / ``.tab``
    write TSV) and writes to ``sys.stdout`` when ``args.output == "-"``.
    """
    if args.output_columns is not None and len(args.output_columns) > 0:
        missing = [col for col in args.output_columns if col not in set(df.columns)]
        if missing:
            raise CommandInputError("Unknown --output-columns %s; valid columns: %s" % (
                ", ".join(missing), ", ".join(df.columns)))
        df = df[args.output_columns]
    sep = table_separator_for_path(args.output)
    target = sys.stdout if str(args.output) == "-" else args.output
    df.to_csv(target, sep=sep, index=False, float_format="%.6g", lineterminator="\n")
