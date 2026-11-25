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
Shared output, logging, and overwrite-control helpers for Isovar CLI commands.

Conventions (#655, #658, #663):
- ``-o`` / ``--output FILE`` always refers to a single output file.
- ``--output-dir DIR`` always refers to an output directory.
- ``-f`` / ``--force`` opts into overwriting existing files or directories on
  commands that refuse to overwrite by default.
- ``-q`` / ``--quiet`` suppresses INFO progress messages on stderr alongside
  ``--log-level``.
- Legacy flag spellings are registered as hidden aliases (``argparse.SUPPRESS``)
  so existing scripts continue to work without cluttering ``--help``.
"""

import argparse
from pathlib import Path
import sys

from .validation import CommandInputError

LOG_LEVELS = ("DEBUG", "INFO", "WARNING", "ERROR")
COMPRESSED_SUFFIXES = frozenset({".gz", ".bz2", ".xz", ".zst"})
TSV_SUFFIXES = frozenset({".tsv", ".tab"})


class _QuietAction(argparse.Action):
    """Set ``quiet=True`` and raise ``log_level`` to at least WARNING."""

    def __init__(self, option_strings, dest, default=False, help=None):
        super().__init__(
            option_strings=option_strings,
            dest=dest,
            nargs=0,
            default=default,
            help=help,
        )

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, True)
        current = getattr(namespace, "log_level", "INFO")
        if current in ("DEBUG", "INFO"):
            setattr(namespace, "log_level", "WARNING")


def add_hidden_alias(parser, *flags, dest, type=str, action="store", default=argparse.SUPPRESS):
    """Register deprecated/legacy flag aliases hidden from ``--help``."""
    kwargs = {"dest": dest, "action": action, "default": default, "help": argparse.SUPPRESS}
    if action in ("store", "append"):
        kwargs["type"] = type
    return parser.add_argument(*flags, **kwargs)


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
    """Add shared ``--log-level`` and ``-q``/``--quiet`` options to ``parser``."""
    parser.add_argument(
        "--log-level",
        type=str.upper,
        choices=LOG_LEVELS,
        default="INFO",
        help="Verbosity of progress messages written to stderr (default: %(default)s)",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action=_QuietAction,
        default=False,
        help="Suppress INFO progress messages on stderr (equivalent to --log-level WARNING)",
    )
    return parser


def add_output_file_arg(
    parser,
    *,
    default=None,
    required=False,
    type=str,
    metavar="FILE",
    help=None,
    description=None,
    dest="output",
    short="-o",
    aliases=(),
    group=None,
):
    """Add canonical ``-o``/``--output FILE`` (always a single output file).

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Parser (or argument group via ``group``) to extend.
    default : Any, optional
        Default output filename when ``required`` is False.
    required : bool, optional
        Whether ``--output`` must be supplied explicitly.
    type : callable, optional
        Value conversion type (e.g., ``str`` or ``Path``).
    metavar : str, optional
        Placeholder shown in ``--help`` (default ``"FILE"``).
    help : str, optional
        Explicit help string. If omitted, derived from ``description``.
    description : str, optional
        Short description used when ``help`` is omitted.
    dest : str, optional
        Parsed attribute name (default ``"output"``).
    short : str or None, optional
        Short flag (default ``"-o"``). Pass ``None`` only when ``-o`` is
        already taken on the parser.
    aliases : sequence of str, optional
        Legacy long flag spellings registered as hidden aliases.
    group : argparse._ArgumentGroup, optional
        Optional argument group to add the option to instead of ``parser``.
    """
    target = group if group is not None else parser
    if help is None:
        base = description or "Output file"
        help = (
            base + " (default: %(default)s)"
            if default is not None and not required
            else base
        )
    flags = [short, "--output"] if short else ["--output"]
    action = target.add_argument(
        *flags,
        dest=dest,
        default=default,
        required=required,
        type=type,
        metavar=metavar,
        help=help,
    )
    for alias in aliases:
        add_hidden_alias(target, alias, dest=dest, type=type)
    return action


def add_output_dir_arg(
    parser,
    *,
    default=None,
    required=False,
    type=str,
    metavar="DIR",
    help="Output directory",
    dest="output_dir",
    aliases=(),
    group=None,
):
    """Add canonical ``--output-dir DIR`` (always an output directory)."""
    target = group if group is not None else parser
    action = target.add_argument(
        "--output-dir",
        dest=dest,
        default=default,
        required=required,
        type=type,
        metavar=metavar,
        help=help,
    )
    for alias in aliases:
        add_hidden_alias(target, alias, dest=dest, type=type)
    return action


def add_force_arg(
    parser,
    *,
    help="Overwrite existing output files or directories if present",
    group=None,
):
    """Add canonical ``-f``/``--force`` flag for overwrite-guarded commands."""
    target = group if group is not None else parser
    return target.add_argument(
        "-f",
        "--force",
        action="store_true",
        default=False,
        help=help,
    )


def add_output_args(parser, filename, description="Output table file"):
    """Add shared table output options (``-o``/``--output``, ``--output-columns``, logging)."""
    add_log_level_arg(parser)
    output_group = parser.add_argument_group("Output")
    add_output_file_arg(
        parser,
        default=filename,
        required=False,
        help=description + " (.tsv/.tab writes TSV, '-' writes stdout; default: %(default)s)",
        group=output_group,
    )
    output_group.add_argument(
        "--output-columns",
        default=None,
        nargs="+",
        help="Subset of columns to write",
    )
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
