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
Reject bad command-line input with a one-line error before any reads are
processed: argparse types for numeric ranges, and loaders which report
missing or unusable inputs as CommandInputError.
"""

from argparse import ArgumentTypeError
import logging
import os

from varcode.cli import variant_collection_from_args as _varcode_variant_collection_from_args

logger = logging.getLogger(__name__)


class CommandInputError(ValueError):
    """A problem with a command's inputs, reported as a usage error."""


def _bounded(cast, check, description):
    def parse(text):
        try:
            value = cast(text)
        except ValueError:
            raise ArgumentTypeError("invalid %s value: %r" % (cast.__name__, text)) from None
        if not check(value):
            raise ArgumentTypeError("%r is not %s" % (text, description))
        return value
    parse.__name__ = cast.__name__
    return parse


positive_int = _bounded(int, lambda value: value > 0, "a positive integer")
non_negative_int = _bounded(int, lambda value: value >= 0, "a non-negative integer")
non_negative_float = _bounded(float, lambda value: value >= 0, "a non-negative number")
fraction = _bounded(float, lambda value: 0 <= value <= 1, "between 0 and 1")
open_fraction = _bounded(float, lambda value: 0 < value < 1, "strictly between 0 and 1")
dpi = _bounded(int, lambda value: value >= 72, "a resolution of at least 72 dpi")


def variant_collection_from_args(args):
    """Load the requested small variants, reporting input problems as CommandInputError.

    Structural-variant records (symbolic or breakend ALTs) are skipped with a
    warning: they need ``isovar sv-source`` and ``isovar sv-audit``.
    """
    try:
        variants = _varcode_variant_collection_from_args(args)
    except (OSError, ValueError) as error:
        message = str(error)
        if "genome" in message.lower() and "--genome" not in message:
            message += " (specify the reference with --genome)"
        raise CommandInputError(message) from error
    structural = [v for v in variants if getattr(v, "is_structural", False)]
    if not structural:
        return variants
    examples = ", ".join(v.short_description for v in structural[:3])
    if len(structural) == len(variants):
        raise CommandInputError(
            "All %d input variants are structural variants (%s); use isovar sv-source and "
            "isovar sv-audit for them" % (len(structural), examples))
    logger.warning("Skipping %d structural variant record(s) (%s%s); use isovar sv-source and "
                   "isovar sv-audit for them", len(structural), examples,
                   ", ..." if len(structural) > 3 else "")
    return variants.filter(lambda v: not getattr(v, "is_structural", False))


def germline_variants_from_args(args, variants):
    """Matched germline variants from --germline-vcf, with the somatic variants' genome."""
    path = getattr(args, "germline_vcf", None)
    if path is None:
        return None
    from varcode import load_vcf
    genome = next(iter(variants)).genome if len(variants) else getattr(args, "genome", None)
    try:
        return load_vcf(path, genome=genome)
    except (OSError, ValueError) as error:
        raise CommandInputError("Cannot load --germline-vcf %s: %s" % (path, error)) from error


def check_output_path(path):
    directory = os.path.dirname(os.path.abspath(path))
    if not os.path.isdir(directory):
        raise CommandInputError("Output directory does not exist: %s" % directory)
    if os.path.isdir(path):
        raise CommandInputError("Output is a directory, not a file: %s" % path)


def _named_files(value):
    values = value if isinstance(value, (list, tuple)) else [value]
    return [str(v) for v in values if isinstance(v, (str, os.PathLike)) and os.path.isfile(v)]


def check_outputs(args, *outputs, directories=()):
    """Check output paths before any work, so a run never fails or clobbers at the end.

    ``outputs`` are file paths (or ``args`` attribute names holding them) and
    ``directories`` output directories. Each file's directory must exist, no
    output may be an existing directory (or a directory an existing file), the
    outputs must be distinct, and none may be a file given in any other option.
    """
    def resolve(value):
        return getattr(args, value) if isinstance(value, str) and hasattr(args, value) else value

    names = {name for name in outputs + tuple(directories) if isinstance(name, str) and hasattr(args, name)}
    files = [str(path) for path in map(resolve, outputs) if path is not None]
    folders = [str(path) for path in map(resolve, directories) if path is not None]
    for path in files:
        check_output_path(path)
    for path in folders:
        if os.path.exists(path) and not os.path.isdir(path):
            raise CommandInputError("Output directory is an existing file: %s" % path)
    absolute = [os.path.abspath(path) for path in files + folders]
    if len(set(absolute)) != len(absolute):
        raise CommandInputError("Outputs must be distinct files: %s" % ", ".join(files + folders))
    inputs = [path for name, value in sorted(vars(args).items()) if name not in names
              for path in _named_files(value)]
    for output in files:
        if os.path.exists(output) and any(os.path.samefile(output, path) for path in inputs):
            raise CommandInputError("Output %s would overwrite an input file" % output)


def input_error(parser, error):
    """Report a problem with a command's inputs: one line, exit status 2, no usage block."""
    parser.exit(2, "%s: error: %s\n" % (parser.prog, error))
