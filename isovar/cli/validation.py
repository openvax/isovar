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
import os

from varcode.cli import variant_collection_from_args as _varcode_variant_collection_from_args


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
dpi = _bounded(int, lambda value: value >= 72, "a resolution of at least 72 dpi")


def variant_collection_from_args(args):
    """Load the requested variants, reporting input problems as CommandInputError."""
    try:
        return _varcode_variant_collection_from_args(args)
    except (OSError, ValueError) as error:
        message = str(error)
        if "genome" in message.lower() and "--genome" not in message:
            message += " (specify the reference with --genome)"
        raise CommandInputError(message) from error


def check_output_path(path):
    directory = os.path.dirname(os.path.abspath(path))
    if not os.path.isdir(directory):
        raise CommandInputError("Output directory does not exist: %s" % directory)
