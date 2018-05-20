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

from .logging import get_logger

logger = get_logger(__name__)


def trim_N_nucleotides(prefix, suffix):
    """
    Drop all occurrences of 'N' from prefix and suffix nucleotide strings
    by trimming.
    """
    if 'N' in prefix:
        # trim prefix to exclude all occurrences of N
        rightmost_index = prefix.rfind('N')
        logger.debug(
            "Trimming %d nucleotides from read prefix '%s'",
            rightmost_index + 1, prefix)
        prefix = prefix[rightmost_index + 1:]

    if 'N' in suffix:
        leftmost_index = suffix.find('N')
        logger.debug(
            "Trimming %d nucleotides from read suffix '%s'",
            len(suffix) - leftmost_index,
            suffix)
        suffix = suffix[:leftmost_index]

    return prefix, suffix

def convert_from_bytes_if_necessary(prefix, suffix):
    """
    Depending on how we extract data from pysam we may end up with either
    a string or a byte array of nucleotides. For consistency and simplicity,
    we want to only use strings in the rest of our code.
    """
    if isinstance(prefix, bytes):
        prefix = prefix.decode('ascii')

    if isinstance(suffix, bytes):
        suffix = suffix.decode('ascii')

    return prefix, suffix
