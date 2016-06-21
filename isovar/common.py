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
import logging

from collections import defaultdict

nucleotides = ["A", "C", "T", "G"]
nucleotide_to_index = {c: i for (i, c) in enumerate(nucleotides)}
index_to_nucleotide = {i: c for (i, c) in enumerate(nucleotides)}

def create_logger(
        name="root",
        level=logging.DEBUG,
        format_string="%(levelname)s: %(message)s (%(filename)s:%(lineno)s - %(funcName)s)"):
    logger = logging.getLogger(name)
    handler = logging.StreamHandler()
    formatter = logging.Formatter(format_string)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(level)
    return logger

logger = create_logger()

def group_unique_sequences(
        allele_reads,
        max_prefix_size=None,
        max_suffix_size=None):
    """
    Given a list of VariantRead objects, extracts all unique
    (prefix, allele, suffix) sequences and associate each with a list
    of reads that contained that sequence.
    """
    groups = defaultdict(set)
    for r in allele_reads:
        prefix = r.prefix
        allele = r.allele
        suffix = r.suffix
        if max_prefix_size and len(prefix) > max_prefix_size:
            prefix = prefix[-max_prefix_size:]
        if max_suffix_size and len(suffix) > max_suffix_size:
            suffix = suffix[:max_suffix_size]
        key = (prefix, allele, suffix)
        groups[key].add(r)
    return groups

def count_unique_sequences(
        allele_reads,
        max_prefix_size=None,
        max_suffix_size=None):
    """
    Given a list of AlleleRead objects, extracts all unique
    (prefix, allele, suffix) sequences and associate each with the number
    of reads that contain that sequence.
    """
    groups = group_unique_sequences(
        allele_reads,
        max_prefix_size=max_prefix_size,
        max_suffix_size=max_suffix_size)
    return {
        seq_tuple: len(read_names)
        for (seq_tuple, read_names) in groups.items()
    }

"""
def group_reads_by_allele(allele_reads, **kwargs):
    Given a sequence of AlleleRead objects, group them into a dictionary
    keyed by allele.

    result = defaultdict(list)
    for (prefix, allele, suffix), reads in group_unique_sequences(allele_reads, **kwargs):
        result[allele].append((prefix, suffix))
    pairs = [(r.prefix, r.suffix) for r in variant_reads]
    return variant_seq, pairs
"""


def list_to_string(list_of_anything, sep=";"):
    """
    Helper function used for building the fields of a printable dataframe
    """
    return sep.join(str(x) for x in list_of_anything)
