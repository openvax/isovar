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

from .common import group_unique_sequences
from .assembly import collapse_substrings

def sequence_counts(variant_reads, context_size=45):
    """
    Returns a dictionary mapping (prefix, variant, suffix) to a pair of
    integers: the first indicating the number of reads which fully support
    the sequence and the second is the number of reads which partially
    support the sequence (weighted by the fraction of nucleotides
    overlapping).
    """
    unique_sequence_groups = group_unique_sequences(
        variant_reads,
        max_prefix_size=context_size,
        max_suffix_size=context_size)

    collapsed_groups = collapse_substrings(unique_sequence_groups)

    print(len(unique_sequence_groups), unique_sequence_groups)
    print(len(collapsed_groups), collapsed_groups)
