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
Sequence alignment helpers
"""

from __future__ import print_function, division, absolute_import


def alignment_score(a, b, min_subsequence_length=1):
    """
    Number of mismatches between all two input sequences, allows
    for trimming of ends of sequences but not insertions or deletions
    within the sequences. Number of trimmed amino acids from each sequence
    count toward the mismatch total.

    Parameters
    ----------
    a : str

    b : str

    min_subsequence_length : int
        Only consider subsequences which are at least this long.
    Returns int
    """

    n_a = len(a)
    n_b = len(b)

    # swap a and b if a is longer since the loops below expect the first
    # sequence to be shorter. This happens because any subsequence of `a`
    # is expected to be of a length that can be found within `b`
    if n_a > n_b:
        a, b = b, a
        n_a, n_b = n_b, n_a

    # compare all subsequences of a and b and count the number of mismatches
    # between
    best_score = n_a + n_b

    # if we need to make sequence of at least length `min_subsequence_length`
    # then we can't start the sequence more than that number of characters
    # from the end of the string.
    # For example, if n_a = 7 and min_subsequence_length = 6 then
    # this loop should only consider start_a = {0, 1} but not any value
    # higher than that.
    for start_a in range(n_a - min_subsequence_length + 1):
        # Similarly, only consider end indices for the subsequence of `a`
        # which span the minimum number of characters required.
        # For example, if n_a = 7, min_subsequence_length = 6, start_a = 1
        # then this loop should only consider end_a = 7.
        for end_a in range(start_a + min_subsequence_length, n_a + 1):
            subseq_a = a[start_a:end_a]
            n_subseq_a = len(subseq_a)
            n_trimmed_a = n_a - n_subseq_a
            # consider all subsequences of the second string of the same length
            # as the subsequence extracted from the first string
            for start_b in range(n_b - n_subseq_a + 1):
                subseq_b = b[start_b:start_b + n_subseq_a]
                n_subseq_b = len(subseq_b)
                assert n_subseq_a == n_subseq_b
                n_trimmed_b = n_b - n_subseq_b
                # now that we have two subsequences of the same length,
                # count up the number of mismatching characters between them
                n_mismatches = sum([
                    a_char != b_char
                    for (a_char, b_char)
                    in zip(subseq_a, subseq_b)])
                score = n_trimmed_a + n_trimmed_b + n_mismatches
                if score < best_score:
                    best_score = score
    return best_score
