# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from isovar.alignment_score import alignment_score

def test_alignment_score_same_strings():
    assert alignment_score("", "") == 0
    assert alignment_score("a", "a") == 0
    assert alignment_score("WUZZLE", "WUZZLE") == 0

def test_alignment_score_strings_of_different_length():
    assert alignment_score("ab", "a") == 1
    assert alignment_score("a", "ab") == 1
    assert alignment_score("WUZZLE", "WUZZLE?") == 1
    assert alignment_score("SNUZZLE", "UZZLE") == 2


def test_alignment_score_totally_different_strings():
    assert alignment_score("", "a") == 1
    assert alignment_score("", "ab") == 2
    assert alignment_score("WUZZLE", "HEAVY") == 6
    assert alignment_score("DOG", "CATCATCAT") == 9
    assert alignment_score("a", "") == 1
    assert alignment_score("ab", "") == 2
    assert alignment_score("WUZZLE", "HEAVY") == 6
    assert alignment_score("CATCATCAT", "DOG") == 9


def test_alignment_min_subsequence_length():
    # if matching subsequence isn't long enough then misalignment score
    # will be sum of two sequence lengths
    assert alignment_score("aaa", "aaa", min_subsequence_length=10) == 6
