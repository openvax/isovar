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
Nucleotide helpers on Python strings, to avoid depending on a bigger library
such as BioPython.
"""

# IUPAC codes, including N, in both cases. Other characters are unchanged.
_COMPLEMENT = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


def complement_dna(seq):
    """
    Complement every base of a DNA sequence, preserving case.

    Parameters
    ----------
    seq : str

    Returns str
    """
    return seq.translate(_COMPLEMENT)


def reverse_complement_dna(seq):
    """
    Reverse complement of a DNA sequence

    Parameters
    ----------
    seq : str

    Returns str
    """
    return complement_dna(seq)[::-1]
