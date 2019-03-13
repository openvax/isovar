# Copyright (c) 2016-2019. Mount Sinai School of Medicine
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

"""
Moving functions here which will be deleted in future versions of Isovar.
"""

from __future__ import print_function, division, absolute_import

from warnings import warn

from .read_creator import ReadCreator

def reads_supporting_variant(
        variant,
        samfile,
        chromosome=None):
    """
    Gather AlleleRead objects supporting variant in a given
    alignments (BAM/SAM) file. Uses default settings for
    gathering and filtering reads.

    Parameters
    ----------
    variant : varcode.Variant

    samfile : pysam.AlignmentFile

    chromosome : str or None

    Returns
    -------
    List of AlleleRead objects

    """
    warn("Use ReadCreator.allele_reads_supporting_variant instead")
    if chromosome is not None:
        raise NotImplementedError(
            "Use of manual chromosome argument no longer supported for this function")
    read_creator = ReadCreator()
    return read_creator.allele_reads_supporting_variant(
        variant=variant,
        alignments=samfile)
