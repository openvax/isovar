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

from collections import namedtuple

OverlappingRead = namedtuple(
    "Read",
    "sequence reference_positions locus_offset name")

def gather_overlapping_reads(
        samfile,
        chromosome,
        base0_start,
        base0_end,
        is_del):
    """
    Generator that yields a sequence of OverlappingRead records for reads which
    overlap the given locus and have a matching alignment on the first position.

    Parameters
    ----------
    samfile : pysam.AlignmentFile

    chromosome : str

    base0_start : int
        Starting position of the locus

    base0_end : int
        End position of the locus

    is_del : bool
        Does this locus originate from a deletion variant?
        (if so, we want to include reads that have deletions at the locus)

    Yields OverlappingRead objects contain the following information:
        - nucleotide sequence of the read
        - list of reference position alignments for each nucleotide
            (where a non-aligned nucleotide is given by None)
        - offset of the locus within the read

    """
    base1_start = base0_start + 1
    # Let pysam pileup the reads covering our location of interest for us
    for column in samfile.pileup(
            chromosome,
            base0_start,
            base0_end):
        if column.pos != base1_start:
            continue
        for i, pileup_element in enumerate(column.pileups):
            if pileup_element.is_refskip:
                # if read sequence doesn't actually align here, skip it
                continue
            elif pileup_element.is_del and not is_del:
                # if read has a deletion at this location and variant isn't a
                # deletion
                continue
            read = pileup_element.alignment
            reference_positions = read.get_reference_positions(
                full_length=False)
            if base0_start in reference_positions:
                offset = reference_positions.index(base0_start)
                yield OverlappingRead(
                    sequence=read.query_sequence,
                    reference_positions=reference_positions,
                    locus_offset=offset,
                    name=read.query_name)
