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

from isovar import (
    unique_counts,
    partition_variant_reads,
    assemble_transcript_fragments
)
from pysam import AlignmentFile

def test_unique_counts():
    samfile = AlignmentFile("data/cancer-wgs-primary.chr12.bam")
    chromosome = "chr12"
    base1_location = 65857041
    ref = "G"
    alt = "C"

    seq_parts = partition_variant_reads(
        samfile=samfile,
        chromosome=chromosome,
        base1_location=base1_location,
        ref=ref,
        alt=alt)

    c = unique_counts(seq_parts)
    # there are some redundant reads, so we expect that the number of
    # unique entries should be less than the total read partitions
    assert len(seq_parts) > len(c)

def test_assemble_transcript_fragments_snv():
    samfile = AlignmentFile("data/cancer-wgs-primary.chr12.bam")
    chromosome = "chr12"
    base1_location = 65857041
    ref = "G"
    alt = "C"

    seq_parts = partition_variant_reads(
        samfile=samfile,
        chromosome=chromosome,
        base1_location=base1_location,
        ref=ref,
        alt=alt)

    for (p, v, s), c in assemble_transcript_fragments(seq_parts):
        print("%s%s%s weight=%d length=%d" % (
            p,
            v,
            s,
            c,
            len(p + v + s)))

if __name__ == "__main__":
    test_unique_counts()
    test_assemble_transcript_fragments_snv()
