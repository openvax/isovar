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
    gather_reads_for_single_variant,
    assemble_transcript_fragments
)
from pysam import AlignmentFile


def test_assemble_transcript_fragments_snv():
    samfile = AlignmentFile("data/cancer-wgs-primary.chr12.bam")
    chromosome = "chr12"
    base1_location = 65857041
    ref = "G"
    alt = "C"

    variant_reads = gather_reads_for_single_variant(
        samfile=samfile,
        chromosome=chromosome,
        base1_location=base1_location,
        ref=ref,
        alt=alt)

    longer_sequences = assemble_transcript_fragments(variant_reads, min_overlap_size=30)
    for (p, v, s), c in longer_sequences:
        print("%s%s%s weight=%d length=%d" % (
            p,
            v,
            s,
            c,
            len(p + v + s)))
    assert len(longer_sequences) == 1, "Expected unique sequence for locus"

if __name__ == "__main__":
    test_assemble_transcript_fragments_snv()
