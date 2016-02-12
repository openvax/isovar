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

from isovar import partition_variant_reads

from pysam import AlignmentFile

def test_partition_variant_reads_snv():
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
    assert len(seq_parts) > 1
    for (prefix, read_alt, suffix) in sorted(seq_parts, key=lambda x: len(x[0])):
        print("%100s%s%-100s" % (prefix, read_alt, suffix))
        assert read_alt == alt

def test_partition_variant_reads_deletion():
    samfile = AlignmentFile("data/cancer-wgs-primary.chr12.bam")
    chromosome = "chr12"
    base1_location = 70091490
    ref = "TTGTAGATGCTGCCTCTCC"
    alt = ""
    seq_parts = partition_variant_reads(
        samfile=samfile,
        chromosome=chromosome,
        base1_location=base1_location,
        ref=ref,
        alt=alt)
    assert len(seq_parts) > 1
    for (prefix, read_alt, suffix) in sorted(seq_parts, key=lambda x: len(x[0])):
        print("%100s%s%-100s" % (prefix, read_alt, suffix))
        assert read_alt == alt


if __name__ == "__main__":
    test_partition_variant_reads_snv()
    test_partition_variant_reads_deletion()
