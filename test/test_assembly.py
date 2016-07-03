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

from isovar.variant_reads import reads_supporting_variant
from isovar.assembly import assemble_reads_into_variant_sequences
from pyensembl import ensembl_grch38
from varcode import Variant
from nose.tools import eq_

from testing_helpers import load_bam

def test_assemble_transcript_fragments_snv():
    samfile = load_bam("data/cancer-wgs-primary.chr12.bam")
    chromosome = "chr12"
    base1_location = 65857041
    ref = "G"
    alt = "C"
    variant = Variant(
        contig=chromosome,
        start=base1_location,
        ref=ref,
        alt=alt,
        ensembl=ensembl_grch38)
    variant_reads = reads_supporting_variant(
        variant=variant,
        samfile=samfile,
        chromosome=chromosome,)

    sequences = assemble_reads_into_variant_sequences(
        variant_reads,
        min_overlap_size=30)

    assert len(sequences) > 0
    max_read_length = max(len(r) for r in variant_reads)
    for s in sequences:
        print("%s%s%s weight=%d length=%d" % (
            s.prefix,
            s.alt,
            s.suffix,
            len(s.reads),
            len(s.sequence)))
        eq_(s.alt, alt)
        assert len(s) > max_read_length, \
            "Expected assembled sequences to be longer than read length (%d)" % (
                max_read_length,)

if __name__ == "__main__":
    test_assemble_transcript_fragments_snv()
