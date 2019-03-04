from __future__ import print_function, division, absolute_import

import pysam


class DummySamFile(object):
    """
    Used instead of real AlignmentFile objects for test.
    """
    def __init__(self, reads):
        self.reads = reads

    def fetch(self, *args, **kwargs):
        return self.reads


def make_read(
        seq,
        cigar,
        mdtag=None,
        name="dummy",
        mapq=10,
        baseq=30,
        reference_start=0,
        reference_id=0):
    read = pysam.AlignedSegment()
    read.seq = seq
    read.cigarstring = cigar
    if mdtag:
        read.set_tag("MD", mdtag)
    read.qname = name
    read.mapq = mapq
    read.reference_start = reference_start
    read.reference_id = reference_id
    qualities_string = pysam.qualities_to_qualitystring([baseq] * len(seq))
    read.qual = qualities_string.encode("ascii")
    return read
