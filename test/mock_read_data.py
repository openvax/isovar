from __future__ import print_function, division, absolute_import

import pysam


class DummyPileupElement(object):
    def __init__(self, alignment, is_refskip, is_del):
        self.alignment = alignment
        self.is_del = is_del
        self.is_refskip = is_refskip


class DummyPileupColumn(object):
    def __init__(self, pos, reads, is_del=False, is_refskip=False):
        self.pos = pos
        self.pileups = [
            DummyPileupElement(read, is_del=is_del, is_refskip=is_refskip)
            for read in reads]


class DummySamFile(object):
    """
    Used instead of real AlignmentFile objects for test.
    """
    def __init__(self, reads):
        self.reads = reads

    def fetch(self, *args, **kwargs):
        return self.reads

    def pileup(self, chromosome, start, end):
        for i in range(start, end + 1):
            yield DummyPileupColumn(pos=i, reads=self.reads)


def make_read(seq, cigar, mdtag=None, name="dummy", mapq=10, baseq=30):
    read = pysam.AlignedSegment()
    read.seq = seq
    read.cigarstring = cigar
    if mdtag:
        read.set_tag("MD", mdtag)
    read.qname = name
    read.mapq = mapq
    qualities_string = pysam.qualities_to_qualitystring([baseq] * len(seq))
    qualities_bytes = qualities_string.encode("ascii")
    read.qual = qualities_bytes
    return read
