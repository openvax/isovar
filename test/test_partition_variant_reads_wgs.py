
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
        assert read_alt == alt

if __name__ == "__main__":
    test_partition_variant_reads_snv()
