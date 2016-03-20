from pysam import AlignmentFile
from isovar import gather_reads_for_single_variant, group_unique_sequences

def test_group_unique_sequences():
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

    groups = group_unique_sequences(variant_reads)
    # there are some redundant reads, so we expect that the number of
    # unique entries should be less than the total read partitions
    assert len(variant_reads) > len(groups)

if __name__ == "__main__":
    test_group_unique_sequences()
