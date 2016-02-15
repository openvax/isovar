from pysam import AlignmentFile
from isovar import sequence_counts, gather_variant_reads

def test_sequence_counts_snv():
    samfile = AlignmentFile("data/cancer-wgs-primary.chr12.bam")
    chromosome = "chr12"
    base1_location = 65857041
    ref = "G"
    alt = "C"

    variant_reads = gather_variant_reads(
        samfile=samfile,
        chromosome=chromosome,
        base1_location=base1_location,
        ref=ref,
        alt=alt)
    print(sequence_counts(variant_reads))

if __name__ == "__main__":
    test_sequence_counts_snv()
