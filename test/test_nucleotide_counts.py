from isovar import most_common_nucleotides, partition_variant_reads
from pysam import AlignmentFile


def test_most_common_nucleotides_for_chr12_deletion():
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
    consensus_sequence, chosen_counts, other_counts = most_common_nucleotides(
        seq_parts)
    print(consensus_sequence)
    print(chosen_counts)
    print(other_counts)
    assert len(chosen_counts) == len(consensus_sequence)
    assert len(other_counts) == len(consensus_sequence)
    assert other_counts.sum() == 0, "Didn't expect disagreement among reads"
    for (p, _, s) in seq_parts:
        assert p in consensus_sequence
        assert s in consensus_sequence

if __name__ == "__main__":
    test_most_common_nucleotides_for_chr12_deletion()
