from nose.tools import eq_
from pysam import AlignmentFile
import skbio

from isovar.variant_sequences import variant_reads_to_sequences
from isovar.variant_reads import gather_reads_for_single_variant


def test_sequence_counts_snv():
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
    result = variant_reads_to_sequences(
        variant_reads,
        context_size=45)
    print(result)
    eq_(result.variant_nucleotides, alt)
    eq_(len(result.combined_sequence_weights), 1)
    eq_(len(result.full_read_counts), 1)

    for ((prefix, suffix), weight) in sorted(
            result.combined_sequence_weights.items(),
            key=lambda x: x[1]):
        eq_(len(prefix), 45)
        eq_(len(suffix), 45)
        variant = result.variant_nucleotides
        print("%s|%s|%s: %f" % (
            prefix,
            variant,
            suffix,
            weight))

        # translate in three reading frames:
        seq = "%s%s%s" % (prefix, variant, suffix)
        for offset in range(3):
            dna = skbio.DNA(seq[offset:])
            print("frame=%d: %s" % (offset, dna.translate()))


if __name__ == "__main__":
    test_sequence_counts_snv()
