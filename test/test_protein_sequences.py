
from pysam import AlignmentFile
from varcode import load_vcf
from isovar import variants_to_protein_sequences

VCF = "/Users/iskander/code/varlens/test/data/CELSR1/vcfs/vcf_1.vcf"

BAM = "/Users/iskander/code/varlens/test/data/CELSR1/bams/bam_1.bam"

GENOME = "hg19"

def test_variants_to_protein_sequences():
    variants = load_vcf(VCF, genome=GENOME)
    samfile = AlignmentFile(BAM)
    result = variants_to_protein_sequences(variants, samfile)
    assert len(result) > 0


if __name__ == "__main__":
    test_variants_to_protein_sequences()
