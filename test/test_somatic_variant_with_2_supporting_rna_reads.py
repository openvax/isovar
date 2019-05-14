from __future__ import print_function, division, absolute_import

from varcode import Variant
from testing_helpers import load_bam
from nose.tools import eq_

from isovar.read_collector import ReadCollector

from genomes_for_testing import grch38


def test_somatic_variant_with_2_supporting_rna_reads():
    variant = Variant("14", 105849746, "G", "A", grch38)
    base_dir = "data/somatic-variant-with-2-supporting-rna-reads/"
    normal_reads = load_bam(base_dir + "normal.14.105849746.G.A.no-alt.sorted.bam")
    tumor_reads = load_bam(base_dir + "tumor.14.105849746.G.A.many-alt.sorted.bam")
    rna_reads = load_bam(base_dir + "rna.14.105849746.G.A.2-alt.sorted.bam")
    read_creator = ReadCollector()
    normal_sample_variant_reads = read_creator.allele_reads_supporting_variant(
        variant=variant,
        alignment_file=normal_reads)
    eq_(len(normal_sample_variant_reads), 0)
    print(normal_sample_variant_reads)

    tumor_sample_variant_reads = read_creator.allele_reads_supporting_variant(
        variant=variant,
        alignment_file=tumor_reads)
    print(tumor_sample_variant_reads)
    eq_(len(tumor_sample_variant_reads), 8)

    rna_sample_variant_reads = read_creator.allele_reads_supporting_variant(
        variant=variant,
        alignment_file=rna_reads)
    print(rna_sample_variant_reads)
    eq_(len(rna_sample_variant_reads), 2)
    # Arun went through the hassle of pulling out the exact read names
    # in IGV
    expected_variant_rna_read_names = {
        "K00193:50:H5NKVBBXX:5:2202:6421:24964",
        "K00193:50:H5NKVBBXX:5:2119:30908:1138",
    }
    for variant_read in rna_sample_variant_reads:
        assert variant_read.name in expected_variant_rna_read_names

if __name__ == "__main__":
    test_somatic_variant_with_2_supporting_rna_reads()
