from __future__ import print_function, division, absolute_import


from isovar.translation import translate_variants
from isovar.allele_read_helpers import reads_supporting_variants

from nose.tools import eq_

from testing_helpers import load_bam, load_vcf


def test_translate_variant_collection():
    variants = load_vcf("data/b16.f10/b16.vcf")
    samfile = load_bam("data/b16.f10/b16.combined.sorted.bam")

    result = list(
        translate_variants(
            reads_supporting_variants(variants, samfile)))
    eq_(
        len(result),
        4,
        "Expected %d translated variants but got %d: %s" % (
            len(variants),
            len(result),
            result))
