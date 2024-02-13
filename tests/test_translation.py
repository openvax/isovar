from __future__ import print_function, division, absolute_import


from isovar import ProteinSequenceCreator, ReadCollector

from nose.tools import eq_

from testing_helpers import load_bam, load_vcf


def test_translate_variant_collection():
    variants = load_vcf("data/b16.f10/b16.vcf")
    samfile = load_bam("data/b16.f10/b16.combined.sorted.bam")
    read_evidence_gen = ReadCollector().read_evidence_generator(
        variants,
        samfile)
    translation_gen = ProteinSequenceCreator().translate_variants(read_evidence_gen)
    translations = list(translation_gen)
    eq_(
        len(translations),
        4,
        "Expected %d translated variants but got %d: %s" % (
            len(variants),
            len(translations),
            translations))
