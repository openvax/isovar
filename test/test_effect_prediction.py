from isovar.effect_prediction import (
    top_varcode_effect,
    reference_coding_transcripts_for_variant
)
from varcode import Variant
from varcode.effects import Intergenic
from nose.tools import eq_

# intergenic variant from error log of using Isovar 1.0.0
intergenic_variant = Variant('1', 30256419, 'G', 'T', 'GRCh37')


def test_top_effect_outside_of_gene():
    top_effect = top_varcode_effect(intergenic_variant)
    assert top_effect is not None
    assert isinstance(top_effect, Intergenic)


def test_reference_coding_transcripts_outside_of_gene():
    transcripts = \
        reference_coding_transcripts_for_variant(
            intergenic_variant)
    eq_(len(transcripts), 0)

