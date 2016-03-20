from varcode import Variant
from pyensembl import ensembl_grch38
from isovar.reference_context import (
    interbase_range_affected_by_variant_on_transcript
)
from nose.tools import eq_


def test_interbase_of_brca2_utr_substitution():
    # rs769125639 is a simple T>A substitution in the 6th nucleotide of
    # BRCA2-001's 5' UTR
    brca2_variant_rs769125639 = Variant("13", 32315479, "T", "A")
    brca2_001 = ensembl_grch38.transcripts_by_name("BRCA2-001")[0]
    (start, end) = interbase_range_affected_by_variant_on_transcript(
        variant=brca2_variant_rs769125639,
        transcript=brca2_001)
    eq_((start, end), (5, 6))
