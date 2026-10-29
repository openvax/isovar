"""Background scope is explicit and source header resolution preserves exclusions."""

import pysam
from varcode import Variant

from isovar.normal_regions import normal_background_regions, resolve_normal_regions


def test_regions_merge_flanks_without_replacing_original_contig_names():
    variants = [Variant("chr1", 100, "A", "G"), Variant("chrM", 10, "A", "AT")]
    assert normal_background_regions(variants, background_flank=2,
        additional_regions=[("chr1", 101, 110), ("chr1", 800, 900)]) == (
            ("chr1", 97, 110), ("chr1", 800, 900), ("chrM", 7, 13))


def test_region_resolution_handles_aliases_ambiguity_and_declared_bounds():
    header = pysam.VariantHeader()
    for contig in ("1", "CHR1", "chrM", "HPV16", "chrUn_KI270442v1"):
        header.contigs.add(contig, length=100)
    result = resolve_normal_regions([("chr1", 0, 10), ("MT", 90, 120), ("HPV16", 100, 110),
        ("HPV16", 10, 20), ("KI270442.1", 20, 30)], header, reference_name="GRCh38")
    assert result.regions == (("HPV16", 10, 20), ("chrM", 90, 100), ("chrUn_KI270442v1", 20, 30))
    assert len(result.unresolved_contigs) == 1
    assert result.clipped_regions == (("MT", 90, 120),)
    assert result.outside_contig_regions == (("HPV16", 100, 110),)
