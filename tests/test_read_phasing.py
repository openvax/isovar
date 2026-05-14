# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from varcode import Variant

from isovar import IsovarReadPhasing
from isovar.allele_read import AlleleRead
from isovar.isovar_result import IsovarResult
from isovar.read_evidence import ReadEvidence

from .common import eq_


def _make_result(variant, alt_read_names, phased_partners=()):
    read_evidence = ReadEvidence(
        trimmed_base1_start=variant.start,
        trimmed_ref=variant.ref,
        trimmed_alt=variant.alt,
        ref_reads=[],
        alt_reads=[
            AlleleRead(prefix="A", allele=variant.alt, suffix="T", name=name)
            for name in alt_read_names
        ],
        other_reads=[],
    )
    return IsovarResult(
        variant=variant,
        read_evidence=read_evidence,
        predicted_effect=None,
        phased_variants_in_supporting_reads=set(phased_partners),
    )


def test_has_evidence_true_for_variant_with_alt_reads():
    v = Variant("1", 10, "A", "C", normalize_contig_names=False)
    phasing = IsovarReadPhasing([_make_result(v, alt_read_names={"r1", "r2"})])
    assert phasing.has_evidence(v)


def test_has_evidence_false_when_no_alt_reads():
    v = Variant("1", 10, "A", "C", normalize_contig_names=False)
    phasing = IsovarReadPhasing([_make_result(v, alt_read_names=set())])
    assert not phasing.has_evidence(v)


def test_has_evidence_false_for_unknown_variant():
    v = Variant("1", 10, "A", "C", normalize_contig_names=False)
    unknown = Variant("1", 99, "G", "T", normalize_contig_names=False)
    phasing = IsovarReadPhasing([_make_result(v, alt_read_names={"r1"})])
    assert not phasing.has_evidence(unknown)


def test_partners_in_cis_returns_phased_variants():
    v1 = Variant("1", 10, "A", "C", normalize_contig_names=False)
    v2 = Variant("1", 11, "G", "T", normalize_contig_names=False)
    v3 = Variant("1", 12, "T", "G", normalize_contig_names=False)
    phasing = IsovarReadPhasing([
        _make_result(v1, alt_read_names={"r1"}, phased_partners={v2, v3}),
        _make_result(v2, alt_read_names={"r1"}, phased_partners={v1}),
        _make_result(v3, alt_read_names={"r1"}, phased_partners={v1}),
    ])
    eq_(phasing.partners_in_cis(v1), (v2, v3))
    eq_(phasing.partners_in_cis(v2), (v1,))


def test_partners_in_cis_empty_when_no_partners():
    v = Variant("1", 10, "A", "C", normalize_contig_names=False)
    phasing = IsovarReadPhasing([_make_result(v, alt_read_names={"r1"})])
    eq_(phasing.partners_in_cis(v), ())


def test_partners_in_cis_empty_for_unknown_variant():
    v = Variant("1", 10, "A", "C", normalize_contig_names=False)
    unknown = Variant("2", 50, "A", "G", normalize_contig_names=False)
    phasing = IsovarReadPhasing([_make_result(v, alt_read_names={"r1"})])
    eq_(phasing.partners_in_cis(unknown), ())


def test_partners_in_cis_is_deterministically_ordered():
    anchor = Variant("1", 10, "A", "C", normalize_contig_names=False)
    later = Variant("1", 20, "A", "C", normalize_contig_names=False)
    earlier = Variant("1", 5, "A", "C", normalize_contig_names=False)
    other_chrom = Variant("2", 1, "A", "C", normalize_contig_names=False)
    phasing = IsovarReadPhasing([
        _make_result(
            anchor,
            alt_read_names={"r1"},
            phased_partners={later, earlier, other_chrom}),
    ])
    eq_(phasing.partners_in_cis(anchor), (earlier, later, other_chrom))
