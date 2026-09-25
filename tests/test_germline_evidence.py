"""Cis and trans between somatic and matched germline variants from co-covering reads (#387)."""

import pysam
import pytest
from varcode import MolecularPhaseResolver, Variant

from isovar import IsovarReadPhasing, ReadCollector, run_isovar
from isovar.germline_evidence import _reference_blocks, germline_read_evidence
from isovar.isovar_result import IsovarResult
from tests.test_sv_rna import record, write_bam

REFERENCE = "CATCGGATCCTAGGCTTACAGTCCGATTGCAGTCAGGATCC"  # chr1:1001-1041
SOMATIC = Variant("1", 1011, REFERENCE[10], "G", normalize_contig_names=False)
GERMLINE = Variant("1", 1031, REFERENCE[30], "C", normalize_contig_names=False)
# Outside every read: never examined, so never phased.
DISTANT = Variant("1", 1500, "G", "T", normalize_contig_names=False)


def read(name, somatic, germline):
    bases = list(REFERENCE[:40])
    bases[10], bases[30] = somatic, germline
    return record(name, "1", 1000, "40M", "".join(bases))


def results_for(tmp_path, reads, germline=(GERMLINE, DISTANT)):
    path = write_bam(tmp_path / "rna.bam", reads)
    with pysam.AlignmentFile(str(path)) as bam:
        return run_isovar([SOMATIC], bam, germline_variants=list(germline))


def test_reads_with_the_germline_reference_allele_are_trans(tmp_path):
    reads = ([read("s%d" % i, "G", REFERENCE[30]) for i in range(3)]    # somatic alt, germline ref
             + [read("g%d" % i, REFERENCE[10], "C") for i in range(2)])  # somatic ref, germline alt
    result, = results_for(tmp_path, reads)
    evidence = result.germline_read_evidence
    assert list(evidence) == [GERMLINE]
    assert evidence[GERMLINE].ref_read_names == {"s0", "s1", "s2"}
    phasing = IsovarReadPhasing([result])
    assert phasing.in_cis(SOMATIC, GERMLINE) is False
    assert phasing.in_cis(GERMLINE, SOMATIC) is False
    assert phasing.in_cis(SOMATIC, DISTANT) is None
    # Varcode's resolver now gets trans for the germline variant.
    assert MolecularPhaseResolver(phasing).in_cis(SOMATIC, GERMLINE) is False


def test_reads_with_both_alt_alleles_are_cis(tmp_path):
    reads = [read("c%d" % i, "G", "C") for i in range(3)] + [read("r0", REFERENCE[10], REFERENCE[30])]
    result, = results_for(tmp_path, reads)
    assert IsovarReadPhasing([result]).in_cis(SOMATIC, GERMLINE) is True


def test_too_few_shared_fragments_stay_unknown(tmp_path):
    result, = results_for(tmp_path, [read("s0", "G", REFERENCE[30])])
    phasing = IsovarReadPhasing([result])
    assert phasing.in_cis(SOMATIC, GERMLINE) is None
    assert IsovarReadPhasing([result], min_shared_fragments_for_phasing=1).in_cis(SOMATIC, GERMLINE) is False


def test_a_germline_variant_in_a_skipped_intron_is_not_examined(tmp_path):
    # Each read aligns 15 bases, skips 20 (over the germline locus), then aligns 5.
    sequence = REFERENCE[:10] + "G" + REFERENCE[11:15] + REFERENCE[35:40]
    reads = [record("n%d" % i, "1", 1000, "15M20N5M", sequence) for i in range(3)]
    result, = results_for(tmp_path, reads)
    assert result.germline_read_evidence == {}
    assert IsovarReadPhasing([result]).in_cis(SOMATIC, GERMLINE) is None


def test_without_germline_variants_nothing_extra_is_collected(tmp_path):
    result, = results_for(tmp_path, [read("s0", "G", REFERENCE[30])], germline=())
    assert result.germline_read_evidence == {}


def test_a_germline_variant_that_is_also_an_input_is_not_collected_twice(tmp_path):
    path = write_bam(tmp_path / "rna.bam", [read("s%d" % i, "G", "C") for i in range(2)])
    with pysam.AlignmentFile(str(path)) as bam:
        collector = ReadCollector()
        results = [IsovarResult(variant=v, read_evidence=collector.read_evidence_for_variant(v, bam),
                                predicted_effect=None) for v in (SOMATIC, GERMLINE)]
        assert germline_read_evidence(results, [GERMLINE], bam, collector) == [{}, {}]


@pytest.mark.parametrize("cigar,blocks", [
    ("40M", [(1000, 1040)]),
    ("10M5I10M", [(1000, 1020)]),
    ("5S10M2D10M", [(1000, 1022)]),
    ("10M100N10M", [(1000, 1010), (1110, 1120)]),
])
def test_reference_blocks_split_only_at_skipped_introns(cigar, blocks):
    assert _reference_blocks(1000, cigar) == blocks


S_REF, G_REF = REFERENCE[10], REFERENCE[30]


@pytest.mark.parametrize("reads,expected", [
    # Heterozygous germline alt in cis at 40% purity: normal cells add many (somatic ref, germline alt).
    ([read("t%d" % i, "G", "C") for i in range(4)] + [read("n%d" % i, S_REF, "C") for i in range(6)]
     + [read("b%d" % i, S_REF, G_REF) for i in range(10)], True),
    # The same sample with the somatic variant on the other copy.
    ([read("t%d" % i, "G", G_REF) for i in range(4)] + [read("n%d" % i, S_REF, "C") for i in range(10)]
     + [read("b%d" % i, S_REF, G_REF) for i in range(6)], False),
    # Homozygous germline alt: every fragment has it, so the somatic variant is cis.
    ([read("a%d" % i, "G", "C") for i in range(3)] + [read("r%d" % i, S_REF, "C") for i in range(7)], True),
    # The canonical scenario from test_read_phasing, with the germline variant supplied as germline.
    ([read("c%d" % i, "G", "C") for i in range(3)] + [read("s%d" % i, "G", G_REF) for i in range(2)]
     + [read("g%d" % i, S_REF, "C") for i in range(2)], True),
])
def test_only_somatic_alt_fragments_decide(tmp_path, reads, expected):
    # Fragments with the somatic reference allele say nothing about its copy (#391 review).
    result, = results_for(tmp_path, reads)
    assert IsovarReadPhasing([result]).in_cis(SOMATIC, GERMLINE) is expected


def test_without_somatic_alt_reads_nothing_is_examined(tmp_path):
    reads = [read("r%d" % i, S_REF, "C") for i in range(3)] + [read("q%d" % i, S_REF, G_REF) for i in range(3)]
    result, = results_for(tmp_path, reads)
    assert result.germline_read_evidence == {}
    assert IsovarReadPhasing([result]).in_cis(SOMATIC, GERMLINE) is None


def test_kept_germline_reads_are_those_on_somatic_alt_fragments(tmp_path):
    reads = [read("t%d" % i, "G", "C") for i in range(2)] + [read("n%d" % i, S_REF, "C") for i in range(5)]
    result, = results_for(tmp_path, reads)
    assert result.germline_read_evidence[GERMLINE].alt_read_names == {"t0", "t1"}


def test_a_germline_variant_in_the_next_exon_is_examined(tmp_path):
    # 15 aligned bases, a 10-base skipped intron, then 15 more over the germline locus.
    sequence = REFERENCE[:10] + "G" + REFERENCE[11:15] + REFERENCE[25:30] + "C" + REFERENCE[31:40]
    reads = [record("x%d" % i, "1", 1000, "15M10N15M", sequence) for i in range(3)]
    result, = results_for(tmp_path, reads)
    assert list(result.germline_read_evidence) == [GERMLINE]
    assert IsovarReadPhasing([result]).in_cis(SOMATIC, GERMLINE) is True


def test_germline_deletions_and_prefixed_contigs_are_examined(tmp_path):
    deletion = Variant("chr1", 1031, REFERENCE[30:32], REFERENCE[30], normalize_contig_names=False)
    sequence = REFERENCE[:10] + "G" + REFERENCE[11:31] + REFERENCE[32:40]
    reads = [record("d%d" % i, "1", 1000, "31M1D8M", sequence) for i in range(3)]
    result, = results_for(tmp_path, reads, germline=(deletion,))
    assert list(result.germline_read_evidence) == [deletion]
    assert IsovarReadPhasing([result]).in_cis(SOMATIC, deletion) is True


def test_a_shared_germline_variant_is_collected_once(tmp_path, monkeypatch):
    second = Variant("1", 1016, REFERENCE[15], "T" if REFERENCE[15] != "T" else "A", normalize_contig_names=False)
    bases = list(REFERENCE[:40])
    bases[10], bases[15] = "G", second.alt
    path = write_bam(tmp_path / "rna.bam", [record("m%d" % i, "1", 1000, "40M", "".join(bases)) for i in range(2)])
    calls = []
    collect = ReadCollector.read_evidence_for_variant
    monkeypatch.setattr(ReadCollector, "read_evidence_for_variant",
                        lambda self, variant, alignment_file: calls.append(variant) or collect(
                            self, variant, alignment_file))
    with pysam.AlignmentFile(str(path)) as bam:
        results = run_isovar([SOMATIC, second], bam, germline_variants=[GERMLINE])
    assert all(list(r.germline_read_evidence) == [GERMLINE] for r in results)
    assert calls.count(GERMLINE) == 1


def _legacy(variant, alt, ref=()):
    from isovar.allele_read import AlleleRead
    from isovar.read_evidence import ReadEvidence
    return ReadEvidence(
        trimmed_base1_start=variant.start, trimmed_ref=variant.ref, trimmed_alt=variant.alt,
        ref_reads=[AlleleRead(prefix="A", allele=variant.ref, suffix="T", name=n) for n in ref],
        alt_reads=[AlleleRead(prefix="A", allele=variant.alt, suffix="T", name=n) for n in alt],
        other_reads=[])


@pytest.mark.parametrize("germline_alt,germline_ref,in_protein,expected", [
    (set(), {"f1", "f2"}, False, False),
    (set(), {"f1", "f2"}, True, None),   # the reads and the assembly disagree
    (set(), set(), True, True),          # no decisive reads: the assembly decides
    ({"f1", "f2"}, set(), True, True),
])
def test_reads_decide_before_the_assembled_protein(germline_alt, germline_ref, in_protein, expected):
    result = IsovarResult(variant=SOMATIC, read_evidence=_legacy(SOMATIC, {"f1", "f2"}), predicted_effect=None,
                          germline_read_evidence={GERMLINE: _legacy(GERMLINE, germline_alt, germline_ref)})
    from types import SimpleNamespace
    shaped = SimpleNamespace(
        variant=SOMATIC, num_alt_fragments=2, phased_variants_in_supporting_reads=set(),
        alt_reads=result.alt_reads, ref_reads=result.ref_reads,
        germline_read_evidence=result.germline_read_evidence,
        germline_variants_in_top_protein_sequence={GERMLINE} if in_protein else set())
    assert IsovarReadPhasing([shaped]).in_cis(SOMATIC, GERMLINE) is expected


def test_results_pickled_before_the_field_existed_still_work():
    import pickle
    result = IsovarResult(variant=SOMATIC, read_evidence=_legacy(SOMATIC, {"f1"}), predicted_effect=None)
    del result.__dict__["germline_read_evidence"]
    restored = pickle.loads(pickle.dumps(result))
    assert restored.to_dict()["germline_read_evidence"] == {} and restored == restored.clone()
    shaped_none = IsovarReadPhasing([restored])
    assert shaped_none.in_cis(SOMATIC, GERMLINE) is None
