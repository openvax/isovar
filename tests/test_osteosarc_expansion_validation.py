"""Offline tests of audit coordinate, protein, provenance and failure contracts."""

from copy import deepcopy
import json
import signal
import time
from types import SimpleNamespace

import pytest

from tests.data.osteosarc.expansion.acquire import assembly_from_header, regions_for_header
from tests.data.osteosarc.expansion.liftover import chain_blocks, lift_variant, map_interval
from tests.data.osteosarc.expansion.metadata import PageSections, cellranger_config
from tests.data.osteosarc.expansion.references import apply_variant, check_translation, minimal_edit, translate
from tests.data.osteosarc.expansion.runner import AuditTimeout, audit_mode, counts_from_evidence, time_limit


def toy_model(strand="+"):
    return dict(contig="1", strand=strand, exons=[(100, 117)], cds_start=0, cds_end=15,
                cdna="ATGAAACCCGGGTTTTAA", protein="MKPGF", genetic_code=1,
                transcript_version="transcript.1", protein_version="protein.1")


@pytest.mark.parametrize("strand,pos,ref,alt", [("+", 103, "A", "G"), ("-", 114, "T", "C")])
def test_independent_plus_and_minus_strand_edit(strand, pos, ref, alt):
    result = apply_variant(dict(chrom="chr1", pos=pos, ref=ref, alt=alt), toy_model(strand))
    assert result["mutant_cdna"] == "ATGGAACCCGGGTTTTAA"
    assert result["mutant_protein"] == "MEPGF"
    assert result["mutant_has_stop"] and not result["frameshift"]
    assert result["variant_offset"] == 3


@pytest.mark.parametrize("ref,alt,protein,frameshift", [
    ("GAAA", "G", "MPGF", False), ("G", "GAAA", "MKKPGF", False),
    ("G", "GT", "M", True), ("GAAA", "GAGG", "MRPGF", False),
])
def test_independent_insertion_deletion_complex_and_stop(ref, alt, protein, frameshift):
    result = apply_variant(dict(chrom="chr1", pos=102, ref=ref, alt=alt), toy_model())
    assert result["mutant_protein"] == protein
    assert result["frameshift"] is frameshift


@pytest.mark.parametrize("record", [dict(chrom="2", pos=103, ref="A", alt="G"),
                                     dict(chrom="1", pos=103, ref="T", alt="G"),
                                     dict(chrom="1", pos=200, ref="A", alt="G"),
                                     dict(chrom="1", pos=115, ref="T", alt="G")])
def test_reference_mismatch_or_unavailable_interval_is_not_silently_edited(record):
    with pytest.raises(ValueError):
        apply_variant(record, toy_model())


def toy_translation():
    expected = apply_variant(dict(chrom="1", pos=103, ref="A", alt="G"), toy_model())
    orf = SimpleNamespace(cdna_sequence="ATGGAACCCGGGTTTTAA", variant_cdna_interval_start=3,
                          variant_cdna_interval_end=4, reference_cdna_sequence_before_variant="ATG",
                          offset_to_first_complete_codon=0)
    translation = SimpleNamespace(variant_orf=orf, frameshift=False, amino_acids="MEPGF",
                                  ends_with_stop_codon=True, mutation_start_idx=1, mutation_end_idx=2,
                                  contains_mutation=True)
    return translation, expected


@pytest.mark.parametrize("length", [0, 49])
def test_independent_oracle_validates_frame_interval_and_stop(length):
    observed, expected = toy_translation()
    result = check_translation(observed, expected, length)
    assert result["matches_expected"] and result["expected_stop"]
    assert result["protein_start_1based"] == 1


@pytest.mark.parametrize("attribute,value", [("amino_acids", "MXXXX"), ("ends_with_stop_codon", False),
                                             ("mutation_start_idx", 0), ("mutation_end_idx", 3),
                                             ("frameshift", True), ("contains_mutation", False)])
def test_independent_oracle_rejects_corrupted_translation(attribute, value):
    observed, expected = toy_translation()
    setattr(observed, attribute, value)
    with pytest.raises(AssertionError):
        check_translation(observed, expected, 49)


def test_independent_oracle_separates_additional_rna_difference_from_translation_error():
    observed, expected = toy_translation()
    observed.variant_orf.cdna_sequence = "ATGGAAGCCGGGTTTTAA"
    observed.amino_acids = "MEAGF"
    result = check_translation(observed, expected, 49)
    assert result["matches_expected"] is False
    assert result["expected_amino_acids"] == "MEPGF"


def test_independent_oracle_caps_sequence_and_clears_stop():
    observed, expected = toy_translation()
    observed.amino_acids, observed.ends_with_stop_codon = "MEP", False
    assert check_translation(observed, expected, 3)["matches_expected"]


def test_reference_mitochondrial_code_changes_translation():
    assert translate("TGAATA", table=2) == ("WM", False)
    assert translate("TGAATA", table=1) == ("", True)
    assert translate("ATGAGAGGG", table=2) == ("M", True)


def test_minimal_edit_preserves_locus_not_repeat_normalization():
    assert minimal_edit(dict(pos=100, ref="AGAA", alt="AGA")) == (103, "A", "")
    with pytest.raises(ValueError, match="Identical"):
        minimal_edit(dict(pos=100, ref="A", alt="A"))


@pytest.mark.parametrize("strand,expected", [("+", (22, 25)), ("-", (75, 78))])
def test_chain_forward_and_reverse_half_open_coordinates(strand, expected):
    blocks = list(chain_blocks([f"chain 100 chr1 100 + 10 20 chr2 100 {strand} 20 30 1\n", "10\n"]))
    mapped = map_interval(blocks, "chr1", 12, 15)
    assert (mapped["start"], mapped["end"]) == expected
    assert mapped["strand"] == strand


def test_chain_rejects_ambiguous_gapped_and_truncated_mapping():
    lines = ["chain 100 chr1 100 + 10 25 chr2 100 + 20 35 1\n", "5 5 5\n", "5\n"]
    blocks = list(chain_blocks(lines))
    with pytest.raises(ValueError, match="found 0"):
        map_interval(blocks, "chr1", 14, 21)
    with pytest.raises(ValueError, match="found 2"):
        map_interval(blocks + blocks, "chr1", 10, 11)
    with pytest.raises(ValueError, match="Truncated"):
        list(chain_blocks(lines[:-1]))


def test_lift_reverse_substitution_complements_alleles():
    blocks = list(chain_blocks(["chain 100 chr1 100 + 10 20 chr2 100 - 20 30 1\n", "10\n"]))
    record = dict(assembly="GRCh38", chrom="chr1", pos=13, ref="AC", alt="GT", variant_id="example")
    mapped = lift_variant(record, blocks)
    assert (mapped["chrom"], mapped["pos"], mapped["ref"], mapped["alt"]) == ("chr2", 77, "GT", "AC")
    with pytest.raises(ValueError, match="reanchoring"):
        lift_variant(dict(record, alt="A"), blocks)


def test_assembly_detection_requires_multiple_consistent_contigs():
    assert assembly_from_header({}) is None
    assert assembly_from_header({"SQ": [dict(SN="chr1", LN=248956422)]}) is None
    header = {"SQ": [dict(SN="chr1", LN=248956422), dict(SN="chr2", LN=242193529)]}
    assert assembly_from_header(header) == "GRCh38"
    header["SQ"].append(dict(SN="3", LN=198022430))
    assert assembly_from_header(header) is None


def test_region_aliasing_and_explicit_missing_or_ambiguous_contigs():
    variants = [dict(chrom="chrM", pos=12994, ref="G", alt="A", variant_id="MT-ND5"),
                dict(chrom="chr1", pos=1, ref="AC", alt="A", variant_id="deletion")]
    regions, missing = regions_for_header(variants, {"SQ": [dict(SN="MT")]})
    assert regions == ["MT:12993-12995"] and missing == ["deletion"]
    with pytest.raises(ValueError, match="Ambiguous"):
        regions_for_header(variants, {"SQ": [dict(SN="MT"), dict(SN="chrM")]})


def test_source_metadata_keeps_compound_definition_and_membership():
    parser = PageSections()
    parser.feed("<h2>Genomic context</h2><p>cDNA change AG&gt;GT</p>"
                "<h2>Vaccines targeting this variant</h2><b>Cure 2024</b>"
                "<script>unrelated script</script><h2>Contact</h2>unrelated")
    assert parser.result() == {"Genomic context": "cDNA change AG>GT", "Vaccines targeting this variant": "Cure 2024"}


def test_library_config_does_not_invent_donor_assignments():
    result = cellranger_config("[gene-expression]\nreference,GRCh38\ntenx-cloud-token-path,/private/token\n"
                               "[libraries]\nfastq_id,fastqs\nRNA-1,/source/reads\n")
    assert "samples" not in result
    assert result["gene-expression"] == [["reference", "GRCh38"]]
    assert result["libraries"][1] == ["RNA-1", "/source/reads"]


def test_alignment_read_objects_and_template_names_are_distinct_units():
    evidence = SimpleNamespace(ref_reads=[], alt_reads=[SimpleNamespace(name="pair", source_read_count=2),
                                                       SimpleNamespace(name="pair", source_read_count=1)], other_reads=[])
    counts = counts_from_evidence(evidence)
    assert counts["reads"]["alt"] == 3
    assert counts["read_objects"]["alt"] == 2
    assert counts["template_names"]["alt"] == 1
    assert counts_from_evidence(None) is None


def test_pipeline_exception_is_not_zero_support(monkeypatch):
    from tests.data.osteosarc.expansion import runner

    def fail(*args, **kwargs):
        raise OSError("unreadable original source")

    monkeypatch.setattr(runner.pysam, "AlignmentFile", fail)
    result = audit_mode("missing.bam", None, {})
    assert result["status"] == "error" and result["counts"] is None
    assert result["error"]["type"] == "OSError"
    assert "proteins" not in result


@pytest.mark.parametrize("seconds", [0, -1, float("inf"), float("nan")])
def test_invalid_audit_budget(seconds):
    with pytest.raises(ValueError, match="finite and positive"):
        with time_limit(seconds):
            pytest.fail("Invalid limit entered the pipeline")


@pytest.mark.skipif(not hasattr(signal, "setitimer"), reason="Unix audit runner")
def test_audit_budget_interrupts_and_restores_signal_handler():
    previous = signal.getsignal(signal.SIGALRM)
    with pytest.raises(AuditTimeout):
        with time_limit(0.01):
            time.sleep(1)
    assert signal.getsignal(signal.SIGALRM) == previous
    assert signal.getitimer(signal.ITIMER_REAL) == (0.0, 0.0)
    with time_limit(1):
        pass
    assert signal.getsignal(signal.SIGALRM) == previous


@pytest.mark.skipif(not hasattr(signal, "setitimer"), reason="Unix audit runner")
def test_audit_does_not_replace_existing_timer():
    previous = signal.getsignal(signal.SIGALRM)
    signal.setitimer(signal.ITIMER_REAL, 100)
    try:
        with pytest.raises(ValueError, match="existing alarm"):
            with time_limit(1):
                pytest.fail("Existing timer was overwritten")
        assert signal.getitimer(signal.ITIMER_REAL)[0] > 90
        assert signal.getsignal(signal.SIGALRM) == previous
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


@pytest.mark.parametrize("after_collection", [False, True])
def test_audit_timeout_keeps_completed_counts_and_stage(monkeypatch, after_collection):
    from contextlib import nullcontext
    from tests.data.osteosarc.expansion import runner

    def fail(*args, read_collector, **kwargs):
        read_collector.trace["stage"] = "rna_sequence_creation" if after_collection else "read_collection"
        if after_collection:
            read_collector.evidence = SimpleNamespace(
                ref_reads=[], alt_reads=[SimpleNamespace(name="real-read", source_read_count=2)], other_reads=[])
        raise AuditTimeout("budget expired")

    monkeypatch.setattr(runner.pysam, "AlignmentFile", lambda path: nullcontext())
    monkeypatch.setattr(runner, "run_isovar", fail)
    result = audit_mode("source.bam", None, {})
    assert result["status"] == "resource_limit" and result["outcome"] == "audit_timeout"
    assert result["error"]["type"] == "AuditTimeout"
    assert result["stages"]["stage"] == ("rna_sequence_creation" if after_collection else "read_collection")
    assert result["counts"] == (dict(reads=dict(ref=0, alt=2, other=0), read_objects=dict(ref=0, alt=1, other=0),
                                      template_names=dict(ref=0, alt=1, other=0), all_template_names=1)
                                if after_collection else None)
    assert "proteins" not in result


def test_uncapped_capture_restores_public_limit_on_exception(monkeypatch):
    from tests.data.osteosarc.expansion import runner

    def fail(creator, *args, **kwargs):
        assert creator.max_protein_sequences_per_variant is None
        raise AuditTimeout("interrupted ranking")

    monkeypatch.setattr(runner.ProteinSequenceCreator, "sorted_protein_sequences_for_variant", fail)
    creator = runner.TracedCreator({}, capture_all_ranked=True)
    with pytest.raises(AuditTimeout):
        creator.sorted_protein_sequences_for_variant(None, None)
    assert creator.max_protein_sequences_per_variant == 1


@pytest.mark.parametrize("cigar,position,expected", [
    ("20M", 109, (100, 120)),
    ("10M10000N10M", 10111, (10110, 10120)),
    ("10M10000N10M", 115, None),
    ("10M1I10M", 111, (100, 120)),
    ("10M2D10M", 115, (100, 122)),
    ("10M51D10M", 162, (161, 171)),
    ("2S10M2S", 105, (100, 110)),
    ("2H10M2H", 105, (100, 110)),
])
def test_numt_diagnostic_does_not_count_an_intron_as_contiguous_mt_context(cigar, position, expected):
    import pysam
    from tests.data.osteosarc.expansion.mitochondrial import variant_block

    read = pysam.AlignedSegment()
    read.reference_start, read.cigarstring = 100, cigar
    assert variant_block(read, position) == expected
