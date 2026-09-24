"""The protein hypothesis export keeps every alternative and scoped RNA evidence."""

import csv
from dataclasses import dataclass
import json

from pysam import AlignmentFile
import pytest
from varcode import Variant

from isovar import (
    ProteinSequence, ProteinSequenceCreator, export_protein_hypotheses, read_group_metadata,
    run_isovar, union_rna_support, write_protein_hypotheses,
)
from isovar.allele_read import AlleleRead
from isovar.cli.commands import run as isovar_cli
from isovar.isovar_result import IsovarResult
from isovar.read_evidence import ReadEvidence
from isovar.reference_context import ReferenceContext
from isovar.translation import Translation
from isovar.variant_orf import VariantORF
from isovar.variant_sequence import VariantSequence

from .testing_helpers import data_path

BAM = data_path("data/b16.f10/b16.combined.sorted.bam")
VCF = data_path("data/b16.f10/b16.vcf")


@dataclass(frozen=True, order=True)
class Gene:
    id: str
    name: str


@dataclass(frozen=True, order=True)
class Transcript:
    id: str
    name: str
    gene: Gene
    strand: str = "+"
    genomic_offset: int = 100

    def spliced_offset(self, dna_pos):
        return dna_pos - self.genomic_offset


VARIANT = Variant("1", 108, "A", "C", normalize_contig_names=False)
TRANSCRIPT = Transcript(id="tx-1", name="TX1", gene=Gene(id="gene-1", name="GENE1"))
# Reference cDNA CC|ATG GCT AAA GGC TAA: 5' UTR, start codon, A>C makes K>Q.
CONTEXT = ReferenceContext(
    strand="+", sequence_before_variant_locus="CCATGGCT", sequence_at_variant_locus="A",
    sequence_after_variant_locus="AAGGCTAA", offset_to_first_complete_codon=2, contains_start_codon=True,
    overlaps_start_codon=True, contains_five_prime_utr=True, amino_acids_before_variant="MA",
    variant=VARIANT, transcripts=[TRANSCRIPT])


def read(name, mate=0, group="rg1", with_identity=True):
    alignments = (((group, name, mate), (0, 100, "17M", False)),) if with_identity else ()
    return AlleleRead(prefix="GCT", allele="C", suffix="AAG", name=name, source_alignments=alignments)


def translation(cdna, offset, variant_start, reference_before, reference_after, amino_acids, stop,
                mutation, reads, mismatches=(0, 0)):
    orf = VariantORF(
        cdna_sequence=cdna, offset_to_first_complete_codon=offset,
        variant_cdna_interval_start=variant_start, variant_cdna_interval_end=variant_start + 1,
        reference_cdna_sequence_before_variant=reference_before,
        reference_cdna_sequence_after_variant=reference_after,
        num_mismatches_before_variant=mismatches[0], num_mismatches_after_variant=mismatches[1])
    sequence = VariantSequence(cdna[:variant_start], "C", cdna[variant_start + 1:], reads)
    return Translation(
        amino_acids=amino_acids, contains_mutation=True, mutation_start_idx=mutation[0],
        mutation_end_idx=mutation[1], ends_with_stop_codon=stop, frameshift=False,
        untrimmed_variant_sequence=sequence, reference_context=CONTEXT, variant_orf=orf)


READS = dict(a1=read("a", 64), a2=read("a", 128), b=read("b"), c=read("c"), d=read("d", group="rg2"))


def synthetic_result(settings=None, reads=READS):
    r = reads
    full = translation("CCATGGCTCAAGGCTAA", 2, 8, "CCATGGCT", "AAGGCTAA", "MAQG", True, (2, 3),
                       [r["a1"], r["a2"], r["b"]])
    # A synonymous GCT>GCC change: another cDNA for the same protein.
    synonymous = translation("CCATGGCCCAAGGCTAA", 2, 8, "CCATGGCT", "AAGGCTAA", "MAQG", True, (2, 3),
                             [r["c"]], mismatches=(1, 0))
    # An adjacent GGC>GAC change gives a different protein.
    alternative = translation("CCATGGCTCAAGACTAA", 2, 8, "CCATGGCT", "AAGGCTAA", "MAQD", True, (2, 3),
                              [r["b"], r["d"]], mismatches=(0, 1))
    # A shorter local window, without the start or stop codon.
    window = translation("GCTCAAGGC", 0, 3, "GCT", "AAGGC", "AQG", False, (1, 2), [r["a1"]])
    proteins = [ProteinSequence.from_translations(ts) for ts in ([full, synonymous], [alternative], [window])]
    evidence = ReadEvidence(trimmed_base1_start=108, trimmed_ref="A", trimmed_alt="C", ref_reads=[],
                            alt_reads=list(reads.values()), other_reads=[])
    return IsovarResult(variant=VARIANT, read_evidence=evidence, predicted_effect=None,
                        sorted_protein_sequences=proteins, protein_sequence_settings=settings)


def evidence(export, support):
    return export["evidence_sets"][support["evidence_set_id"]]


def test_synonymous_translations_share_one_protein_with_distinct_nucleotide_provenance():
    export = export_protein_hypotheses([synthetic_result()], sample_id="s", source="rna.bam")
    event, = export["events"]
    first, alternative, window = event["protein_hypotheses"]
    assert [p["isovar_rank"] for p in event["protein_hypotheses"]] == [1, 2, 3]
    assert len(first["translations"]) == 2
    assert len({t["nucleotide_sequence_id"] for t in first["translations"]}) == 2
    assert len({t["translation_id"] for t in first["translations"]}) == 2
    synonymous, = [t for t in first["translations"] if t["mismatches_before_variant"] == 1]
    assert {(e["cdna_interval"][0], e["alt_bases"], e["origin"]) for e in synonymous["observed_edits"]} == {
        (8, "C", "nominated_variant"), (7, "C", "unexplained")}
    # Each translation keeps its own reads; the protein counts their union once.
    assert [(t["rna_support"]["segments"], t["rna_support"]["fragments"]) for t in sorted(
        first["translations"], key=lambda t: t["mismatches_before_variant"])] == [(3, 2), (1, 1)]
    assert (first["rna_support"]["segments"], first["rna_support"]["fragments"]) == (4, 3)
    assert (alternative["rna_support"]["segments"], alternative["rna_support"]["fragments"]) == (2, 2)
    assert first["protein_sequence_id"] != alternative["protein_sequence_id"]
    assert first["mutation_interval"] == [2, 3] and window["mutation_interval"] == [1, 2]
    assert json.loads(json.dumps(export)) == export


def test_termini_distinguish_stop_complete_and_partial_windows():
    first, _, window = export_protein_hypotheses(
        [synthetic_result()], sample_id="s", source="rna.bam")["events"][0]["protein_hypotheses"]
    assert (first["n_terminus"], first["c_terminus"]) == ("annotated_start_codon", "stop_codon")
    assert (window["n_terminus"], window["c_terminus"]) == ("local_window", "no_stop_in_context")
    assert all(t["starts_at_annotated_start_codon"] for t in first["translations"])
    assert first["translations"][0]["translated_interval"] == [2, 17]


def test_shorter_windows_are_marked_covered_without_losing_them():
    first, alternative, window = export_protein_hypotheses(
        [synthetic_result()], sample_id="s", source="rna.bam")["events"][0]["protein_hypotheses"]
    assert (first["representative"], alternative["representative"], window["representative"]) == (
        True, True, False)
    covering, = window["covered_by"]
    assert covering["hypothesis_id"] == first["hypothesis_id"]
    start, end = covering["protein_interval"]
    assert first["amino_acids"][start:end] == window["amino_acids"]
    assert window["translations"] and window["rna_support"]["segments"] == 1


def test_combining_hypotheses_counts_shared_reads_once():
    export = export_protein_hypotheses([synthetic_result()], sample_id="s", source="rna.bam")
    event, = export["events"]
    first, alternative, _ = event["protein_hypotheses"]
    combined = union_rna_support([evidence(export, first["rna_support"]),
                                  evidence(export, alternative["rna_support"])])
    # Read b supports both proteins; adding the counts would give 6 and 5.
    assert (combined["segments"], combined["fragments"]) == (5, 4)
    assert combined["evidence_set_id"] == event["allele_support"]["alt"]["evidence_set_id"]


def test_exports_of_the_same_reads_share_ids_and_other_scopes_cannot_combine():
    one = export_protein_hypotheses([synthetic_result()], sample_id="s", source="rna.bam")
    two = export_protein_hypotheses([synthetic_result()], sample_id="s", source="rna.bam")
    other = export_protein_hypotheses([synthetic_result()], sample_id="s", source="reprocessed.bam")
    support = one["events"][0]["protein_hypotheses"][0]["rna_support"]
    assert support["evidence_set_id"] in two["evidence_sets"]
    assert support["evidence_set_id"] not in other["evidence_sets"]
    same = union_rna_support([evidence(one, support), evidence(two, support)])
    assert same["segments"] == support["segments"]
    with pytest.raises(ValueError, match="different evidence scopes"):
        union_rna_support([evidence(one, support), next(iter(other["evidence_sets"].values()))])
    with pytest.raises(ValueError, match="No RNA supports"):
        union_rna_support([])
    with pytest.raises(ValueError, match="evidence_sets"):
        union_rna_support([support])


def test_unknown_library_metadata_and_reads_without_identity_stay_unknown():
    export = export_protein_hypotheses([synthetic_result()], sample_id="s", source="rna.bam")
    assert export["read_group_metadata"] == "not_supplied"
    assert export["read_groups"] == {"rg1": None, "rg2": None}
    header = {"RG": [dict(ID="rg1", SM="tumor", LB="lib1"), dict(ID="dup"), dict(ID="dup")]}
    assert read_group_metadata(header) == {"rg1": dict(sample="tumor", library="lib1", platform=None)}
    known = export_protein_hypotheses([synthetic_result()], sample_id="s", source="rna.bam",
                                      alignment_header=header)
    assert known["read_groups"] == {"rg1": dict(sample="tumor", library="lib1", platform=None), "rg2": None}
    legacy = {name: read(name, with_identity=False) for name in ("a1", "a2", "b", "c", "d")}
    unscoped = export_protein_hypotheses([synthetic_result(reads=legacy)], sample_id="s", source="rna.bam")
    support = unscoped["events"][0]["protein_hypotheses"][0]["rna_support"]
    assert support["evidence_set_id"] is None and support["segments"] == 4
    assert unscoped["evidence_sets"] == {}


def test_completeness_follows_the_recorded_protein_limit():
    def complete(settings):
        event, = export_protein_hypotheses([synthetic_result(settings)], sample_id="s",
                                           source="rna.bam")["events"]
        return event["protein_sequence_limit"], event["protein_hypotheses_complete"]
    assert complete(None) == (None, None)
    assert complete(dict(max_protein_sequences_per_variant=0)) == (0, True)
    assert complete(dict(max_protein_sequences_per_variant=5)) == (5, True)
    assert complete(dict(max_protein_sequences_per_variant=3)) == (3, False)


def test_inconsistent_translation_is_rejected():
    result = synthetic_result()
    result.sorted_protein_sequences[0].translations[0].variant_orf.offset_to_first_complete_codon = 1
    with pytest.raises(ValueError, match="disagree"):
        export_protein_hypotheses([result], sample_id="s", source="rna.bam")


def test_real_rna_export_keeps_every_hypothesis_and_records_settings(tmp_path):
    creator = ProteinSequenceCreator(max_protein_sequences_per_variant=0)
    results = run_isovar(VCF, BAM, protein_sequence_creator=creator)
    assert all(r.protein_sequence_settings == creator.settings() for r in results)
    export = export_protein_hypotheses(results, sample_id="b16", source="b16.bam",
                                       alignment_header=AlignmentFile(BAM).header)
    assert export["schema"] == "isovar.protein_hypotheses.v1"
    assert export["read_groups"]["Tumor_B16_F10_0810_127A"]["sample"] == "Tumor_B16_F10_0810_127A"
    events = [e for e in export["events"] if e["protein_hypotheses"]]
    assert events and all(e["protein_hypotheses_complete"] for e in export["events"])
    for event, result in zip(export["events"], results):
        assert len(event["protein_hypotheses"]) == len(result.sorted_protein_sequences)
        alt = event["allele_support"]["alt"]
        assert (alt["segments"], alt["fragments"]) == (result.num_alt_reads, result.num_alt_fragments)
        if event["protein_hypotheses"]:
            assert any(p["representative"] for p in event["protein_hypotheses"])
            union = union_rna_support([evidence(export, p["rna_support"]) for p in event["protein_hypotheses"]])
            assert union["segments"] <= alt["segments"]
    paths = write_protein_hypotheses(export, tmp_path / "hypotheses.json")
    assert json.loads(paths["json"].read_text()) == export
    with paths["tsv"].open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == sum(len(p["translations"]) for e in events for p in e["protein_hypotheses"])


def test_two_reconstruction_modes_of_the_same_reads_combine_without_inflation():
    exports = [export_protein_hypotheses(
        run_isovar(VCF, BAM, protein_sequence_creator=ProteinSequenceCreator(
            max_protein_sequences_per_variant=0, variant_sequence_assembly=assembly)),
        sample_id="b16", source="b16.bam") for assembly in (True, False)]
    for events in zip(*(e["events"] for e in exports)):
        supports = [evidence(export, p["rna_support"]) for export, event in zip(exports, events)
                    for p in event["protein_hypotheses"]]
        if supports:
            combined = union_rna_support(supports)
            assert combined["segments"] <= events[0]["allele_support"]["alt"]["segments"]
            assert combined["segments"] < sum(s["segments"] for s in supports)


def test_write_rejects_other_schemas(tmp_path):
    with pytest.raises(ValueError, match="protein_hypotheses.v1"):
        write_protein_hypotheses({"schema": "isovar.sv_rna_orfs.v3"}, tmp_path / "x.json")


def test_command_keeps_all_proteins_by_default(tmp_path):
    output = tmp_path / "b16.json"
    isovar_cli(["protein-hypotheses", "--vcf", VCF, "--bam", BAM, "--sample-id", "b16",
                "--output", str(output), "--log-level", "WARNING"])
    export = json.loads(output.read_text())
    assert export["source"] == BAM and (tmp_path / "b16.tsv").exists()
    assert {e["protein_sequence_limit"] for e in export["events"]} == {0}
    assert max(len(e["protein_hypotheses"]) for e in export["events"]) > 1


def test_sv_and_small_variant_exports_of_one_read_set_share_read_ids():
    from isovar import export_sv_rna_orfs
    from tests.test_sv_rna_orf_export import reconstruction

    sv = export_sv_rna_orfs(reconstruction())["candidates"][0]["rna_support"]
    assert sv["evidence_scope"] == ["sample", "source"]
    shared = dict(READS, b=read("full", group="library"))
    small = export_protein_hypotheses([synthetic_result(reads=shared)], sample_id="sample", source="source")
    support = evidence(small, small["events"][0]["protein_hypotheses"][0]["rna_support"])
    assert set(sv["segment_ids"]) <= set(support["segment_ids"])
    assert union_rna_support([sv, support])["segments"] == support["segments"]


def test_command_reports_unusable_output_as_a_usage_error(tmp_path, capsys):
    with pytest.raises(SystemExit) as exit_info:
        isovar_cli(["protein-hypotheses", "--vcf", VCF, "--bam", BAM, "--sample-id", "b16",
                    "--output", str(tmp_path / "missing" / "b16.json")])
    assert exit_info.value.code == 2
    assert "Output directory does not exist" in capsys.readouterr().err
