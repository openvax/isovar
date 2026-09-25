"""Assembled cDNA edits are attributed to supplied somatic and germline variants (#297)."""

from dataclasses import dataclass
import json

import pysam
import pytest
from varcode import Variant

from isovar import ProteinSequence, ReadCollector, export_protein_hypotheses, run_isovar
from isovar.known_variants import KnownVariants
from isovar.reference_context import ReferenceContext
from isovar.transcript_edit_helpers import categorize_transcript_assembly_edits_from_translation
from isovar.translation import Translation
from isovar.variant_orf import VariantORF
from isovar.variant_sequence import VariantSequence


@dataclass(frozen=True)
class Gene:
    id: str
    name: str


@dataclass(frozen=True)
class Transcript:
    id: str
    name: str
    gene: Gene
    contig: str = "1"
    start: int = 100
    end: int = 200
    strand: str = "+"

    def spliced_offset(self, position):
        if not self.start <= position <= self.end:
            raise ValueError("Invalid position: %d" % position)
        return position - self.start


def variant(position, ref, alt):
    return Variant("1", position, ref, alt, normalize_contig_names=False)


TRANSCRIPT = Transcript("tx-1", "TX1", Gene("gene-1", "GENE1"))
# Transcript offsets 0-16: CC|ATG GCT AAA GGC TAA (M A K G *); A>C at offset 8 gives K>Q.
FOCAL = variant(108, "A", "C")
SILENT_GERMLINE = variant(107, "T", "C")        # GCT>GCC keeps Ala
CHANGING_GERMLINE = variant(112, "G", "A")      # GGC>GAC gives Asp
CO_SOMATIC = variant(113, "C", "T")             # adjacent to the germline change
UNSUPPLIED = variant(100, "C", "G")             # a 5' UTR change nobody supplied


def translation(prefix, suffix, amino_acids, reference_prefix="CCATGGCT", reference_suffix="AAGGCTAA",
                focal=FOCAL, context_before="CCATGGCT", context_after="AAGGCTAA"):
    context = ReferenceContext(
        strand="+", sequence_before_variant_locus=context_before, sequence_at_variant_locus=focal.ref,
        sequence_after_variant_locus=context_after, offset_to_first_complete_codon=2,
        contains_start_codon=True, overlaps_start_codon=True, contains_five_prime_utr=True,
        amino_acids_before_variant="MA", variant=focal, transcripts=[TRANSCRIPT])
    cdna = prefix + focal.alt + suffix
    orf = VariantORF(
        cdna_sequence=cdna, offset_to_first_complete_codon=2, variant_cdna_interval_start=len(prefix),
        variant_cdna_interval_end=len(prefix) + len(focal.alt),
        reference_cdna_sequence_before_variant=reference_prefix,
        reference_cdna_sequence_after_variant=reference_suffix,
        num_mismatches_before_variant=sum(a != b for a, b in zip(prefix, reference_prefix)),
        num_mismatches_after_variant=sum(a != b for a, b in zip(suffix, reference_suffix)))
    return Translation(
        amino_acids=amino_acids, contains_mutation=True, mutation_start_idx=2, mutation_end_idx=3,
        ends_with_stop_codon=True, frameshift=False,
        untrimmed_variant_sequence=VariantSequence(prefix, focal.alt, suffix, []),
        reference_context=context, variant_orf=orf)


# The observed cDNA carries every change: GCC|CAA|GAT, and a UTR G.
OBSERVED = translation("GCATGGCC", "AAGATTAA", "MAQD")


def categories(known, t=OBSERVED):
    edits = categorize_transcript_assembly_edits_from_translation(t, TRANSCRIPT, known)
    return {category: {(e.cdna_start, e.cdna_end, e.alt_bases, e.source_variant) for e in group}
            for category, group in edits.items()}


def test_germline_somatic_and_unmatched_edits_are_separated():
    known = KnownVariants(somatic=[FOCAL, CO_SOMATIC], germline=[SILENT_GERMLINE, CHANGING_GERMLINE])
    assert categories(known) == dict(
        known_somatic={(8, 9, "C", FOCAL), (13, 14, "T", CO_SOMATIC)},
        known_germline={(7, 8, "C", SILENT_GERMLINE), (12, 13, "A", CHANGING_GERMLINE)},
        unexplained={(0, 1, "G", None)})


def test_without_supplied_variants_every_other_edit_stays_unexplained():
    assert categories(None) == dict(
        known_somatic={(8, 9, "C", FOCAL)}, known_germline=set(),
        unexplained={(0, 1, "G", None), (7, 8, "C", None), (12, 14, "AT", None)})
    # The co-somatic variant alone explains only its own base of the two-base run.
    assert categories(KnownVariants(somatic=[FOCAL, CO_SOMATIC]))["unexplained"] == {
        (0, 1, "G", None), (7, 8, "C", None), (12, 13, "A", None)}


def test_run_split_keeps_leftover_bases_unexplained_and_overlapping_candidates_unused():
    # Only the first base of the GC>AT run is supplied: its second base stays unexplained.
    assert categories(KnownVariants(somatic=[FOCAL], germline=[CHANGING_GERMLINE]))["unexplained"] == {
        (0, 1, "G", None), (7, 8, "C", None), (13, 14, "T", None)}
    # An SNV and an MNV over the same base would each explain it, so neither is used.
    overlapping = KnownVariants(somatic=[FOCAL], germline=[CHANGING_GERMLINE, variant(112, "GC", "AT")])
    assert (12, 14, "AT", variant(112, "GC", "AT")) in categories(overlapping)["known_germline"]
    # In a three-base run, an SNV and an overlapping MNV both fit; neither is used.
    run = translation("CCATGGCT", "AAGATGAA", "MAQDE")
    both = KnownVariants(somatic=[FOCAL], germline=[variant(112, "G", "A"), variant(112, "GC", "AT")])
    assert categories(both, run)["unexplained"] == {(12, 15, "ATG", None)}
    assert repr(both) == "KnownVariants(somatic=1, germline=2)"


def test_protein_exposes_attribution_and_isovar_result_counts_it():
    from isovar.isovar_result import IsovarResult
    from isovar.read_evidence import ReadEvidence

    known = KnownVariants(somatic=[FOCAL, CO_SOMATIC], germline=[SILENT_GERMLINE, CHANGING_GERMLINE])
    protein = ProteinSequence.from_translations([OBSERVED])
    assert protein.known_variants is None and len(protein.unexplained_transcript_edits) == 3
    attributed = protein.with_known_variants(known)
    assert {e.source_variant for e in attributed.known_germline_transcript_edits} == {
        SILENT_GERMLINE, CHANGING_GERMLINE}
    assert attributed.subsequence(1, 3).known_variants is known
    evidence = ReadEvidence(trimmed_base1_start=108, trimmed_ref="A", trimmed_alt="C",
                            ref_reads=[], alt_reads=[], other_reads=[])
    result = IsovarResult(variant=FOCAL, read_evidence=evidence, predicted_effect=None,
                          sorted_protein_sequences=[attributed])
    assert result.co_somatic_variants_in_top_protein_sequence == {CO_SOMATIC}
    assert result.germline_variants_in_top_protein_sequence == {SILENT_GERMLINE, CHANGING_GERMLINE}
    assert (result.num_co_somatic_variants_in_top_protein_sequence,
            result.num_germline_variants_in_top_protein_sequence,
            result.num_unexplained_edits_in_top_protein_sequence) == (1, 2, 1)
    assert {"num_germline_variants_in_top_protein_sequence",
            "num_unexplained_edits_in_top_protein_sequence"} <= set(IsovarResult.record_columns())


def test_repeat_shifted_deletion_matches_the_supplied_representation():
    # Suffix GAAAAC loses one A; the diff and the supplied variant place it differently.
    focal = variant(108, "A", "C")
    t = translation("CCATGGCT", "GAAAC", "MAQ", reference_suffix="GAAAAC", focal=focal,
                    context_after="GAAAAC")
    diff, = categorize_transcript_assembly_edits_from_translation(t, TRANSCRIPT, None)["unexplained"]
    supplied = variant(110, "A", "")
    assert (diff.cdna_start, diff.cdna_end) != (10, 11)
    germline = categories(KnownVariants(somatic=[focal], germline=[supplied]), t)["known_germline"]
    assert germline == {(10, 11, "", supplied)}


def test_equivalent_candidates_of_conflicting_origin_stay_unexplained():
    focal = variant(108, "A", "C")
    t = translation("CCATGGCT", "GAAAC", "MAQ", reference_suffix="GAAAAC", focal=focal,
                    context_after="GAAAAC")
    known = KnownVariants(somatic=[focal, variant(111, "A", "")], germline=[variant(110, "A", "")])
    assert categories(known, t)["known_germline"] == set()
    assert len(categories(known, t)["unexplained"]) == 1


def test_known_variants_rejects_contradictory_inputs():
    with pytest.raises(ValueError, match="both somatic and germline"):
        KnownVariants(somatic=[FOCAL], germline=[FOCAL])
    other_reference = Variant("1", 150, "A", "G", "GRCh37", normalize_contig_names=False)
    with pytest.raises(ValueError, match="different references"):
        KnownVariants(somatic=[FOCAL], germline=[other_reference])
    known = KnownVariants(somatic=[FOCAL], germline=[SILENT_GERMLINE])
    assert known == KnownVariants(somatic=[FOCAL], germline=[SILENT_GERMLINE])
    assert hash(known) == hash(KnownVariants(somatic=[FOCAL], germline=[SILENT_GERMLINE]))
    assert known.overlapping("2", 1, 10**9) == [] and len(known) == 2


def test_export_labels_edit_origins():
    from isovar.isovar_result import IsovarResult
    from isovar.read_evidence import ReadEvidence

    known = KnownVariants(somatic=[FOCAL, CO_SOMATIC], germline=[SILENT_GERMLINE, CHANGING_GERMLINE])
    t = translation("GCATGGCC", "AAGATTAA", "MAQD")
    evidence = ReadEvidence(trimmed_base1_start=108, trimmed_ref="A", trimmed_alt="C",
                            ref_reads=[], alt_reads=[], other_reads=[])
    result = IsovarResult(variant=FOCAL, read_evidence=evidence, predicted_effect=None,
                          sorted_protein_sequences=[ProteinSequence.from_translations([t], known)])
    event, = export_protein_hypotheses([result], sample_id="s", source="rna.bam")["events"]
    assert event["edit_attribution"] == dict(somatic_variants=2, germline_variants=2)
    edits = event["protein_hypotheses"][0]["translations"][0]["observed_edits"]
    assert sorted(e["origin"] for e in edits) == [
        "co_somatic_variant", "germline_variant", "germline_variant", "nominated_variant", "unexplained"]
    assert all((e["source_event_id"] is None) == (e["origin"] == "unexplained") for e in edits)


# Original osteosarc RNA: NTF3's A>G is read as a compound AG>GT, and NR2F2's
# assemblies keep a downstream 3-nt deletion. Which variant set they belong to is
# the caller's statement; Isovar only attributes them.

def case(case_id):
    from examples import osteosarc_assembly_figures as figures
    from tests.data.osteosarc.expansion.references import load_reference

    manifest = json.loads((figures.CORPUS / "manifest.json").read_text())
    selected = next(c for c in manifest["cases"] if c["case_id"] == case_id)
    reference_dir = figures.CORPUS / "references" / selected.get("reference", "GRCh38")
    reference_manifest, _ = load_reference(reference_dir)
    return selected, reference_dir, reference_manifest["variant_transcripts"][selected["variant"]["variant_id"]]


@pytest.fixture(scope="module")
def ntf3(tmp_path_factory):
    from tests.data.osteosarc.expansion.references import reference_genome
    from examples import osteosarc_assembly_figures as figures

    selected, reference_dir, transcripts = case("28-NTF3-chr12-5494381-53f498a544883d51")
    genome = reference_genome(reference_dir, tmp_path_factory.mktemp("ntf3-reference"))
    r = selected["variant"]
    focal = Variant("12", r["pos"], r["ref"], r["alt"], ensembl=genome)
    adjacent = Variant("12", r["pos"] + 1, "G", "T", ensembl=genome)

    def run(**kwargs):
        with pysam.AlignmentFile(str(figures.CORPUS / selected["primary_bam"])) as bam:
            return run_isovar(read_collector=ReadCollector(use_secondary_alignments=False),
                              alignment_file=bam, transcript_id_whitelist=set(transcripts), **kwargs)
    return focal, adjacent, run


def test_ntf3_adjacent_change_is_unexplained_germline_or_co_somatic_as_supplied(ntf3):
    focal, adjacent, run = ntf3
    alone, = run(variants=[focal])
    assert alone.num_unexplained_edits_in_top_protein_sequence == 1
    assert alone.germline_variants_in_top_protein_sequence == set()
    germline, = run(variants=[focal], germline_variants=[adjacent])
    assert germline.germline_variants_in_top_protein_sequence == {adjacent}
    assert germline.num_unexplained_edits_in_top_protein_sequence == 0
    # The protein is the same RNA sequence either way; only the attribution changes.
    assert germline.top_protein_sequence.amino_acids == alone.top_protein_sequence.amino_acids
    results = {r.variant: r for r in run(variants=[focal, adjacent])}
    assert results[focal].co_somatic_variants_in_top_protein_sequence == {adjacent}
    assert results[focal].num_unexplained_edits_in_top_protein_sequence == 0


def test_nr2f2_downstream_deletion_is_attributed_when_supplied(tmp_path):
    from examples import osteosarc_assembly_figures as figures
    from tests.data.osteosarc.expansion.references import reference_genome

    selected, reference_dir, transcripts = case("48-NR2F2-chr15-96332299-fbeb5d8a415e4d18")
    genome = reference_genome(reference_dir, tmp_path)
    r = selected["variant"]
    focal = Variant("15", r["pos"], r["ref"], r["alt"], ensembl=genome)
    # The 3-nt CIGAR deletion the RNA templates carry (0-based [96875576, 96875579)).
    deletion = Variant("15", 96875577, "GTG", "", ensembl=genome)
    with pysam.AlignmentFile(str(figures.CORPUS / selected["primary_bam"])) as bam:
        alone, = run_isovar([focal], bam, read_collector=ReadCollector(use_secondary_alignments=False),
                            transcript_id_whitelist=set(transcripts))
        germline, = run_isovar([focal], bam, read_collector=ReadCollector(use_secondary_alignments=False),
                               transcript_id_whitelist=set(transcripts), germline_variants=[deletion])
    assert alone.num_unexplained_edits_in_top_protein_sequence == 1
    assert germline.germline_variants_in_top_protein_sequence == {deletion}
    assert germline.num_unexplained_edits_in_top_protein_sequence == 0


def test_command_loads_germline_vcf_and_reports_attribution_columns(tmp_path, capsys):
    import pandas as pd
    from isovar.cli.commands import run as isovar_cli
    from .testing_helpers import data_path

    vcf, bam = data_path("data/b16.f10/b16.vcf"), data_path("data/b16.f10/b16.combined.sorted.bam")
    header = [line for line in open(vcf) if line.startswith("#")]
    germline = tmp_path / "germline.vcf"
    germline.write_text("".join(header) + "chr9\t82927200\t.\tA\tG\t.\tPASS\t.\n")
    output = tmp_path / "results.csv"
    isovar_cli(["run", "--vcf", vcf, "--bam", bam, "--germline-vcf", str(germline),
                "--output", str(output), "--log-level", "WARNING"])
    table = pd.read_csv(output)
    assert {"num_germline_variants_in_top_protein_sequence", "num_co_somatic_variants_in_top_protein_sequence",
            "num_unexplained_edits_in_top_protein_sequence"} <= set(table.columns)
    with pytest.raises(SystemExit) as exit_info:
        isovar_cli(["run", "--vcf", vcf, "--bam", bam, "--germline-vcf", str(tmp_path / "missing.vcf"),
                    "--output", str(output)])
    assert exit_info.value.code == 2 and "--germline-vcf" in capsys.readouterr().err


def test_run_isovar_rejects_a_variant_that_is_both_somatic_and_germline():
    from .testing_helpers import data_path

    vcf, bam = data_path("data/b16.f10/b16.vcf"), data_path("data/b16.f10/b16.combined.sorted.bam")
    with pytest.raises(ValueError, match="both somatic and germline"):
        run_isovar(vcf, bam, germline_variants=vcf)
