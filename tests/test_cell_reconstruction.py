"""Discovery can pool cells while attribution retains sparse, ambiguous evidence."""

import json
from types import SimpleNamespace

import pysam
import pytest
from pyensembl import Genome
from varcode import Variant

from isovar import IsovarResult, ProteinSequenceCreator, ReadCollector, reconstruct_cell_groups
from isovar.allele_read import AlleleRead
from isovar.cell_reconstruction import _attribute, _comparison, _protein_assignment
from isovar.dna import reverse_complement_dna
from tests.test_sv_rna import record, write_bam

CDNA = "ATG" + "GCT" * 4 + "AAA" + "GGC" * 4 + "TAA"
MUTANT = CDNA[:15] + "C" + CDNA[16:]
HEADER = pysam.AlignmentHeader.from_dict(dict(
    SQ=[dict(SN="1", LN=20000), dict(SN="2", LN=20000)],
    RG=[dict(ID="a", SM="sample", LB="lib1"), dict(ID="b", SM="sample", LB="lib2"),
        dict(ID="c", SM="other-sample", LB="lib1"), dict(ID="unknown", SM="sample")]))


@pytest.fixture(params=["+", "-"])
def reference(tmp_path, request):
    strand = request.param
    gtf = tmp_path / "reference.gtf"
    attrs = 'gene_id "g"; gene_name "G"; transcript_id "t"; transcript_name "T"; exon_number "1"; exon_id "e"; protein_id "p"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";'
    features = [("gene", 1001, 1033), ("transcript", 1001, 1033), ("exon", 1001, 1033),
                ("CDS", 1001, 1030), ("start_codon", 1001, 1003), ("stop_codon", 1031, 1033)]
    if strand == "-":
        features = [(feature, 2034 - end, 2034 - start) for feature, start, end in features]
    gtf.write_text("".join('1\ttest\t%s\t%d\t%d\t.\t%s\t0\t%s\n' % (feature, start, end, strand, attrs)
                           for feature, start, end in features))
    fasta = tmp_path / "cdna.fa"
    fasta.write_text(">t\n" + CDNA + "\n")
    protein = tmp_path / "protein.fa"
    protein.write_text(">p transcript:t\nMAAAAKGGGG\n")
    genome = Genome(reference_name="cell-test-" + strand, annotation_name="test", annotation_version=1,
                    gtf_path_or_url=str(gtf), transcript_fasta_paths_or_urls=[str(fasta)],
                    protein_fasta_paths_or_urls=[str(protein)], cache_directory_path=str(tmp_path / "cache"))
    genome.index()
    variant = Variant("1", 1016 if strand == "+" else 1018, "A" if strand == "+" else "T",
                      "C" if strand == "+" else "G", ensembl=genome)
    return strand, variant


def reads_for(reference, specifications):
    strand, _ = reference
    records = []
    for i, (sequence, cell, library) in enumerate(specifications):
        sequence = sequence if strand == "+" else reverse_complement_dna(sequence)
        read = record("read%d" % i, "1", 1000, "%dM" % len(sequence), sequence, group=library)
        if cell is not None:
            read.set_tag("CB", cell)
            read.set_tag("UB", "same-umi")
        records.append(read)
    return records


def reconstruct(tmp_path, reference, specifications, **kwargs):
    _, variant = reference
    path = write_bam(tmp_path / "cells.bam", reads_for(reference, specifications), header=HEADER)
    with pysam.AlignmentFile(path) as bam:
        collector = ReadCollector()
        reads = collector.read_evidence_for_variant(variant, bam)
        creator = ProteinSequenceCreator(protein_sequence_length=20, min_transcript_prefix_length=3,
                                         min_variant_sequence_coverage=2, max_protein_sequences_per_variant=1)
        result = IsovarResult(variant, reads, None, protein_sequence_settings=creator.settings())
        return reconstruct_cell_groups([result], bam, sample_id="sample-T1", source="input", **kwargs)


def test_one_read_per_cell_can_support_a_sequence_requiring_pooled_coverage(tmp_path, reference):
    specs = [(MUTANT, "C%d" % i, "a") for i in range(4)]
    pooled = reconstruct(tmp_path, reference, specs)
    event, = pooled["events"]
    pool, = event["pools"]
    protein, = pool["protein_hypotheses"]
    assert protein["amino_acids"] == "MAAAAQGGGG"  # Explicit standard-code expectation.
    assert protein["rna_support"]["cells"] == 4
    candidate, = event["candidates"]
    assert len(candidate["unique_cell_ids"]) == candidate["num_unique_cells"] == 4
    assert event["protein_attribution"][0]["num_unique_cells"] == 4
    assert candidate["sequence"] == (MUTANT if reference[0] == "+" else reverse_complement_dna(MUTANT))
    assert all(c["observations"][0]["rna_support"]["reads"] == 1 for c in event["cells"])
    assert all(c["observations"][0]["status"] == "unique_within_catalog" for c in event["cells"])
    assert event["reconstruction_settings"]["min_variant_sequence_coverage"] == 2
    assert json.loads(json.dumps(pooled)) == pooled
    separate = reconstruct(tmp_path, reference, specs, mode="cell")["events"][0]
    assert len(separate["pools"]) == 4 and not separate["candidates"]
    assert all(p["no_result_reason"] == "no_protein_at_configured_coverage_and_reference_settings"
               for p in separate["pools"])
    assert all(c["observations"][0]["status"] == "no_candidate" for c in separate["cells"])


def test_noisy_singleton_attribution_is_separate_from_discovery(tmp_path, reference):
    noisy = MUTANT[:6] + "A" + MUTANT[7:]
    specs = [(MUTANT, "C1", "a"), (MUTANT, "C2", "a"), (noisy, "noisy", "a"), (CDNA, "ref", "a")]
    event = reconstruct(tmp_path, reference, specs)["events"][0]
    candidate, = event["candidates"]
    assert len(candidate["unique_cell_ids"]) == 3
    cells = {c["cell_barcode"]: c for c in event["cells"]}
    noisy_observation, = cells["noisy"]["observations"]
    assert noisy_observation["status"] == "unique_within_catalog"
    assert noisy_observation["comparisons"][0][0]["edit_distance"] == 1
    assert cells["ref"]["observations"][0]["status"] == "reference_or_other_allele"
    # The noisy cell did not contribute to the reconstructed sequence.
    assert event["pools"][0]["protein_hypotheses"][0]["rna_support"]["reads"] == 2
    strict = reconstruct(tmp_path, reference, specs, max_edits=0)["events"][0]
    assert next(c for c in strict["cells"] if c["cell_barcode"] == "noisy")["observations"][0]["status"] == "outside_tolerance"


def test_library_collision_missing_labels_and_repeated_umis(tmp_path, reference):
    specs = [(MUTANT, "C1", "a"), (MUTANT, "C1", "a"), (MUTANT, "C1", "b"),
             (MUTANT, "C1", "b"), (MUTANT, None, "a"), (MUTANT, "C1", "unknown")]
    event = reconstruct(tmp_path, reference, specs, mode="library")["events"][0]
    assert len(event["cells"]) == 3
    assert len({c["cell_id"] for c in event["cells"]}) == 3
    assert event["unresolved_cell_support"]["reads"] == 1
    assert event["input_support"]["reads"] == 6
    assert sum(p["rna_support"]["reads"] for p in event["pools"]) == 6
    lib1 = next(p for p in event["pools"] if p["description"]["scope"]["library"] == "lib1")
    assert lib1["rna_support"]["reads"] == 3 and lib1["rna_support"]["umis"] == 1
    assert not next(c for c in event["cells"] if c["scope"]["library"] is None)["library_scope_known"]


def test_explicit_groups_and_samples_do_not_cross_pool_boundaries(tmp_path, reference):
    specs = [(MUTANT, "C1", "a"), (MUTANT, "C2", "a"), (MUTANT, "C3", "a"), (MUTANT, "C1", "c")]
    event = reconstruct(tmp_path, reference, specs)["events"][0]
    assert len(event["pools"]) == 2
    foreign = next(c for c in event["cells"] if c["scope"]["header_sample"] == "other-sample")
    assert foreign["observations"][0]["status"] == "no_candidate"
    groups = {c["cell_id"]: "chosen" for c in event["cells"] if c["cell_barcode"] in ("C1", "C2")}
    grouped = reconstruct(tmp_path, reference, specs, mode="group", cell_groups=groups)["events"][0]
    assert grouped["excluded_from_pool"]["reads"] == 1
    assert len(grouped["pools"]) == 2
    assert sum(p["rna_support"]["reads"] for p in grouped["pools"]) == 3
    # Exclusion from discovery doesn't erase the cell's original evidence.
    assert len(grouped["cells"]) == 4


TX = SimpleNamespace(id="tx", exon_intervals=[(1, 20000)])


def candidate(key, prefix="AAACCC", alt="G", suffix="TTTAAA"):
    return dict(candidate_id=key, sequence=prefix + alt + suffix, variant_interval=[len(prefix), len(prefix) + len(alt)])


def observation(prefix="AAACCC", allele="G", suffix="TTTAAA"):
    return AlleleRead(prefix, allele, suffix, "read")


def assign(reads, candidates, **kwargs):
    return _attribute(reads, candidates, {c["candidate_id"]: [TX] for c in candidates}, "A",
                      kwargs.get("max_edits", 1), kwargs.get("min_overlap", 3))


def test_partial_read_does_not_choose_between_distant_differences():
    candidates = [candidate("a"), candidate("b", suffix="TTTCCC")]
    result = assign([observation(suffix="TTT")], candidates)
    assert result["status"] == "ambiguous" and result["candidate_ids"] == ["a", "b"]
    assert all(not c["spans_candidate"] for c in result["comparisons"][0])
    assert assign([observation()], candidates)["candidate_ids"] == ["a"]


def test_shorter_candidate_does_not_lose_due_to_unobserved_flanks():
    result = assign([observation()], [candidate("full"), candidate("short", prefix="CCC", suffix="TTT")])
    assert result["status"] == "ambiguous"


def test_error_at_focal_base_cannot_be_rescued_by_pooled_support():
    for allele in ("A", "T"):
        result = assign([observation(allele=allele)], [candidate("a")])
        assert result["status"] == "reference_or_other_allele" and not result["candidate_ids"]


def test_alternative_placements_do_not_become_independent_or_unique_support():
    result = assign([observation(), observation(allele="A")], [candidate("a")])
    assert result["status"] == "ambiguous"


def test_insufficient_overlap_and_non_acgt_are_unassessed():
    for read in (observation(prefix="", suffix=""), observation(prefix="NNNNNN")):
        assert assign([read], [candidate("a")])["status"] == "insufficient_coverage_or_incompatible_path"


def test_indel_competitor_and_observed_splice_path():
    insertion = candidate("insertion", alt="GG")
    assert _comparison(observation(allele="GG"), insertion, "", [TX], 1, 3)["reference_edit_distance"] == 2
    incompatible = AlleleRead("AAACCC", "G", "TTTAAA", "splice", reference_blocks=((0, 6, 10, 16), (6, 13, 100, 107)),
                              splice_junctions=((16, 100),))
    assert assign([incompatible], [candidate("a")])["status"] == "insufficient_coverage_or_incompatible_path"


@pytest.mark.parametrize("options", [dict(mode="bad"), dict(max_edits=-1), dict(max_edits=True),
                                    dict(min_overlap=0), dict(mode="group"), dict(cell_groups={}),
                                    dict(mode="group", cell_groups={"c": ""})])
def test_invalid_options_fail_before_reading_alignments(options):
    with pytest.raises(ValueError):
        reconstruct_cell_groups([], None, sample_id="s", source="b", **options)


def test_alternative_sequences_survive_pooling_and_a_cell_can_support_both(tmp_path, reference):
    alternative = MUTANT[:21] + "TTC" + MUTANT[24:]
    specs = [(MUTANT, "both", "a"), (MUTANT, "first", "a"),
             (alternative, "both", "a"), (alternative, "second", "a")]
    event = reconstruct(tmp_path, reference, specs)["events"][0]
    assert {p["amino_acids"] for p in event["pools"][0]["protein_hypotheses"]} == {"MAAAAQGGGG", "MAAAAQGFGG"}
    assert len(event["candidates"]) == 2
    both = next(c for c in event["cells"] if c["cell_barcode"] == "both")
    assert all(both["cell_id"] in c["unique_cell_ids"] for c in event["candidates"])
    assert len(both["observations"]) == 2


def test_empty_and_unlabelled_inputs_are_accounted_for(tmp_path, reference):
    empty = reconstruct(tmp_path, reference, [])["events"][0]
    assert not empty["pools"] and not empty["cells"] and empty["input_support"]["reads"] == 0
    unlabelled = reconstruct(tmp_path, reference, [(MUTANT, None, "a")] * 2, mode="cell")["events"][0]
    assert not unlabelled["pools"] and not unlabelled["cells"]
    assert unlabelled["excluded_from_pool"]["reads"] == unlabelled["unresolved_cell_support"]["reads"] == 2


def test_command_export_is_opt_in_and_contains_scoped_attribution(tmp_path):
    from isovar.cli.commands import run
    from tests.test_cell_evidence import ALT, REF, HEADER as COUNTS_HEADER

    bam = write_bam(tmp_path / "rna.bam", ALT + REF, header=COUNTS_HEADER)
    output = tmp_path / "proteins.json"
    args = ["protein-hypotheses", "--variant", "1", "1021", "A", "G", "--genome", "GRCh38",
            "--bam", str(bam), "--sample-id", "sample", "--output", str(output), "--log-level", "ERROR"]
    run(args)
    assert "cell_reconstruction" not in json.loads(output.read_text())
    run(args + ["--cell-reconstruction", "library"])
    export = json.loads(output.read_text())["cell_reconstruction"]
    assert export["schema"] == "isovar.cell_reconstruction.v1" and export["mode"] == "library"
    assert export["events"][0]["input_support"]["reads"] == 6
    assert len(export["events"][0]["cells"]) == 3
    for extra in (["--cell-max-edits", "1"], ["--cell-reconstruction", "group"],
                  ["--cell-reconstruction", "cell", "--cell-min-overlap", "0"]):
        with pytest.raises(SystemExit) as error:
            run(args + extra)
        assert error.value.code == 2


def test_original_ont_single_cell_records_and_pinned_reference(tmp_path, monkeypatch):
    import socket
    from examples import osteosarc_assembly_figures as figures
    from isovar import run_isovar
    from tests.data.osteosarc.expansion.references import load_reference, reference_genome, apply_variant

    manifest = json.loads((figures.CORPUS / "manifest.json").read_text())
    case = next(c for c in manifest["cases"] if c["case_id"] == "26-NME1-chr17-51154443-53f498a544883d51")
    bam_path = figures.corpus_file(case["primary_bam"])
    def offline(*args, **kwargs):
        raise AssertionError("Cell reconstruction must run offline")
    monkeypatch.setattr(socket, "socket", offline)
    reference_dir = figures.CORPUS / "references" / "GRCh38"
    reference_manifest, models = load_reference(reference_dir)
    genome = reference_genome(reference_dir, tmp_path)
    record = case["variant"]
    variant = Variant("17", record["pos"], record["ref"], record["alt"], ensembl=genome)
    whitelist = set(reference_manifest["variant_transcripts"][record["variant_id"]])
    collector = ReadCollector(use_secondary_alignments=False)
    with pysam.AlignmentFile(str(bam_path)) as bam:
        results = run_isovar([variant], bam, read_collector=collector, transcript_id_whitelist=whitelist)
        export = reconstruct_cell_groups(results, bam, sample_id="sid-T1", source="nme1-original",
                                         read_collector=collector, transcript_id_whitelist=whitelist)
    event, = export["events"]
    assert event["input_support"]["reads"] == 68
    assert event["unresolved_cell_support"]["reads"] == 0
    assert sum(o["rna_support"]["reads"] for c in event["cells"] for o in c["observations"]) == 68
    assert all(not c["library_scope_known"] for c in event["cells"])
    proteins = [p for pool in event["pools"] for p in pool["protein_hypotheses"]]
    assert proteins
    for protein in proteins:
        # Independent transcript editing and translation from the pinned source
        # reference, not the output of the reconstruction algorithm.
        expected = [apply_variant(record, models[t])["mutant_protein"] for t in protein["transcript_ids"]]
        assert all(protein["amino_acids"] in sequence for sequence in expected)
    assert any(c["unique_cell_ids"] or c["ambiguous_cell_ids"] for c in event["candidates"])


def test_synonymous_rna_ambiguity_can_identify_one_protein():
    candidates = [dict(candidate_id="a", hypothesis_ids=["protein"]),
                  dict(candidate_id="b", hypothesis_ids=["protein"])]
    attribution = dict(status="ambiguous", candidate_ids=["a", "b"], placement_ambiguous=False)
    assert _protein_assignment(attribution, candidates) == dict(hypothesis_ids=["protein"],
                                                              protein_status="unique_within_catalog")
    attribution["placement_ambiguous"] = True
    assert _protein_assignment(attribution, candidates)["protein_status"] == "ambiguous"
    attribution["placement_ambiguous"] = False
    candidates[1]["hypothesis_ids"].append("other-frame")
    assert _protein_assignment(attribution, candidates)["protein_status"] == "ambiguous"


def test_quality_assessments_preserve_pooled_orfs_and_sparse_cell_attribution(tmp_path, reference):
    from isovar import BaseQualityPolicy, export_protein_hypotheses

    noisy = MUTANT[:6] + "A" + MUTANT[7:]
    specs = [(MUTANT, "high", "a"), (MUTANT, "boundary", "a"),
             (noisy, "low-noisy", "a"), (MUTANT, "unknown", "a")]
    records = reads_for(reference, specs)
    for read, quality in zip(records, (40, 20, 0, None)):
        read.query_qualities = None if quality is None else [quality] * len(read.query_sequence)
    path = write_bam(tmp_path / "quality.bam", records, header=HEADER)
    _, variant = reference
    with pysam.AlignmentFile(path) as bam:
        reads = ReadCollector().read_evidence_for_variant(variant, bam)
        creator = ProteinSequenceCreator(protein_sequence_length=20, min_transcript_prefix_length=3,
                                         min_variant_sequence_coverage=2)
        proteins = creator.sorted_protein_sequences_for_variant(variant, reads)
        result = IsovarResult(variant, reads, proteins, protein_sequence_settings=creator.settings())
        raw = reconstruct_cell_groups([result], bam, sample_id="s", source="input")
        assessed = reconstruct_cell_groups([result], bam, sample_id="s", source="input",
                                           base_quality_policy=BaseQualityPolicy(20))
        sample_raw = export_protein_hypotheses([result], sample_id="s", source="input")
        sample_assessed = export_protein_hypotheses([result], sample_id="s", source="input",
                                                   base_quality_policy=BaseQualityPolicy(20))
    event, = assessed["events"]
    candidate, = event["candidates"]
    assert candidate["num_unique_cells"] == 4
    assert [candidate["base_quality"][s]["num_unique_cells"] for s in ("passed", "failed", "unassessed")] == [2, 1, 1]
    assert event["protein_attribution"][0]["base_quality"] == candidate["base_quality"]
    protein, = event["pools"][0]["protein_hypotheses"]
    assert protein["amino_acids"] == "MAAAAQGGGG"
    assert protein["rna_support"]["reads"] == 3  # noisy singleton only attributed afterward
    assert protein["rna_support"]["base_quality"]["unassessed"]["reads"] == 1
    assert protein["rna_support"]["base_quality"]["failed"]["reads"] == 0
    for cell in event["cells"]:
        observation, = cell["observations"]
        assert observation["status"] == "unique_within_catalog"
        support = observation["rna_support"]
        assert support["reads"] == 1
        for status in ("passed", "failed", "unassessed"):
            subset = support["base_quality"][status]
            evidence = assessed["evidence_sets"][subset["evidence_set_id"]]
            assert len(evidence["read_ids"]) == subset["reads"]
    assert json.loads(json.dumps(assessed)) == assessed

    def without_quality(value):
        if isinstance(value, dict):
            return {key: without_quality(item) for key, item in value.items()
                    if key not in ("base_quality", "base_quality_policy", "evidence_sets")}
        return [without_quality(v) for v in value] if isinstance(value, list) else value

    assert without_quality(raw) == without_quality(assessed)
    assert without_quality(sample_raw) == without_quality(sample_assessed)
    # No QUAL still gives the same ORFs and sparse assignments through one pipeline.
    for read in records:
        read.query_qualities = None
    missing_path = write_bam(tmp_path / "missing-quality.bam", records, header=HEADER)
    with pysam.AlignmentFile(missing_path) as bam:
        missing_reads = ReadCollector().read_evidence_for_variant(variant, bam)
        missing_result = IsovarResult(variant, missing_reads, None, protein_sequence_settings=creator.settings())
        missing = reconstruct_cell_groups([missing_result], bam, sample_id="s", source="input",
                                          base_quality_policy=BaseQualityPolicy(20))
    assert without_quality(missing) == without_quality(raw)
    assert missing["events"][0]["candidates"][0]["base_quality"]["unassessed"]["num_unique_cells"] == 4
