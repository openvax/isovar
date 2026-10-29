"""Known transcript mixtures, original fragment identity and coding sequences."""

from copy import deepcopy
from dataclasses import replace
from hashlib import md5
import json
from pathlib import Path
import random

import pysam
import pytest
from pyensembl import Genome

from isovar import IsoformTarget, reconstruct_antigen_isoforms
from isovar.antigen_reference import (
    gdc_genome, gdc_reference_dictionary, isoform_targets_from_genome, validate_alignment_reference,
)
from isovar.cli.commands import run
from isovar.dna import reverse_complement_dna
from isovar.fusion import FusionReference
from isovar.isoform_quantification import add_isoform_abundances, _augmented_transcriptome, write_isoform_fastas
from tests.test_sv_rna import HEADER, aligned, record, write_bam


def target_models(strand="+"):
    rng = random.Random(465)
    codons = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in {"TAA", "TAG", "TGA"}]
    first = "A" * 30 + "ATG" + "".join(rng.choice(codons) for _ in range(49))
    middle = "".join(rng.choice(codons) for _ in range(30))
    last = "".join(rng.choice(codons) for _ in range(49)) + "TAA" + "C" * 30
    exons = ((100, 280), (400, 490), (600, 780))
    sequences = [first + middle + last, first + last]
    paths = [exons, (exons[0], exons[2])]
    targets, positions = [], []
    for i, (seq, blocks) in enumerate(zip(sequences, paths)):
        if strand == "-":
            blocks = tuple(sorted((1000 - b, 1000 - a) for a, b in blocks))
        reference = FusionReference("T%d" % (i + 1), "test", "test:1", "1", strand, blocks, seq, 30, len(seq) - 30)
        targets.append(IsoformTarget("G", "CTA", reference))
        mapping = [("1", p, strand) for a, b in blocks for p in range(a, b)]
        positions.append(mapping[::-1] if strand == "-" else mapping)
    return targets, positions


def reconstruct(tmp_path, reads, targets=None, **kwargs):
    targets = target_models()[0] if targets is None else targets
    path = write_bam(tmp_path / "reads.bam", reads)
    with pysam.AlignmentFile(path) as bam:
        return reconstruct_antigen_isoforms(bam, targets=targets, sample_id="sample", source="library",
                                           scratch_dir=tmp_path, **kwargs)


def estimates(tmp_path, result, weights):
    path = tmp_path / "quant.sf"
    with path.open("w") as out:
        out.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
        for row, weight in zip(result["isoforms"], weights):
            out.write("%s\t%d\t%d\t%g\t%g\n" % (row["id"], len(row["sequence"]), len(row["sequence"]) - 50,
                                                      weight * 1000000, weight * 100))
    return path


@pytest.mark.parametrize("strand", ["+", "-"])
def test_two_isoform_mixture_reconstructs_exact_sequences_and_counts_once(tmp_path, strand):
    targets, paths = target_models(strand)
    reads = [read for i, n in enumerate((8, 2)) for j in range(n)
             for read in aligned("isoform%d-%d" % (i, j), targets[i].reference.sequence, paths[i])]
    reads += [read for j in range(3) for read in aligned("shared%d" % j, targets[0].reference.sequence[:100], paths[0][:100])]
    result = reconstruct(tmp_path, reads, targets)
    assert {r["sequence"] for r in result["isoforms"]} == {t.reference.sequence for t in targets}
    assert not result["partial_sequences"]
    by_id = {r["annotated_transcript_ids"][0]: r for r in result["isoforms"]}
    assert [by_id[t]["unique_fragments"] for t in ("T1", "T2")] == [8, 2]
    assert [by_id[t]["compatible_fragments"] for t in ("T1", "T2")] == [11, 5]
    assert result["support"]["target_fragments"] == 13
    assert sum(c["fragments"] for c in result["support"]["equivalence_classes"]) == 13
    for row in by_id.values():
        assert row["sequence_status"] == "whole_path_observed"
        assert row["coding"][0]["status"] == "complete_CDS"
        # Independent BioPython translation, including the observed start/stop.
        from Bio.Seq import Seq
        assert row["coding"][0]["amino_acids"] == str(Seq(row["sequence"][30:-30]).translate(to_stop=True))
    quant = estimates(tmp_path, result, [0.8 if r["id"] == by_id["T1"]["id"] else 0.2 for r in result["isoforms"]])
    add_isoform_abundances(result, quant, provenance={"fixture": "known-mixture"})
    gene, = result["quantification"]["genes"]
    assert gene["dominant_isoform_id"] == by_id["T1"]["id"]
    assert gene["isoform_fractions"][by_id["T1"]["id"]] == pytest.approx(0.8)
    # A missing annotated competitor prevents a resolved call, even when the
    # available transcripts have distinguishing reads and unequal estimates.
    result["target_selection"] = dict(genes=["CTA"], unavailable_transcripts=[
        dict(gene_id="G", transcript_id="T3", reason="ambiguous reference bases")])
    add_isoform_abundances(result, quant, provenance={"fixture": "missing-competitor"})
    assert result["quantification"]["genes"][0]["dominant_isoform_id"] is None


def test_partial_reads_do_not_fill_reference_or_resolve_ambiguous_dominance(tmp_path):
    targets, paths = target_models()
    observed = targets[0].reference.sequence[:120]
    result = reconstruct(tmp_path, aligned("partial", observed, paths[0][:120]))
    partial, = result["partial_sequences"]
    assert partial["sequence"] == observed
    assert all(r["origin"] == "reference_model" and r["sequence_status"] == "partially_covered_reference"
               for r in result["isoforms"])
    quant = estimates(tmp_path, result, [0.99, 0.01])
    add_isoform_abundances(result, quant, provenance={"fixture": "ambiguous"})
    gene, = result["quantification"]["genes"]
    assert gene["dominant_isoform_id"] is None
    assert len(gene["groups"]) == 1
    assert gene["groups"][0]["fraction"] == 1


def test_observed_sequence_allele_keeps_exact_bases_and_original_quality_metadata(tmp_path):
    targets, paths = target_models()
    original = targets[0].reference.sequence
    mutant = original[:75] + ("A" if original[75] != "A" else "T") + original[76:]
    reads = aligned("reference", original, paths[0]) + aligned("allele", mutant, paths[0])
    reads[-1].query_qualities = None
    reads[-1].set_tag("CB", "cell")
    reads[-1].set_tag("UB", "umi")
    result = reconstruct(tmp_path, reads)
    novel, = [r for r in result["isoforms"] if r["origin"] == "observed_assembly"]
    assert novel["sequence"] == mutant
    assert novel["whole_path_segments"] == [["a", "allele", 0]]
    assert novel["unique_fragments"] == 1
    witness, = [w for w in result["support"]["witnesses"] if w["fragment"][1] == "allele"]
    assert witness["records"][0]["base_qualities_available"] is False
    assert witness["records"][0]["tags"]["CB"] == "cell"


def test_discovery_limit_never_caps_support_and_prevents_dominance_claim(tmp_path):
    targets, paths = target_models()
    reads = [r for i in range(20) for r in aligned("r%02d" % i, targets[0].reference.sequence, paths[0])]
    result = reconstruct(tmp_path, reads, max_discovery_reads=1)
    assert result["discovery"]["truncated"]
    row, = [r for r in result["isoforms"] if "T1" in r["annotated_transcript_ids"]]
    assert row["unique_fragments"] == 20
    assert result["support"]["target_fragments"] == 20
    add_isoform_abundances(result, estimates(tmp_path, result, [0.5, 0.5]), provenance={"fixture": "limited"})
    assert result["quantification"]["genes"][0]["dominant_isoform_id"] is None


def test_paralog_secondary_alignments_and_read_groups_preserve_ambiguity(tmp_path):
    targets, paths = target_models()
    target = targets[0]
    paralog = IsoformTarget("P", "PARALOG", replace(target.reference, transcript_id="P1", contig="2"))
    other_path = [("2", p, s) for _, p, s in paths[0]]
    reads = aligned("same", target.reference.sequence, paths[0])
    reads += aligned("same", target.reference.sequence, other_path, flag=256)
    reads += aligned("same", target.reference.sequence, paths[0], group="b")
    result = reconstruct(tmp_path, reads, [target, paralog])
    assert result["support"]["target_fragments"] == 2
    first = next(w for w in result["support"]["witnesses"] if w["fragment"][0] == "a")
    assert len(first["compatible_isoforms"]) == 2
    assert len(first["records"]) == 2
    add_isoform_abundances(result, estimates(tmp_path, result, [0.8, 0.2]), provenance={"fixture": "paralog"})
    assert all(g["dominant_isoform_id"] is None for g in result["quantification"]["genes"])


def test_mates_and_repeated_payloads_do_not_inflate_fragments(tmp_path):
    targets, paths = target_models()
    seq = targets[0].reference.sequence
    reads = aligned("pair", seq[:210], paths[0][:210], flag=65)
    reads += aligned("pair", seq[-210:], paths[0][-210:], flag=129)
    result = reconstruct(tmp_path, reads + reads)
    assert result["input"]["records"] == 4 and result["input"]["unique_records"] == 2
    assert result["support"]["target_fragments"] == 1
    row = next(r for r in result["isoforms"] if "T1" in r["annotated_transcript_ids"])
    assert row["unique_fragments"] == 1


def test_unreported_competing_alignment_prevents_unique_attribution(tmp_path):
    targets, paths = target_models()
    reads = aligned("multi", targets[0].reference.sequence, paths[0])
    reads[0].set_tag("NH", 2)
    result = reconstruct(tmp_path, reads)
    assert all(r["unique_fragments"] == 0 for r in result["isoforms"])
    assert result["support"]["witnesses"][0]["competing_placement_outside_targets"]
    add_isoform_abundances(result, estimates(tmp_path, result, [0.9, 0.1]), provenance={"fixture": "NH"})
    assert result["quantification"]["genes"][0]["dominant_isoform_id"] is None


def test_identical_proteins_sum_transcript_estimates_once(tmp_path):
    targets, paths = target_models()
    original = targets[0]
    r = original.reference
    shorter = IsoformTarget("G", "CTA", replace(r, transcript_id="T3", sequence=r.sequence[15:],
        exons=((115, 280), *r.exons[1:]), cds_start=15, cds_end=r.cds_end - 15))
    reads = aligned("full", r.sequence, paths[0]) + aligned("short", r.sequence[15:], paths[0][15:])
    result = reconstruct(tmp_path, reads, [original, shorter])
    add_isoform_abundances(result, estimates(tmp_path, result, [0.8, 0.2]), provenance={"fixture": "UTR"})
    protein, = result["quantification"]["proteins"]
    assert protein["tpm"] == 1000000
    assert protein["estimated_fragments"] == 100
    assert len(protein["isoform_ids"]) == 2
    observed, = result["protein_groups"]
    assert observed["compatible_fragments"] == 2


def test_complete_cds_does_not_require_observed_utrs(tmp_path):
    targets, paths = target_models()
    r = targets[0].reference
    result = reconstruct(tmp_path, aligned("cds", r.sequence[30:-30], paths[0][30:-30]), [targets[0]])
    model, = result["isoforms"]
    assert model["whole_path_segments"] == []
    coding, = model["coding"]
    assert coding["all_bases_observed"] and coding["whole_cds_segments"] == [["a", "cds", 0]]
    write_isoform_fastas(result, tmp_path / "rna.fa", tmp_path / "protein.fa")
    assert coding["amino_acids"] in (tmp_path / "protein.fa").read_text()


def test_uncovered_reference_cds_is_not_exported_as_observed_protein(tmp_path):
    targets, paths = target_models()
    r = targets[0].reference
    result = reconstruct(tmp_path, aligned("utr", r.sequence[:20], paths[0][:20]), [targets[0]])
    write_isoform_fastas(result, tmp_path / "rna.fa", tmp_path / "protein.fa")
    assert (tmp_path / "protein.fa").read_text() == ""
    add_isoform_abundances(result, estimates(tmp_path, result, [1.0]), provenance={"fixture": "UTR-only"})
    protein, = result["quantification"]["proteins"]
    assert protein["observed_coding_tpm"] == 0
    assert protein["uncovered_reference_coding_tpm"] == 1000000


def test_start_loss_and_missing_cds_never_acquire_an_invented_start(tmp_path):
    targets, paths = target_models()
    r = targets[0].reference
    lost = r.sequence[:30] + "CCC" + r.sequence[33:]
    result = reconstruct(tmp_path, aligned("lost", lost, paths[0]), [targets[0]])
    novel, = [r for r in result["isoforms"] if r["origin"] == "observed_assembly"]
    assert all(c["status"] == "partial_CDS" for c in novel["coding"])
    unknown = IsoformTarget("G", "CTA", replace(r, cds_start=None, cds_end=None))
    assert unknown.coding_status == "no_complete_annotated_CDS"
    result = reconstruct(tmp_path, aligned("no-cds", r.sequence, paths[0]), [unknown])
    assert result["protein_groups"] == [] and all(not row["coding"] for row in result["isoforms"])


def test_novel_splice_path_preserves_observed_coding_sequence(tmp_path):
    targets, paths = target_models()
    # A newly observed internal splice removes 30 coding bases without using
    # any reference fill or introducing a nominated variant.
    r = targets[0].reference
    seq = r.sequence[:210] + r.sequence[240:]
    positions = paths[0][:210] + paths[0][240:]
    result = reconstruct(tmp_path, aligned("novel-splice", seq, positions), [targets[0]])
    novel, = [row for row in result["isoforms"] if row["origin"] == "observed_assembly"]
    assert novel["sequence"] == seq and novel["reaches_annotated_transcript_ends"]
    from Bio.Seq import Seq
    assert novel["coding"][0]["amino_acids"] == str(Seq(seq[30:-30]).translate(to_stop=True))


@pytest.mark.parametrize("bad", ["nan", "-1", "inf"])
def test_invalid_abundance_does_not_mutate_result(tmp_path, bad):
    targets, paths = target_models()
    result = reconstruct(tmp_path, aligned("r", targets[0].reference.sequence, paths[0]))
    quant = estimates(tmp_path, result, [0.8, 0.2])
    quant.write_text(quant.read_text().replace("800000", bad))
    with pytest.raises(ValueError, match="Invalid Salmon estimate"):
        add_isoform_abundances(result, quant, provenance={"fixture": "bad"})
    assert result["quantification"] is None


def test_background_export_replaces_exact_target_models_and_keeps_competitors(tmp_path):
    targets, paths = target_models()
    result = reconstruct(tmp_path, aligned("r", targets[0].reference.sequence, paths[0]))
    source = tmp_path / "all.fa"
    source.write_text("".join(">%s|gencode\n%s\n" % (t.reference.transcript_id, t.reference.sequence) for t in targets)
                      + ">unrelated\nACGT\n")
    output = tmp_path / "augmented.fa"
    _augmented_transcriptome(result, source, output)
    with pysam.FastxFile(output) as fasta:
        records = {r.name: r.sequence for r in fasta}
    assert records.pop("unrelated") == "ACGT"
    assert records == {r["id"]: r["sequence"] for r in result["isoforms"]}
    source.write_text(source.read_text().replace(targets[0].reference.sequence, "ACGT"))
    with pytest.raises(ValueError, match="differs"):
        _augmented_transcriptome(result, source, output)


def reference_fixture(tmp_path):
    targets, paths = target_models()
    bases = list("N" * 20000)
    for t, mapping in zip(targets, paths):
        for base, (_, p, _) in zip(t.reference.sequence, mapping):
            bases[p] = base
    sequence = "".join(bases)
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + sequence + "\n>chr2\n" + "N" * 20000 + "\n")
    gtf = tmp_path / "reference.gtf"
    with gtf.open("w") as out:
        def feature(kind, start, end, attrs):
            out.write('chr1\ttest\t%s\t%d\t%d\t.\t+\t.\t%s\n' % (kind, start, end, attrs))
        base = 'gene_id "G"; gene_name "CTA"; gene_biotype "protein_coding";'
        feature("gene", 101, 780, base)
        for target in targets:
            attrs = base + ' transcript_id "%s"; transcript_biotype "protein_coding";' % target.reference.transcript_id
            feature("transcript", 101, 780, attrs)
            for i, (a, b) in enumerate(target.reference.exons):
                feature("exon", a + 1, b, attrs + ' exon_number "%d"; exon_id "%s%d";' % (i + 1, target.reference.transcript_id, i))
                lo, hi = max(a, 130), min(b, 747)
                feature("CDS", lo + 1, hi, attrs)
            feature("start_codon", 131, 133, attrs)
            feature("stop_codon", 748, 750, attrs)
    header = HEADER.to_dict()
    header["SQ"] = [dict(SN="chr1", LN=20000, M5=md5(sequence.encode()).hexdigest()),
                    dict(SN="chr2", LN=20000)]
    bam = write_bam(tmp_path / "input.bam", aligned("read", targets[0].reference.sequence, paths[0]),
                    pysam.AlignmentHeader.from_dict(header))
    return fasta, gtf, bam


def test_pyensembl_adapter_and_cli_use_same_genome_dna_and_exon_coordinates(tmp_path):
    fasta, gtf, bam_path = reference_fixture(tmp_path)
    with gdc_genome(genome_fasta=fasta, gtf=gtf, annotation_version="test", cache_directory=tmp_path / "cache") as genome:
        genome.index()
        targets, missing = isoform_targets_from_genome(genome, ["CTA"], contig_names=["chr1", "chr2"])
        assert not missing and len(targets) == 2
        assert targets[0].reference.sequence == target_models()[0][0].reference.sequence
        assert targets[0].reference.cds_start == 30
        with pysam.AlignmentFile(bam_path) as bam:
            evidence = validate_alignment_reference(bam, genome, contigs=["1"])
            assert evidence["checked_sequence_contigs"] == ["chr1"]
            with pytest.raises(ValueError, match="dictionary"):
                validate_alignment_reference(bam, genome, profile="GRCh38.d1.vd1")
    output = tmp_path / "output"
    run(["isoforms", "reconstruct", "--bam", str(bam_path), "--gene", "CTA", "--sample-id", "sample",
         "--genome-fasta", str(fasta), "--gtf", str(gtf), "--annotation-version", "test",
         "--alignment-reference", "custom", "--reference-cache", str(tmp_path / "cli-cache"), "--output-dir", str(output)])
    result = json.loads((output / "isoforms.json").read_text())
    assert len(result["isoforms"]) == 2
    assert result["reference"]["files"]["gtf"]["sha256"]
    assert (output / "proteins.fa").is_file()


@pytest.mark.parametrize("bam_name,fasta_name", [
    ("1", "chr1"), ("KI270728.1", "chr16_KI270728v1_random"),
])
def test_aliases_reach_transcript_sequence_extraction_and_cli(tmp_path, bam_name, fasta_name):
    fasta, gtf, original = reference_fixture(tmp_path)
    fasta.write_text(fasta.read_text().replace(">chr1\n", ">%s\n" % fasta_name))
    gtf.write_text(gtf.read_text().replace("chr1\t", bam_name + "\t"))
    bam_path = tmp_path / "renamed.bam"
    with pysam.AlignmentFile(original) as source:
        header = source.header.to_dict()
        header["SQ"][0]["SN"] = bam_name
        with pysam.AlignmentFile(bam_path, "wb", header=header) as destination:
            for read in source:
                destination.write(read)
    output = tmp_path / "aliased-output"
    run(["isoforms", "reconstruct", "--bam", str(bam_path), "--gene", "CTA", "--sample-id", "sample",
         "--genome-fasta", str(fasta), "--gtf", str(gtf), "--annotation-version", "test",
         "--alignment-reference", "custom", "--reference-cache", str(tmp_path / "alias-cache"),
         "--output-dir", str(output)])
    result = json.loads((output / "isoforms.json").read_text())
    assert result["support"]["target_fragments"] == 1
    assert {target["reference"]["contig"] for target in result["targets"]} == {bam_name}
    assert {row["sequence"] for row in result["isoforms"]} == {
        target.reference.sequence for target in target_models()[0]}
    assert result["reference"]["checked_sequence_contigs"] == [bam_name]


def test_gdc_factory_remains_a_normal_pyensembl_transcript_provider(tmp_path):
    fasta, gtf, _ = reference_fixture(tmp_path)
    cdna = tmp_path / "transcripts.fa"
    cdna.write_text("".join(">%s|GENCODE\n%s\n" % (t.reference.transcript_id, t.reference.sequence)
                             for t in target_models()[0]))
    with gdc_genome(genome_fasta=fasta, gtf=gtf, annotation_version="test", transcript_fasta=cdna,
                    cache_directory=tmp_path / "cache") as genome:
        genome.index()
        assert isinstance(genome, Genome)
        transcript = genome.transcript_by_id("T1")
        assert transcript.sequence == target_models()[0][0].reference.sequence
        assert transcript.coding_sequence == transcript.sequence[30:-30]


def test_ensembl_separate_transcript_versions_match_versioned_fasta_ids(tmp_path):
    fasta, gtf, _ = reference_fixture(tmp_path)
    gtf.write_text(gtf.read_text().replace('transcript_biotype "protein_coding";',
                                         'transcript_biotype "protein_coding"; transcript_version "2";'))
    with gdc_genome(genome_fasta=fasta, gtf=gtf, annotation_version="test", cache_directory=tmp_path / "cache") as genome:
        genome.index()
        targets, unavailable = isoform_targets_from_genome(genome, ["CTA"], contig_names=["chr1", "chr2"])
    assert not unavailable
    assert {t.reference.transcript_id for t in targets} == {"T1.2", "T2.2"}


def test_identical_annotations_share_cache_but_changed_contents_never_reuse_old_models(tmp_path):
    fasta, gtf, _ = reference_fixture(tmp_path)
    kwargs = dict(genome_fasta=fasta, gtf=gtf, annotation_version="test", cache_directory=tmp_path / "cache")
    with gdc_genome(**kwargs) as first:
        first.index()
        first_cache = first.cache_directory_path
        assert first.gene_by_id("G").name == "CTA"
    with gdc_genome(**kwargs) as same:
        assert same.cache_directory_path == first_cache
    gtf.write_text(gtf.read_text().replace('gene_name "CTA"', 'gene_name "CHANGED"'))
    with gdc_genome(**kwargs) as changed:
        changed.index()
        assert changed.cache_directory_path != first_cache
        assert changed.gene_by_id("G").name == "CHANGED"


def test_reference_checksum_mismatch_is_rejected_even_when_lengths_match(tmp_path):
    fasta, gtf, bam_path = reference_fixture(tmp_path)
    with pysam.AlignmentFile(bam_path) as source:
        header = source.header.to_dict()
    header["SQ"][0]["M5"] = "0" * 32
    wrong = tmp_path / "wrong.bam"
    with pysam.AlignmentFile(wrong, "wb", header=header):
        pass
    with gdc_genome(genome_fasta=fasta, gtf=gtf, annotation_version="test", cache_directory=tmp_path / "cache") as genome:
        with pysam.AlignmentFile(wrong) as bam, pytest.raises(ValueError, match="checksum mismatch"):
            validate_alignment_reference(bam, genome, contigs=["chr1"])


def test_multiple_declared_samples_are_not_pooled(tmp_path):
    header = HEADER.to_dict()
    header["RG"] = [dict(ID="a", SM="sample-A"), dict(ID="b", SM="sample-B")]
    path = tmp_path / "mixed.bam"
    with pysam.AlignmentFile(path, "wb", header=header):
        pass
    with pysam.AlignmentFile(path) as bam, pytest.raises(ValueError, match="multiple samples"):
        reconstruct_antigen_isoforms(bam, targets=target_models()[0], sample_id="pool", source="mixed")


def test_unresolved_partial_allele_blocks_dominance(tmp_path):
    targets, paths = target_models()
    r = targets[0].reference
    altered = r.sequence[:75] + ("A" if r.sequence[75] != "A" else "T") + r.sequence[76:120]
    reads = aligned("complete", r.sequence, paths[0]) + aligned("new-partial", altered, paths[0][:120])
    result = reconstruct(tmp_path, reads)
    assert any(not row["compatible_models"] for row in result["partial_sequences"])
    add_isoform_abundances(result, estimates(tmp_path, result, [0.9, 0.1]), provenance={"fixture": "partial-allele"})
    assert result["quantification"]["genes"][0]["dominant_isoform_id"] is None


def test_published_gdc_dictionary_keeps_human_and_viral_contigs():
    manifest = gdc_reference_dictionary()
    assert manifest["source_md5"] == "f64be73587a7f376c0d8353f1636dca7"
    assert len(manifest["contigs"]) == 2779
    assert next(r for r in manifest["contigs"] if r["name"] == "chr1")["length"] == 248956422
    assert {"chrEBV", "HPV16", "HPV18"}.issubset({row["name"] for row in manifest["contigs"]})
