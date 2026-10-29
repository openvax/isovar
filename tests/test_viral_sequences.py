"""Viral recovery must work without inventing transcripts or coding models."""

from hashlib import md5
import gzip
import json
from pathlib import Path
import random

import pysam
import pytest

from isovar import gdc_viral_panel, reconstruct_viral_sequences, write_viral_fasta
from isovar.antigen_reference import GDC_REFERENCE, gdc_reference_dictionary
from isovar.cli.commands import run


@pytest.fixture
def inputs(tmp_path):
    rng = random.Random(465)
    dna = {name: "".join(rng.choice("ACGT") for _ in range(length)) for name, length in
           (("chr1", 1400), ("HPV16", 900), ("chrEBV", 1100), ("HPV18", 800))}
    fasta = tmp_path / "reference.fa"
    fasta.write_text("".join(">%s\n%s\n" % item for item in dna.items()))
    header = pysam.AlignmentHeader.from_dict(dict(
        SQ=[dict(SN=n, LN=len(s), M5=md5(s.encode()).hexdigest()) for n, s in dna.items()],
        RG=[dict(ID="a", SM="tumor"), dict(ID="b", SM="tumor")]))
    return dna, fasta, header


def read(inputs, name, contig, start, length=100, *, sequence=None, cigar=None, flag=0, rg="a", mapq=60):
    dna, _, header = inputs
    r = pysam.AlignedSegment(header)
    r.query_name, r.flag = name, flag
    r.reference_id, r.reference_start = header.get_tid(contig), start
    r.cigarstring = cigar or "%dM" % length
    r.query_sequence = sequence if sequence is not None else dna[contig][start:start + length]
    r.query_qualities = [30] * len(r.query_sequence)
    r.mapping_quality = mapq
    r.set_tag("RG", rg)
    return r


def reconstruct(tmp_path, inputs, reads, **kwargs):
    _, fasta, header = inputs
    bam_path = tmp_path / "rna.bam"
    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for r in reads:
            bam.write(r)
    with pysam.AlignmentFile(bam_path) as bam:
        return reconstruct_viral_sequences(
            bam, genome_fasta=fasta, sample_id="tumor", source="RNA-1",
            contigs=kwargs.pop("contigs", ["HPV16", "chrEBV", "HPV18"]),
            reference_profile=kwargs.pop("reference_profile", None),
            reference_cache=str(tmp_path / "cache"), scratch_dir=tmp_path, **kwargs)


def test_published_panel_is_explicit_and_reconciles_all_dictionary_entries():
    panel = gdc_viral_panel()
    expected = {r["name"]: r for r in gdc_reference_dictionary()["contigs"]}
    assert panel["source_sha256"] == "f221e758f47e4784714b4a89198b89536ddacdc905da9c778a6802d3e1e51c56"
    assert len(panel["contigs"]) == len({r["name"] for r in panel["contigs"]}) == 200
    assert {"chrEBV", "HPV16", "HPV18", "HBV", "HCV-1", "KSHV", "MCV"} <= {r["name"] for r in panel["contigs"]}
    assert "chr1" not in {r["name"] for r in panel["contigs"]}
    for row in panel["contigs"]:
        assert (row["length"], row["md5"]) == (expected[row["name"]]["length"], expected[row["name"]]["md5"])


def test_custom_fasta_alias_reaches_coding_reference_lookup(tmp_path, inputs):
    dna, fasta, _ = inputs
    # Random synthetic sequence: this tests names, not a viral coding model.
    fasta.write_text(fasta.read_text().replace(">chrEBV\n", ">EBV\n"))
    result = reconstruct(tmp_path, inputs, [read(inputs, "observed", "chrEBV", 100)], contigs=["chrEBV"])
    row, = result["sequences"]
    assert row["sequence"] == dna["chrEBV"][100:200]
    assert row["compatible_fragments"] == 1
    coding_reference, = result["coding_reference"]["contigs"]
    assert coding_reference["contig"] == "chrEBV"
    assert coding_reference["status"] == "unavailable"
    assert result["reference"]["checked_sequence_contigs"] == ["chrEBV"]


def test_real_ebv_reference_matches_gdc_and_recovers_exact_observed_rna(tmp_path):
    fixture = Path(__file__).parent / "data/viral/ebv-NC_007605.1.fa.gz"
    text = gzip.decompress(fixture.read_bytes()).decode()
    sequence = "".join(text.splitlines()[1:]).upper()
    ebv = next(r for r in gdc_viral_panel()["contigs"] if r["name"] == "chrEBV")
    assert len(sequence) == ebv["length"] and md5(sequence.encode()).hexdigest() == ebv["md5"]
    fasta = tmp_path / "ebv.fa"
    fasta.write_text(">chrEBV\n" + sequence + "\n")
    header = pysam.AlignmentHeader.from_dict(dict(SQ=[dict(SN="chrEBV", LN=len(sequence), M5=ebv["md5"])],
                                                 RG=[dict(ID="a", SM="tumor")]))
    inputs = {"chrEBV": sequence}, fasta, header
    result = reconstruct(tmp_path, inputs, [read(inputs, "left", "chrEBV", 90000),
                                          read(inputs, "right", "chrEBV", 90050)], contigs=["chrEBV"])
    row, = result["sequences"]
    assert row["sequence"] == sequence[90000:90150]
    assert row["coverage"] == [1] * 50 + [2] * 50 + [1] * 50
    assert row["compatible_fragments"] == 2


def test_overlaps_recover_observed_sequence_without_gtf_gap_fill_or_duplicate_orientations(tmp_path, inputs):
    dna, _, _ = inputs
    reads = [read(inputs, "r1", "HPV16", 100), read(inputs, "r2", "HPV16", 160, flag=16),
             read(inputs, "distant", "HPV16", 500), read(inputs, "background", "chr1", 100)]
    result = reconstruct(tmp_path, inputs, reads)
    assert {r["sequence"] for r in result["sequences"]} == {dna["HPV16"][100:260], dna["HPV16"][500:600]}
    assert all(not r["coding"] and r["rna_strand"] == r["transcript_completeness"] == "unknown"
               for r in result["sequences"])
    assert all(all(r["coverage"]) for r in result["sequences"])
    assert result["support"]["target_fragments"] == 3
    assert sum(c["fragments"] for c in result["support"]["equivalence_classes"]) == 3
    counts = {r["contig"]: r for r in result["contigs"]}
    assert counts["HPV16"]["aligned_fragments"] == counts["HPV16"]["unambiguous_contig_fragments"] == 3
    assert counts["chrEBV"]["status"] == "no_eligible_alignments"
    assert result["reference"]["annotation"]["name"] == "unannotated"
    output = tmp_path / "viral.fa"
    write_viral_fasta(result, output)
    with pysam.FastxFile(str(output)) as records:
        assert {r.sequence for r in records} == {r["sequence"] for r in result["sequences"]}


def test_splice_and_sequence_alleles_remain_distinct(tmp_path, inputs):
    dna, _, _ = inputs
    spliced = dna["HPV16"][100:160] + dna["HPV16"][300:360]
    variant = spliced[:40] + ("A" if spliced[40] != "A" else "C") + spliced[41:]
    result = reconstruct(tmp_path, inputs, [
        read(inputs, "splice", "HPV16", 100, sequence=spliced, cigar="60M140N60M"),
        read(inputs, "variant", "HPV16", 100, sequence=variant, cigar="60M140N60M"),
        read(inputs, "retained", "HPV16", 100, length=260)])
    assert {r["sequence"] for r in result["sequences"]} == {spliced, variant, dna["HPV16"][100:360]}
    for row in result["sequences"]:
        assert row["whole_path_segments"]
        assert row["unique_fragments"] == 1


def test_contained_seed_cannot_truncate_competing_sequence_at_five_prime_end(tmp_path, inputs):
    reference = "T" * 50 + "A" * 50 + "C" * 50
    allele = reference[:120] + "G" + reference[121:]
    result = reconstruct(tmp_path, inputs, [
        read(inputs, "complete-allele", "HPV16", 100, sequence=allele, length=150),
        read(inputs, "left", "HPV16", 100, sequence=reference[:100]),
        read(inputs, "right", "HPV16", 150, sequence=reference[50:])])
    assert {row["sequence"] for row in result["sequences"]} == {reference, allele}
    assert result["support"]["target_fragments"] == 3
    assert sum(row["fragments"] for row in result["support"]["equivalence_classes"]) == 3
    assert all(row["unique_fragments"] == 1 for row in result["sequences"])


def test_discovery_is_per_contig_and_later_reads_still_count(tmp_path, inputs):
    reads = [read(inputs, "hpv%d" % i, "HPV16", 100) for i in range(8)]
    reads.append(read(inputs, "ebv", "chrEBV", 200))
    result = reconstruct(tmp_path, inputs, reads, max_discovery_reads=1)
    assert {c for row in result["sequences"] for c in row["viral_contigs"]} == {"HPV16", "chrEBV"}
    counts = {r["contig"]: r for r in result["contigs"]}
    assert counts["HPV16"]["aligned_fragments"] == counts["HPV16"]["sequence_assigned_fragments"] == 8
    assert counts["HPV16"]["status"] == "discovery_truncated"
    assert counts["chrEBV"]["status"] == "sequences_reconstructed"
    assert result["support"]["target_fragments"] == 9


@pytest.mark.parametrize("other", ["chr1", "chrEBV"])
def test_competing_host_or_viral_placements_never_become_unique(tmp_path, inputs, other):
    first = read(inputs, "ambiguous", "HPV16", 100)
    alternate = read(inputs, "ambiguous", other, 100, flag=256, sequence=first.query_sequence)
    result = reconstruct(tmp_path, inputs, [first, alternate])
    counts = {r["contig"]: r for r in result["contigs"]}
    assert counts["HPV16"]["aligned_fragments"] == 1
    assert counts["HPV16"]["unambiguous_contig_fragments"] == 0
    assert all(r["unique_fragments"] == 0 for r in result["sequences"])
    assert result["support"]["target_fragments"] == 1
    witness, = result["support"]["witnesses"]
    assert len(witness["records"]) == 2
    assert {r["contig"] for r in witness["records"]} == {"HPV16", other}


def test_contig_assignment_can_be_unambiguous_while_sequence_placement_is_not(tmp_path, inputs):
    first = read(inputs, "same-virus", "HPV16", 100)
    alternate = read(inputs, "same-virus", "HPV16", 400, flag=256, sequence=first.query_sequence)
    for r in (first, alternate):
        r.set_tag("NH", 2)
    result = reconstruct(tmp_path, inputs, [first, alternate])
    count = next(c for c in result["contigs"] if c["contig"] == "HPV16")
    assert count["unambiguous_contig_fragments"] == 1
    assert len(result["sequences"]) == 2 and all(row["unique_fragments"] == 0 for row in result["sequences"])


@pytest.mark.parametrize("mapq,nh", [(0, None), (255, None), (60, 2)])
def test_missing_placement_information_remains_ambiguous(tmp_path, inputs, mapq, nh):
    r = read(inputs, "unknown", "HPV16", 100, mapq=mapq)
    if nh:
        r.set_tag("NH", nh)
    result = reconstruct(tmp_path, inputs, [r])
    assert result["sequences"][0]["unique_fragments"] == 0
    assert next(c for c in result["contigs"] if c["contig"] == "HPV16")["unambiguous_contig_fragments"] == 0


@pytest.mark.parametrize("declaration,reason", [
    ("chr1,201,+,100M,60,0;", "missing_SA_placement"),
    ("not-an-SA-entry", "unparseable_SA_tag"),
])
def test_missing_supplementary_records_retain_sequence_without_unique_attribution(tmp_path, inputs, declaration, reason):
    r = read(inputs, "missing-piece", "HPV16", 100)
    r.set_tag("SA", declaration)
    result = reconstruct(tmp_path, inputs, [r])
    assert result["sequences"][0]["compatible_fragments"] == 1
    assert result["sequences"][0]["unique_fragments"] == 0
    assert reason in result["support"]["witnesses"][0]["placement_uncertainty"]


def test_missing_mapped_mate_stays_ambiguous(tmp_path, inputs):
    result = reconstruct(tmp_path, inputs, [read(inputs, "missing-mate", "HPV16", 100, flag=99)])
    assert result["sequences"][0]["unique_fragments"] == 0
    assert "missing_mapped_mate" in result["support"]["witnesses"][0]["placement_uncertainty"]


def test_actual_host_viral_split_keeps_both_sequence_pieces_and_witnesses(tmp_path, inputs):
    dna, _, _ = inputs
    sequence = dna["HPV16"][100:160] + dna["chr1"][200:260]
    viral = read(inputs, "split", "HPV16", 100, sequence=sequence, cigar="60M60S")
    host = read(inputs, "split", "chr1", 200, sequence=sequence, cigar="60S60M", flag=2048)
    viral.set_tag("SA", "chr1,201,+,60S60M,60,0;")
    host.set_tag("SA", "HPV16,101,+,60M60S,60,0;")
    result = reconstruct(tmp_path, inputs, [viral, host])
    row, = result["sequences"]
    assert row["sequence"] == sequence and row["other_contigs"] == ["chr1"]
    witness, = result["support"]["witnesses"]
    assert not witness["placement_uncertainty"]
    assert {r["contig"] for r in witness["records"]} == {"HPV16", "chr1"}


def test_original_mates_read_groups_duplicate_payloads_and_metadata(tmp_path, inputs):
    r1 = read(inputs, "pair", "HPV16", 100, flag=99)
    r2 = read(inputs, "pair", "HPV16", 140, flag=147)
    other = read(inputs, "pair", "HPV16", 100, rg="b")
    other.query_qualities = None
    other.set_tag("CB", "cell-1")
    other.set_tag("UB", "umi-1")
    result = reconstruct(tmp_path, inputs, [r1, r2, other, other])
    row, = result["sequences"]
    assert row["compatible_fragments"] == 2
    assert row["supporting_segments"] == 3
    assert result["input"]["records"] == 4 and result["input"]["unique_records"] == 3
    assert result["support"]["target_fragments"] == 2
    count = next(c for c in result["contigs"] if c["contig"] == "HPV16")
    assert count["aligned_fragments"] == 2 and count["aligned_segments"] == 3
    witness = next(w for w in result["support"]["witnesses"] if w["fragment"][0] == "b")
    assert witness["records"][0]["tags"]["CB"] == "cell-1"
    assert witness["records"][0]["base_qualities_available"] is False


def test_unassemblable_read_is_reported_as_aligned_not_as_zero_support(tmp_path, inputs):
    result = reconstruct(tmp_path, inputs, [read(inputs, "ambiguous-bases", "HPV16", 100, sequence="N" * 100)])
    assert not result["sequences"]
    count = next(c for c in result["contigs"] if c["contig"] == "HPV16")
    assert count["status"] == "aligned_but_not_reconstructed"
    assert count["aligned_fragments"] == count["unassigned_fragments"] == 1
    assert result["support"]["observation_rejections"]["ambiguous_bases"] > 0


def test_empty_input_and_excluded_reads_are_explicit(tmp_path, inputs):
    result = reconstruct(tmp_path, inputs, [read(inputs, "qc-fail", "HPV16", 100, flag=512)])
    assert not result["sequences"] and not result["support"]["target_fragments"]
    assert result["input"]["excluded_records"]
    assert all(row["status"] == "no_eligible_alignments" for row in result["contigs"])


def test_equal_length_wrong_reference_is_rejected(tmp_path, inputs):
    _, fasta, _ = inputs
    fasta.write_text(fasta.read_text().replace("A", "T"))
    with pytest.raises(ValueError, match="checksum mismatch"):
        reconstruct(tmp_path, inputs, [])


def test_gdc_profile_does_not_accept_custom_or_human_only_alignments(tmp_path, inputs):
    with pytest.raises(ValueError, match="dictionary does not match"):
        reconstruct(tmp_path, inputs, [], reference_profile=GDC_REFERENCE)
    with pytest.raises(ValueError, match="missing from the BAM"):
        reconstruct(tmp_path, inputs, [], contigs=["CMV"])
    with pytest.raises(ValueError, match="not in GDC"):
        reconstruct(tmp_path, inputs, [], reference_profile=GDC_REFERENCE, contigs=["chr1"])


def test_default_gdc_selection_includes_ebv_and_all_declared_viral_contigs(tmp_path, inputs, monkeypatch):
    import isovar.antigen_reference as reference
    import isovar.viral_sequences as viral
    dna, _, _ = inputs
    dictionary = [dict(name=n, length=len(s), md5=md5(s.encode()).hexdigest()) for n, s in dna.items()]
    panel = dict(profile=GDC_REFERENCE, source="offline fixture", source_sha256="fixture",
                 contigs=[r for r in dictionary if r["name"] != "chr1"])
    monkeypatch.setattr(reference, "gdc_reference_dictionary", lambda: dict(contigs=dictionary))
    monkeypatch.setattr(viral, "gdc_viral_panel", lambda: panel)
    result = reconstruct(tmp_path, inputs, [read(inputs, "ebv", "chrEBV", 100)],
                         reference_profile=GDC_REFERENCE, contigs=None)
    assert result["selected_contigs"] == ["HPV16", "HPV18", "chrEBV"]
    assert result["reference"]["checked_sequence_contigs"] == result["selected_contigs"]
    assert result["sequences"][0]["viral_contigs"] == ["chrEBV"]


def test_installed_command_shape_exports_without_any_gene_or_gtf(tmp_path, inputs):
    _, fasta, _ = inputs
    reconstruct(tmp_path, inputs, [read(inputs, "viral", "HPV16", 100)])
    output = tmp_path / "result"
    args = ["viral", "--bam", str(tmp_path / "rna.bam"), "--genome-fasta", str(fasta),
            "--sample-id", "tumor", "--contig", "HPV16", "--alignment-reference", "custom",
            "--reference-cache", str(tmp_path / "cache"), "--output-dir", str(output)]
    run(args)
    report = json.loads((output / "viral.json").read_text())
    assert report["schema"] == "isovar.viral_sequences.v1"
    assert report["input"]["sha256"]
    assert (output / "sequences.fa").read_text().count(">") == 1
    with pytest.raises(SystemExit) as error:
        run(args)
    assert error.value.code == 2
