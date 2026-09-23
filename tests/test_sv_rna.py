"""#305/#306: RNA paths, frames and event linkage for a nominated SV, from BAM records.

Synthetic cases place reads on explicit transcript models; no genome sequence
is needed because reconstruction never reads reference bases. Real cases use
the pinned, unchanged original records of the fusion corpora.
"""
from dataclasses import replace
import gzip
from hashlib import sha256
import json
from pathlib import Path
import random
import time

import pysam
import pytest

from isovar import reconstruct_sv_rna
from isovar.cli import commands
from isovar.fusion import FusionBreakpoint, FusionReference, fusion_from_dict, reconstruct_fusion
from isovar.read_collector import ReadCollector
from isovar.read_end_inference import Adapter, ReadEndProfile, reverse_complement
from isovar.sv_rna import RnaObservation, _extend_both, _ObservationIndex, sv_rna_input_from_dict
from tests.data.osteosarc.expansion.references import translate

HEADER = pysam.AlignmentHeader.from_dict(dict(
    SQ=[dict(SN="1", LN=20000), dict(SN="2", LN=20000)], RG=[dict(ID="a"), dict(ID="b")]))
STOPS = {"TAA", "TAG", "TGA"}
FUSIONS = Path(__file__).parent / "data/fusions"


def antisense(sequence, placements):
    return reverse_complement(sequence), [(c, p, "-" if s == "+" else "+") for c, p, s in reversed(placements)]


def record(name, contig, start, cigar, sequence, flag=0, group="a", qualities=True, mapq=60):
    read = pysam.AlignedSegment(HEADER)
    read.query_name = name
    read.reference_id = HEADER.get_tid(contig)
    read.reference_start = start
    read.cigarstring = cigar
    read.flag = flag
    read.query_sequence = sequence
    read.query_qualities = [30] * len(sequence) if qualities else None
    read.mapping_quality = mapq
    read.set_tag("RG", group)
    return read


def link(records):
    for read in records:
        read.set_tag("SA", "".join("%s,%d,%s,%s,%d,0;" % (
            other.reference_name, other.reference_start + 1, "-" if other.is_reverse else "+",
            other.cigarstring.replace("H", "S"), other.mapping_quality) for other in records if other is not read))
    return records


def aligned(name, sequence, placements, flag=0, group="a", qualities=True):
    """Records for one segment: a primary and hard-clipped supplementary pieces."""
    pieces = []
    for q, (contig, position, strand) in enumerate(placements):
        if pieces and pieces[-1][0] == (contig, strand):
            pieces[-1][1].append((q, position))
        else:
            pieces.append(((contig, strand), [(q, position)]))
    primary = max(range(len(pieces)), key=lambda k: len(pieces[k][1]))
    n, records = len(sequence), []
    for k, ((contig, strand), bases) in enumerate(pieces):
        q0, q1 = bases[0][0], bases[-1][0] + 1
        reverse = strand == "-"
        genomic = sorted(p for _, p in bases)
        operations = []
        for previous, position in zip([None] + genomic, genomic):
            if previous is not None and position > previous + 1:
                operations.append([position - previous - 1, "N"])
            if operations and operations[-1][1] == "M" and previous == position - 1:
                operations[-1][0] += 1
            else:
                operations.append([1, "M"])
        left, right = (n - q1, q0) if reverse else (q0, n - q1)
        clip = "S" if k == primary else "H"
        cigar = "".join("%d%s" % tuple(op) for op in
                        ([[left, clip]] if left else []) + operations + ([[right, clip]] if right else []))
        stored = reverse_complement(sequence) if reverse else sequence
        if k != primary:
            stored = stored[left:n - right]
        records.append(record(name, contig, genomic[0], cigar, stored,
                              flag | (16 if reverse else 0) | (0 if k == primary else 2048), group, qualities))
    return link(records) if len(records) > 1 else records


def write_bam(path, records, header=HEADER):
    unsorted = str(path) + ".unsorted"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for read in records:
            out.write(read)
    pysam.sort("-o", str(path), unsorted)
    pysam.index(str(path))
    return path


class Scenario:
    """Donor D-201 (chr 1, +) exons 1-2 spliced to acceptor A-201 (chr 2) exons 2-3."""

    def __init__(self, acceptor_strand="+", seed=0):
        self.rng = random.Random(seed)
        self.donor_ref, self.donor_positions = self.gene(
            "D-201", "1", "+", [(1000, 1100), (2000, 2150), (3000, 3300)], 20, 150)
        exons = [(5000, 5100), (6000, 6300), (7000, 7400)] if acceptor_strand == "+" else [
            (5000, 5400), (6000, 6300), (7000, 7100)]
        self.acceptor_ref, self.acceptor_positions = self.gene("A-201", "2", acceptor_strand, exons, 30, 200)
        self.sequence = self.donor_ref.sequence[:250] + self.acceptor_ref.sequence[100:]
        self.positions = self.donor_positions[:250] + self.acceptor_positions[100:]
        # DNA breakpoints in donor intron 2 and acceptor intron 1: splicing
        # removes both from the RNA.
        self.donor = FusionBreakpoint("1", 2500, "+")
        self.acceptor = FusionBreakpoint("2", 5500, "+") if acceptor_strand == "+" else FusionBreakpoint("2", 6500, "-")
        # The literal adjacency at the RNA junction.
        self.exact = (FusionBreakpoint("1", 2150, "+"),
                      FusionBreakpoint("2", 6000, "+") if acceptor_strand == "+" else FusionBreakpoint("2", 6300, "-"))

    def random(self, n):
        return "".join(self.rng.choice("ACGT") for _ in range(n))

    def gene(self, name, contig, strand, exons, cds_start, codons):
        coding = "ATG"
        while len(coding) < 3 * codons + 3:
            codon = self.random(3)
            coding += "" if codon in STOPS else codon
        sequence = self.random(cds_start) + coding + "TAA"
        sequence += self.random(sum(b - a for a, b in exons) - len(sequence))
        positions = [(contig, p, strand) for a, b in exons for p in range(a, b)]
        positions = positions[::-1] if strand == "-" else positions
        return FusionReference(name, "test", "synthetic 1", contig, strand, exons, sequence,
                               cds_start, cds_start + len(coding) + 3), positions

    def tile(self, prefix, sequence, positions, length=100, step=10, errors=None, **kwargs):
        """Alternate sense/antisense reads, as from an unstranded library."""
        records = []
        for i, s in enumerate(range(0, len(sequence) - length + 1, step)):
            read = sequence[s:s + length]
            if errors is not None:
                rng, rate = errors
                read = "".join(rng.choice([b for b in "ACGT" if b != x]) if rng.random() < rate else x for x in read)
            read = (read, positions[s:s + length])
            records += aligned("%s%d" % (prefix, i), *(antisense(*read) if i % 2 else read), **kwargs)
        return records

    def reads(self, **kwargs):
        return (self.tile("f", self.sequence, self.positions, **kwargs)
                + self.tile("d", self.donor_ref.sequence, self.donor_positions, step=13)
                + self.tile("a", self.acceptor_ref.sequence, self.acceptor_positions, step=17))

    def run(self, bam, donor=None, acceptor=None, references=None, regions=(), **kwargs):
        with pysam.AlignmentFile(str(bam)) as alignments:
            return reconstruct_sv_rna(
                alignments, event_id="D--A", reference_name="test", donor=donor or self.donor,
                acceptor=acceptor or self.acceptor, regions=regions,
                references=references or [self.donor_ref, self.acceptor_ref], sample_id="sample",
                source=str(bam), event_provenance=dict(caller="synthetic DNA call"), **kwargs)


def spanning(result, relation=None):
    return [j for p in result["paths"] for j in p["junctions"] if relation in (None, j["relation"])]


@pytest.mark.parametrize("acceptor_strand", ["+", "-"])
def test_spliced_fusion_is_reconstructed_and_translated_in_the_donor_frame(tmp_path, acceptor_strand):
    s = Scenario(acceptor_strand)
    result = s.run(write_bam(tmp_path / "rna.bam", s.reads()))
    assert result["status"] == "event_linked_candidates" and not result["limitations"]
    path, = result["paths"]
    assert path["sequence"] == s.sequence  # Assembled through overlaps, no reference padding.
    assert path["unplaced_intervals"] == []
    junction, = path["junctions"]
    assert junction["query_interval"] == [249, 250] and junction["kinds"] == ["split"]
    assert junction["relation"] == path["event_linkage"]["status"] == "event_compatible_junction"
    assert junction["breakpoint_assignment"] == [350, 500 if acceptor_strand == "+" else 200]
    assert not junction["annotated"] and not junction["forward_splice_geometry"]
    starts = range(0, len(s.sequence) - 99, 10)
    assert junction["direct_fragments"] == sum(a <= 249 and 250 < a + 100 for a in starts)
    assert not path["event_linkage"]["somatic_causation_proven"]

    assert path["frame_status"] == "translated"
    translation, = path["translations"]
    expected, stop = translate(s.sequence[20:], annotated_start=True)  # Independent NCBI table 1.
    assert (translation["amino_acids"], translation["ends_with_stop_codon"]) == (expected, stop)
    assert translation["translation_start"] == 20 and translation["transcript_ids"] == ["D-201"]
    evidence, = translation["frame_evidence"]
    assert evidence["complete_5prime"] and not evidence["upstream_frame_assumed"]
    assert evidence["departure"] == 250 and evidence["departure_relation"] == "event_compatible_junction"
    assert path["reference_readings"] == ["A-201"]  # The acceptor's own frame is not a novel protein.
    reference_proteins = [translate(r.sequence[r.cds_start:], annotated_start=True)[0]
                          for r in (s.donor_ref, s.acceptor_ref)]
    junction_codon = (250 - 20) // 3
    for peptide in translation["candidate_peptides"]:
        a, b = peptide["protein_interval"]
        assert expected[a:b] == peptide["sequence"] and b > junction_codon
        assert not any(peptide["sequence"] in protein for protein in reference_proteins)
    assert any(a <= junction_codon < b for a, b in (p["protein_interval"] for p in translation["candidate_peptides"]))


def test_exact_breakpoint_and_observed_junction_homology(tmp_path):
    s = Scenario()
    read = s.sequence[200:300]  # Junction after read base 49.
    records = []
    for name in ("h1", "h2"):
        # Both pieces align read bases 47-49: donor 2147-2149 and acceptor 5997-5999.
        records += link([record(name, "1", 2100, "50M50S", read),
                         record(name, "2", 5997, "47H53M", read[47:], 2048)])
    bam = write_bam(tmp_path / "rna.bam", records)
    for donor, acceptor, expected in [
            (2150, 6000, ("breakpoint_junction", [3, 0])), (2147, 5997, ("breakpoint_junction", [0, 3])),
            (2151, 6000, ("event_compatible_junction", [4, 0]))]:
        result = s.run(bam, FusionBreakpoint("1", donor, "+"), FusionBreakpoint("2", acceptor, "+"))
        path, = result["paths"]
        junction, = path["junctions"]
        assert path["sequence"] == read and junction["query_interval"] == [46, 50]
        assert junction["unplaced_bases"] == read[47:50]  # Placed twice, so by neither piece.
        assert (junction["relation"], junction["breakpoint_assignment"]) == expected
        assert junction["direct_fragments"] == 2


def test_one_read_is_enough_to_report_an_event_junction(tmp_path):
    s = Scenario()
    records = aligned("only", s.sequence[200:300], s.positions[200:300])
    variant = s.sequence[200:249] + ("A" if s.sequence[249] != "A" else "C") + s.sequence[250:300]
    records += aligned("other", variant, s.positions[200:300])  # A base variant at the junction.
    result = s.run(write_bam(tmp_path / "rna.bam", records), *s.exact)
    assert sorted(p["sequence"] for p in result["paths"]) == sorted([s.sequence[200:300], variant])
    # Both reads make the breakpoint join; a base elsewhere does not split its support.
    assert all(p["junctions"][0]["direct_fragments"] == 2 for p in result["paths"])
    assert not any(seed.get("minor_variant_of_seeded_junction") for seed in result["seeds"])
    # Against twenty reads, the single-read variant of the same join is pruned.
    records += [r for i in range(20) for r in aligned("m%d" % i, s.sequence[200:300], s.positions[200:300])]
    result = s.run(write_bam(tmp_path / "many.bam", records), *s.exact)
    assert [p["sequence"] for p in result["paths"]] == [s.sequence[200:300]]
    assert [seed.get("minor_variant_of_seeded_junction", False) for seed in result["seeds"]] == [False, True]


def test_sa_tag_without_its_observed_record_creates_no_junction(tmp_path):
    s = Scenario()
    primaries = [r for r in s.tile("f", s.sequence, s.positions) if not r.is_supplementary]
    bam = write_bam(tmp_path / "rna.bam", primaries)
    result = s.run(bam)
    assert result["status"] == "no_candidate_paths" and result["paths"] == []
    assert result["segment_path_notes"]["SA_declared_piece_not_linked"] > 0
    # A soft clip at the nominated breakpoint is unplaced sequence, not a partner.
    result = s.run(bam, *s.exact, read_collector=ReadCollector(use_soft_clipped_bases=True))
    assert result["status"] == "event_linked_candidates"
    assert {p["event_linkage"]["status"] for p in result["paths"]} == {"breakpoint_clip_partner_unplaced"}
    path, = [p for p in result["paths"] if p["junctions"][0]["left"]]
    junction, = path["junctions"]
    assert junction["right"] is None and junction["left"] == ["1", 2149, "+"]
    assert [path["unplaced_intervals"][-1][1]] == [len(path["sequence"])]
    assert path["blocks"][-1]["contig"] == "1"  # No acceptor mapping is invented.
    translation, = path["translations"]
    assert translation["amino_acids"] in translate(s.sequence[20:], annotated_start=True)[0]


@pytest.mark.parametrize("side", ["donor", "acceptor"])
@pytest.mark.parametrize("flag", [0, 16])
def test_breakpoint_soft_clips_inside_terminal_hard_clips(tmp_path, side, flag):
    s = Scenario()
    if side == "donor":
        sequence = s.sequence[200:250] + "T" * 20
        read = record("clip", "1", 2100, "50M20S10H", sequence, flag)
    else:
        sequence = "T" * 20 + s.sequence[250:300]
        read = record("clip", "2", 6000, "10H20S50M", sequence, flag)
    bam = write_bam(tmp_path / "clip.bam", [read])
    assert s.run(bam, *s.exact)["paths"] == []
    result = s.run(bam, *s.exact, read_collector=ReadCollector(use_soft_clipped_bases=True))
    path, = result["paths"]
    assert path["sequence"] == sequence
    junction, = path["junctions"]
    assert junction["relation"] == "breakpoint_clip_partner_unplaced"
    assert junction["direct_fragments"] == 1
    assert junction["right" if side == "donor" else "left"] is None


def junction_reads(s, step=10, length=100):
    return ["f%d" % i for i, a in enumerate(range(0, len(s.sequence) - length + 1, step))
            if a < s.junction <= a + length - 1]


def test_distant_pieces_and_mates_are_retrieved_but_mates_are_not_joined(tmp_path):
    s = Scenario()
    s.junction = 250
    first = aligned("pair", s.sequence[100:200], s.positions[100:200], flag=1 | 64)
    second = aligned("pair", *antisense(s.sequence[300:400], s.positions[300:400]), flag=1 | 128)
    for read, mate in ((first[0], second[0]), (second[0], first[0])):
        read.next_reference_id, read.next_reference_start = mate.reference_id, mate.reference_start
    bam = write_bam(tmp_path / "rna.bam", first + second + s.tile("f", s.sequence, s.positions))
    # Only donor exons and a narrow breakpoint window are searched: partner
    # records are reached only through observed SA and mate locations.
    result = s.run(bam, references=[s.donor_ref], breakpoint_window=10)
    hops = result["acquisition"]["hop_queries"]
    assert hops and all(h["fetched"] and h["contig"] == "2" for h in hops)
    path, = result["paths"]
    # The retrieved mate overlaps junction reads, so it extends the path to its
    # own end; nothing past it was retrieved. Alone, the pair bridges nothing.
    assert path["sequence"] == s.sequence[:400] and path["blocks"][-1]["contig"] == "2"
    direct = {result["original_records"][rid].split("\t")[0] for key in path["junctions"][0]["direct_observations"]
              for rid in result["observations"][key]["records"]}
    assert direct == set(junction_reads(s))  # The discordant pair is not junction evidence.
    alone = s.run(write_bam(tmp_path / "pair.bam", first + second), references=[s.donor_ref], breakpoint_window=10)
    assert alone["acquisition"]["records"] == 2 and alone["paths"] == []


def test_alternative_placements_and_duplicate_records_add_no_support(tmp_path):
    s = Scenario()
    s.junction = 250
    reads = s.tile("f", s.sequence, s.positions)
    baseline = s.run(write_bam(tmp_path / "base.bam", reads))
    junction, = spanning(baseline)
    assert junction["direct_fragments"] == len(junction_reads(s))
    primaries = {r.query_name: r for r in reads if not r.is_supplementary}
    # Secondary placements elsewhere are other hypotheses for those segments,
    # never additional supporting pieces.
    secondary = [record(name, "2", 4600, "100M", primaries[name].get_forward_sequence(), 256)
                 for name in junction_reads(s)[:2]]
    other_group = [r.__copy__() for r in reads if r.query_name in junction_reads(s)[:3]]
    for read in other_group:
        read.set_tag("RG", "b")  # The same QNAME in another read group is another fragment.
    result = s.run(write_bam(tmp_path / "alt.bam", reads + reads[:6] + secondary + other_group))
    path, = result["paths"]
    assert path["sequence"] == baseline["paths"][0]["sequence"]
    assert path["junctions"][0]["direct_fragments"] == junction["direct_fragments"] + 3
    assert path["sequence_evidence"]["secondary_segments"] == 0
    dropped = s.run(write_bam(tmp_path / "drop.bam", reads + secondary),
                    read_collector=ReadCollector(use_secondary_alignments=False))
    assert dropped["excluded_records"] == {"secondary": 2}
    assert dropped["paths"][0]["sequence"] == path["sequence"]


def test_missing_qualities_are_kept_unless_policy_requires_them(tmp_path):
    s = Scenario()
    bam = write_bam(tmp_path / "rna.bam", s.tile("f", s.sequence, s.positions, qualities=False))
    result = s.run(bam)
    path, = result["paths"]
    assert path["sequence"] == s.sequence and path["sequence_evidence"]["missing_quality_segments"] > 0
    assert all(o["missing_qualities"] for o in result["observations"].values())
    strict = s.run(bam, read_collector=ReadCollector(use_reads_without_base_qualities=False))
    assert strict["paths"] == [] and set(strict["excluded_records"]) == {"qualities_unavailable"}


def test_assembly_extends_beyond_junction_reads_only_when_enabled(tmp_path):
    s = Scenario()
    bam = write_bam(tmp_path / "rna.bam", s.tile("f", s.sequence, s.positions, step=40))
    assembled, = s.run(bam)["paths"]
    direct, = s.run(bam, assemble=False)["paths"]
    assert assembled["sequence"] == s.sequence[:max(range(0, len(s.sequence) - 99, 40)) + 100]
    assert direct["sequence"] in s.sequence and len(direct["sequence"]) < len(s.sequence)
    junction, = direct["junctions"]
    assert junction["linked_interval"] == [0, len(direct["sequence"])]  # Every base co-observed with the join.
    assert assembled["junctions"][0]["linked_interval"][1] < len(s.sequence)
    assert assembled["sequence_evidence"]["assembly_phase"] == "hypothesis_not_proven_long_range_phase"
    assert assembled["reconstruction_scopes"] == ["regional", "seed_spanning"]
    assert direct["reconstruction_scopes"] == ["seed_spanning"]


def test_overlap_extension_still_stops_before_revisiting_a_placement():
    sequence = "ACGTAC"
    positions = tuple(("1", q, "+") for q in range(6))
    observation = RnaObservation(
        key="cycle", identity=("rg", "cycle", 0), records=("cycle",),
        sequence=sequence + sequence, positions=positions + positions, breaks=((5, 6, "split"),),
        reverse=False, query_interval=(0, 12), missing_qualities=False, secondary=False)
    notes = set()
    paths = list(_extend_both(sequence, positions, _ObservationIndex([observation]), 3,
                              (2, 0.1, 0.5), 8, [], notes))
    assert len(paths) == 1 and paths[0][:2] == (sequence, positions)
    assert notes == {"repeated_genomic_position"}


def test_deep_noisy_coverage_prunes_error_branches_without_path_explosion(tmp_path):
    s = Scenario()
    errors = (random.Random(5), 0.004)
    reads = (s.tile("f", s.sequence, s.positions, 150, 2, errors)
             + s.tile("d", s.donor_ref.sequence, s.donor_positions, 150, 3, errors)
             + s.tile("a", s.acceptor_ref.sequence, s.acceptor_positions, 150, 3, errors))
    bam = write_bam(tmp_path / "deep.bam", reads)
    started = time.perf_counter()
    result = s.run(bam)
    assert time.perf_counter() - started < 30
    assert "path_limit" not in result["limitations"]
    path, = result["paths"]
    assert path["sequence"] == s.sequence  # Consensus of observed bases, not one read.
    assert result["pruned_branches"]


def variant(s):
    return s.sequence[:600] + ("A" if s.sequence[600] != "A" else "C") + s.sequence[601:]


def test_undistinguished_alternatives_fork_and_strong_ones_prune(tmp_path):
    s = Scenario()
    reads = aligned("x", s.sequence, s.positions) + aligned("y", variant(s), s.positions)
    two = s.run(write_bam(tmp_path / "two.bam", reads))
    assert sorted(p["sequence"] for p in two["paths"]) == sorted([s.sequence, variant(s)])
    assert not two["pruned_branches"]  # One read each: neither can outweigh the other.
    many = [r for i in range(20) for r in aligned("m%d" % i, s.sequence, s.positions)] + reads[2:]
    pruned = s.run(write_bam(tmp_path / "many.bam", many))
    path, = pruned["paths"]
    assert path["sequence"] == s.sequence
    branch, = pruned["pruned_branches"]
    assert branch["fragments"] == 1 and branch["direction"] == "3prime"
    limited = s.run(tmp_path / "two.bam", max_paths=1)
    assert len(limited["paths"]) == 1 and "path_limit" in limited["limitations"]


def test_local_errors_are_resolved_by_majority_but_splice_choices_fork(tmp_path):
    s = Scenario()
    base = s.sequence[:600] + ("A" if s.sequence[600] != "A" else "C") + s.sequence[601:]
    reads = [r for i in range(7) for r in aligned("m%d" % i, s.sequence, s.positions)]
    reads += [r for i in range(3) for r in aligned("v%d" % i, base, s.positions)]  # 30%: like an ONT error.
    bam = write_bam(tmp_path / "local.bam", reads)
    majority = s.run(bam)
    assert [p["sequence"] for p in majority["paths"]] == [s.sequence]
    assert [b["fragments"] for b in majority["pruned_branches"]] == [3]
    both = s.run(bam, min_local_variant_fraction=0.3)
    assert sorted(p["sequence"] for p in both["paths"]) == sorted([s.sequence, base])
    # A different acceptor exon (a splice choice) at the same 30% still forks.
    skipped = s.sequence[:350] + s.sequence[650:]
    reads = [r for i in range(7) for r in aligned("m%d" % i, s.sequence, s.positions)]
    reads += [r for i in range(3) for r in aligned("x%d" % i, skipped, s.positions[:350] + s.positions[650:])]
    forks = s.run(write_bam(tmp_path / "splice.bam", reads))
    assert sorted(p["sequence"] for p in forks["paths"]) == sorted([s.sequence, skipped])


def test_short_cigar_n_gaps_are_deletions_not_introns(tmp_path):
    s = Scenario()
    positions = s.donor_positions[100:160] + s.donor_positions[165:250]  # A 5-base "intron" in exon 2.
    sequence = s.donor_ref.sequence[100:160] + s.donor_ref.sequence[165:250]
    reads = [r for i in range(3) for r in aligned("n%d" % i, sequence, positions)]
    assert reads[0].cigarstring == "60M5N85M"
    result = s.run(write_bam(tmp_path / "rna.bam", reads))
    assert result["status"] == "no_candidate_paths" and result["seeds"] == []


def test_noisy_reads_count_as_direct_support_for_their_own_junction(tmp_path):
    s = Scenario()
    reads = s.tile("f", s.sequence, s.positions, 150, 7, (random.Random(3), 0.06))  # ONT-like 6% errors.
    spans = sum(a < 250 <= a + 149 for a in range(0, len(s.sequence) - 149, 7))
    result = s.run(write_bam(tmp_path / "rna.bam", reads))
    junction, = [j for p in result["paths"] for j in p["junctions"] if j["relation"] == "event_compatible_junction"][:1]
    # Almost no read matches the consensus end to end, but each makes the join.
    assert junction["direct_fragments"] == junction["direct_segments"] == spans
    assert result["paths"][0]["sequence_evidence"]["voting_fragments"] > 0


def test_junction_bases_placed_past_the_breakpoint_are_the_same_adjacency(tmp_path):
    s = Scenario()
    donor, acceptor = s.exact  # Donor exon ends at 2150, acceptor exon starts at 6000.
    exact = (s.sequence[200:300], s.positions[200:300])
    shifted = s.positions[200:250] + [("1", p, "+") for p in (2150, 2151, 2152)] + s.positions[253:300]
    reads = [r for i in range(3) for r in aligned("e%d" % i, *exact)]
    reads += [r for i in range(2) for r in aligned("s%d" % i, exact[0], shifted)]  # Aligner kept 3 homologous bases.
    bam = write_bam(tmp_path / "rna.bam", reads)
    result = s.run(bam, donor, acceptor)
    junctions = {tuple(j["breakpoint_assignment"]): j for p in result["paths"] for j in p["junctions"]}
    assert set(junctions) <= {(0, 0), (-3, 3)} and all(j["relation"] == "breakpoint_junction" for j in junctions.values())
    assert all(j["direct_fragments"] == 5 for j in junctions.values())
    assert [s["minor_variant_of_seeded_junction"] for s in result["seeds"] if "minor_variant_of_seeded_junction" in s] \
        in ([], [True])
    # Beyond the shift tolerance it is a different, unsupported adjacency.
    narrow = s.run(bam, donor, acceptor, max_breakpoint_shift=2)
    assert {tuple(j["breakpoint_assignment"]) for p in narrow["paths"] for j in p["junctions"]
            if j["relation"] == "breakpoint_junction"} == {(0, 0)}
    assert narrow["paths"][0]["junctions"][0]["direct_fragments"] == 3


def test_breakpoint_join_between_annotated_splice_sites_is_read_through_ambiguous(tmp_path):
    s = Scenario()
    sequence = s.acceptor_ref.sequence
    positions = [("1", p, "+") for a, b in [(10000, 10100), (11000, 11300), (12000, 12400)] for p in range(a, b)]
    downstream = replace(s.acceptor_ref, transcript_id="B-201", contig="1",
                         exons=((10000, 10100), (11000, 11300), (12000, 12400)))
    fused = s.donor_ref.sequence[:250] + sequence[100:]
    placed = s.donor_positions[:250] + positions[100:]
    bam = write_bam(tmp_path / "rna.bam", s.tile("f", fused, placed))
    references = [s.donor_ref, downstream]
    # D exon 2 (ends 2150) spliced to B exon 2 (starts 11000), downstream on
    # the same strand: read-through splicing makes the same RNA.
    at_sites = s.run(bam, FusionBreakpoint("1", 2150, "+"), FusionBreakpoint("1", 11000, "+"), references=references)
    junction, = at_sites["paths"][0]["junctions"]
    assert junction["relation"] == "splice_ambiguous_event_junction" and junction["breakpoint_assignment"] == [0, 0]
    assert junction["forward_splice_geometry"] and at_sites["status"] == "splice_ambiguous_candidates"
    # Where B's annotated exon starts elsewhere, the join is not an annotated splice.
    moved = replace(downstream, exons=((10000, 10100), (11001, 11300), (12000, 12401)))
    other = s.run(bam, FusionBreakpoint("1", 2150, "+"), FusionBreakpoint("1", 11000, "+"),
                  references=[s.donor_ref, moved])
    assert other["paths"][0]["junctions"][0]["relation"] == "breakpoint_junction"
    assert other["status"] == "event_linked_candidates"


def test_unannotated_joins_near_annotated_ones_are_aligner_wobble(tmp_path):
    s = Scenario()
    for shift, seeded in ((3, False), (8, True)):
        positions = s.donor_positions[:100] + [("1", p + shift, "+") for _, p, _ in s.donor_positions[100:200]]
        reads = [r for i in range(3) for r in aligned("w%d" % i, s.donor_ref.sequence[:200], positions)]
        result = s.run(write_bam(tmp_path / ("w%d.bam" % shift), reads))
        assert bool(result["seeds"]) == seeded
        assert ("join_near_annotated_junction" in result["segment_path_notes"]) != seeded


def test_deep_long_read_locus_builds_only_what_paths_need(tmp_path):
    s = Scenario()
    reads = [r for i in range(4000) for r in aligned("d%d" % i, s.donor_ref.sequence, s.donor_positions)]
    reads += [r for i in range(4000) for r in aligned("a%d" % i, s.acceptor_ref.sequence, s.acceptor_positions)]
    reads += [r for i in range(12) for r in aligned("f%d" % i, s.sequence, s.positions)]
    started = time.perf_counter()
    result = s.run(write_bam(tmp_path / "deep.bam", reads), max_records=20000)
    assert time.perf_counter() - started < 60
    path, = result["paths"]
    assert path["sequence"] == s.sequence and path["junctions"][0]["direct_fragments"] == 12
    counts = result["observation_counts"]
    assert counts["eligible_segments"] == 8012 and counts["built_segments"] < 1000
    assert "extension_segment_limit" in result["limitations"]


def test_lazily_queued_joins_beyond_the_path_budget_are_reported_unexplored(tmp_path):
    s = Scenario()
    skip = s.donor_ref.sequence[:100] + s.donor_ref.sequence[250:]
    skip_positions = s.donor_positions[:100] + s.donor_positions[250:]
    reads = s.tile("f", s.sequence, s.positions) + s.tile("s", skip, skip_positions, step=40)
    result = s.run(write_bam(tmp_path / "rna.bam", reads), max_paths=1)
    assert len(result["paths"]) == 1 and "path_limit" in result["limitations"]
    unexplored, = [row for row in result["seeds"] if row.get("unexplored")]
    assert unexplored["relation"] == "regional_novel_junction" and unexplored["fragments"] == 2
    assert unexplored["left"] == ["1", 1099, "+"] and unexplored["paths"] == []
    # The queued join's reads were never built.
    assert result["observation_counts"]["built_segments"] < result["observation_counts"]["eligible_segments"]


def genome_base(s, position):
    """Transcript bases where annotated; a fixed pseudo-random base elsewhere."""
    if not hasattr(s, "bases"):
        s.bases = dict(zip(s.donor_positions, s.donor_ref.sequence))
        s.bases.update(zip(s.acceptor_positions, s.acceptor_ref.sequence))
    return s.bases.get(position) or random.Random("%s:%d" % position[:2]).choice("ACGT")


def breakpoint_read(s, name, d, a, flank=60):
    """A split read placed ``d`` / ``a`` bases from the exact adjacency (negative: past it)."""
    placements = ([("1", p, "+") for p in range(2150 - d - flank, 2150 - d)]
                  + [("2", p, "+") for p in range(6000 + a, 6000 + a + flank)])
    return aligned(name, "".join(genome_base(s, p) for p in placements), placements)


def test_review_regressions_for_adjacency_placements(tmp_path):
    s = Scenario()
    # A single split read at an intronic-breakpoint fusion still seeds a path.
    single = s.reads(step=10000) + aligned("f", s.sequence[200:300], s.positions[200:300])
    result = s.run(write_bam(tmp_path / "single.bam", single))
    assert result["status"] == "event_linked_candidates" and len(result["paths"]) == 1
    # Singleton near-misses are one event, not one path each.
    reads = [r for k, (d, a) in enumerate([(1, 2), (2, 1), (3, 3), (0, 4), (4, 0), (2, 2)])
             for r in breakpoint_read(s, "n%d" % k, d, a)]
    near = s.run(write_bam(tmp_path / "near.bam", reads), *s.exact)
    assert len(near["paths"]) == 1
    # The strongest placement is the reference: a weak exact join and a
    # 3-read near-miss are pruned against 20 reads at (2, 3).
    reads = [r for k in range(20) for r in breakpoint_read(s, "m%d" % k, 2, 3)]
    reads += [r for k in range(3) for r in breakpoint_read(s, "w%d" % k, 4, 4)] + breakpoint_read(s, "x", 0, 0)
    strongest = s.run(write_bam(tmp_path / "strong.bam", reads), *s.exact)
    assert len(strongest["paths"]) == 1
    assert sum(bool(row.get("minor_variant_of_seeded_junction")) for row in strongest["seeds"]) == 2
    # A join with one side just past a breakpoint is a near-miss of the event.
    past = s.run(write_bam(tmp_path / "past.bam", [r for k in range(4) for r in breakpoint_read(s, "p%d" % k, -2, 5)]),
                 *s.exact)
    assert past["status"] == "event_linked_candidates"
    assert past["paths"][0]["junctions"][0]["relation"] == "event_compatible_junction"
    # Reference bases from past both breakpoints are another adjacency.
    both = s.run(write_bam(tmp_path / "both.bam", [r for k in range(3) for r in breakpoint_read(s, "b%d" % k, -10, -10)]),
                 *s.exact)
    assert not any(j["relation"] == "breakpoint_junction" for p in both["paths"] for j in p["junctions"])


def test_breakpoint_join_with_inserted_bases_in_a_same_strand_event(tmp_path):
    s = Scenario()
    placements = [("1", p, "+") for p in range(2050, 2100)] + [("1", p, "+") for p in range(3052, 3102)]
    sequence = "".join(genome_base(s, p) for p in placements[:50]) + "TTT" + "".join(
        genome_base(s, p) for p in placements[50:])
    read = record("ins", "1", 2050, "50M3I952N50M", sequence)
    result = s.run(write_bam(tmp_path / "rna.bam", [read]), FusionBreakpoint("1", 2100, "+"),
                   FusionBreakpoint("1", 3050, "+"))
    junction, = result["paths"][0]["junctions"]
    assert junction["relation"] == "breakpoint_junction" and junction["breakpoint_assignment"] == [0, 2]
    assert junction["unplaced_bases"] == "TTT" and result["status"] == "event_linked_candidates"


def test_read_through_is_a_property_of_the_adjacency_not_of_one_placement(tmp_path):
    s = Scenario()
    positions = [("1", p, "+") for a, b in [(10000, 10100), (11000, 11300), (12000, 12400)] for p in range(a, b)]
    downstream = replace(s.acceptor_ref, transcript_id="B-201", contig="1",
                         exons=((10000, 10100), (11000, 11300), (12000, 12400)))
    exact = s.donor_positions[200:250] + positions[100:150]
    shifted = s.donor_positions[200:250] + [("1", 2150, "+"), ("1", 2151, "+")] + positions[102:150]
    sequence = s.donor_ref.sequence[200:250] + s.acceptor_ref.sequence[100:150]
    reads = [r for k in range(9) for r in aligned("e%d" % k, sequence, exact)]
    reads += [r for k in range(2) for r in aligned("s%d" % k, sequence, shifted)]
    result = s.run(write_bam(tmp_path / "rna.bam", reads), FusionBreakpoint("1", 2150, "+"),
                   FusionBreakpoint("1", 11000, "+"), references=[s.donor_ref, downstream])
    assert result["status"] == "splice_ambiguous_candidates"
    assert not any(j["relation"] == "breakpoint_junction" for p in result["paths"] for j in p["junctions"])


def test_multi_record_reads_share_a_lazily_seeded_join_with_single_records(tmp_path):
    s = Scenario()
    skip = s.donor_ref.sequence[:100] + s.donor_ref.sequence[250:]
    skip_positions = s.donor_positions[:100] + s.donor_positions[250:]
    reads = [r for k in range(3) for r in aligned("s%d" % k, skip[40:200], skip_positions[40:200])]
    reads.append(record("s0", "1", 2600, "40M", skip[300:340], 256))  # A secondary placement of s0.
    result = s.run(write_bam(tmp_path / "rna.bam", reads))
    seeded = [row for row in result["seeds"] if row["relation"] == "regional_novel_junction" and row["paths"]]
    assert [row["fragments"] for row in seeded] == [3] and len(result["paths"]) == 1


def test_extension_cap_keeps_reads_that_can_extend(tmp_path):
    s = Scenario()
    bam = write_bam(tmp_path / "rna.bam", s.tile("f", s.sequence, s.positions, 100, 2))
    for cap in (1000, 5):
        path, = s.run(bam, max_extension_segments=cap)["paths"]
        assert path["sequence"] == s.sequence


def test_one_long_read_does_not_decide_a_path_the_short_reads_contradict(tmp_path):
    s = Scenario()
    skipped = s.sequence[:450] + s.sequence[650:]
    long_read = aligned("long", skipped, s.positions[:450] + s.positions[650:])
    result = s.run(write_bam(tmp_path / "rna.bam", s.tile("f", s.sequence, s.positions, 100, 10) + long_read))
    sequences = [p["sequence"] for p in result["paths"]]
    assert s.sequence in sequences  # The majority isoform, not only the long read's.
    majority = result["paths"][sequences.index(s.sequence)]
    junction = majority["junctions"][0]
    # The long read makes the junction but skips path bases: it adds no linked interval.
    assert junction["linked_interval"][1] <= 350


@pytest.mark.parametrize("long_insertion", ["", "AAA", "TTT"])
def test_linked_interval_requires_the_observed_junction_insertion(tmp_path, long_insertion):
    s = Scenario()
    short_sequence = s.sequence[200:250] + "TTT" + s.sequence[250:300]
    reads = [r for k in range(7) for r in link([
        record("short%d" % k, "1", 2100, "50M3I50S", short_sequence),
        record("short%d" % k, "2", 6000, "53H50M", s.sequence[250:300], 2048)])]
    reads += s.tile("donor", s.donor_ref.sequence, s.donor_positions)
    reads += s.tile("acceptor", s.acceptor_ref.sequence, s.acceptor_positions)
    long_sequence = s.sequence[:250] + long_insertion + s.sequence[250:]
    reads += link([
        record("long", "1", 1000, "100M900N150M" + ("3I" if long_insertion else "") + "700S", long_sequence),
        record("long", "2", 6000, "%dH300M700N400M" % (250 + len(long_insertion)), s.sequence[250:], 2048)])
    result = s.run(write_bam(tmp_path / "insertion.bam", reads), *s.exact)
    path, = [p for p in result["paths"] if p["sequence"] == s.sequence[:250] + "TTT" + s.sequence[250:]]
    junction, = path["junctions"]
    assert junction["direct_fragments"] == 8  # All eight reads make the join.
    assert junction["linked_interval"] == ([0, 953] if long_insertion == "TTT" else [200, 303])


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("insertion_cigar", ["3I1900N", "1900N3I", "1I1900N2I"])
def test_unbuilt_junction_reads_preserve_inserted_sequence(tmp_path, strand, insertion_cigar):
    s = Scenario()
    sequence = s.donor_ref.sequence[40:100] + "TTT" + s.donor_ref.sequence[250:350]
    reads = [record("ins%03d" % k, "1", 1040, "10H7S60M" + insertion_cigar + "100M5S8H",
                    "G" * 7 + sequence + "C" * 5, flag=16 if k % 2 else 0) for k in range(250)]
    reference = s.donor_ref
    if strand == "-":
        reference = replace(reference, strand="-", sequence=reverse_complement(reference.sequence),
                            cds_start=None, cds_end=None)
    result = s.run(write_bam(tmp_path / "deep-insertion.bam", reads), references=[reference])
    junction, = spanning(result, "regional_novel_junction")
    assert junction["direct_segments"] == junction["direct_fragments"] == 250
    assert junction["direct_read_lineage"]["unresolved_segments"] == 250
    assert result["observation_counts"]["built_segments"] == 200
    assert "seed_segment_limit" in result["limitations"]
    assert junction["direct_junction_sequences"] == [["TTT" if strand == "+" else "AAA", 250]]


def test_reads_whose_build_fails_are_not_direct_support(tmp_path):
    s = Scenario()
    skip = s.donor_ref.sequence[:100] + s.donor_ref.sequence[250:]
    skip_positions = s.donor_positions[:100] + s.donor_positions[250:]
    reads = [r for k in range(3) for r in aligned("c%d" % k, skip[40:200], skip_positions[40:200])]
    ambiguous = skip[40:120] + "N" + skip[121:200]
    reads += [r for k in range(4) for r in aligned("n%d" % k, ambiguous, skip_positions[40:200])]
    result = s.run(write_bam(tmp_path / "rna.bam", reads))
    junction, = [j for p in result["paths"] for j in p["junctions"]]
    assert junction["direct_segments"] == len(junction["direct_observations"]) == 3
    assert result["segment_path_notes"]["ambiguous_bases"] == 4


def test_record_and_query_limits_are_reported(tmp_path):
    s = Scenario()
    bam = write_bam(tmp_path / "rna.bam", s.reads())
    assert "record_limit" in s.run(bam, max_records=5)["limitations"]
    assert "query_limit" in s.run(bam, max_queries=1)["limitations"]


def test_novel_noncoding_continuation_is_translated_in_the_donor_frame(tmp_path):
    s = Scenario()
    tail = s.random(400)
    sequence = s.donor_ref.sequence[:250] + tail
    positions = s.donor_positions[:250] + [("2", 12000 + i, "+") for i in range(400)]
    bam = write_bam(tmp_path / "rna.bam", s.tile("f", sequence, positions))
    result = s.run(bam, s.exact[0], FusionBreakpoint("2", 12000, "+"), references=[s.donor_ref])
    assert result["reference_models"] == [dict(
        transcript_id="D-201", annotation="synthetic 1", contig="1", strand="+",
        sequence_sha256=__import__("hashlib").sha256(s.donor_ref.sequence.encode()).hexdigest())]
    path, = result["paths"]
    assert path["sequence"] == sequence and path["event_linkage"]["status"] == "breakpoint_junction"
    translation, = path["translations"]
    assert (translation["amino_acids"], translation["ends_with_stop_codon"]) == translate(
        sequence[20:], annotated_start=True)
    assert result["acquisition"]["searched_regions"][-1] == ["2", 11000, 13000]


def test_unattributed_junctions_need_support_and_annotation_competes(tmp_path):
    s = Scenario()
    skip = s.donor_ref.sequence[:100] + s.donor_ref.sequence[250:]
    skip_positions = s.donor_positions[:100] + s.donor_positions[250:]
    one = s.run(write_bam(tmp_path / "one.bam", aligned("s0", skip[:200], skip_positions[:200])))
    assert one["status"] == "no_candidate_paths"
    two = s.run(write_bam(tmp_path / "two.bam", s.tile("s", skip, skip_positions, step=40)))
    assert two["status"] == "regional_candidates_only"
    path, = two["paths"]
    assert path["event_linkage"]["status"] == "regional_novel_junction"
    translation, = path["translations"]  # Exon 2 skipping is in frame: 50 codons removed.
    assert translation["amino_acids"] == translate(skip[20:], annotated_start=True)[0]
    assert translation["departure_relations"] == ["regional_novel_junction"]
    # The same exon skip across a nominated deletion is also ordinary splicing
    # of an annotated isoform: it competes, and is not a candidate.
    isoform = FusionReference("D-202", "test", "synthetic 1", "1", "+", [(1000, 1100), (3000, 3300)],
                              skip, 20, s.donor_ref.cds_end - 150)
    result = s.run(tmp_path / "two.bam", FusionBreakpoint("1", 1500, "+"), FusionBreakpoint("1", 2600, "+"),
                   references=[s.donor_ref, isoform, s.acceptor_ref])
    assert result["status"] == "no_candidate_paths"
    competing, = result["competing_annotated_junctions"]
    assert competing["left"] == ["1", 1099, "+"] and competing["right"] == ["1", 3000, "+"]
    unannotated = s.run(tmp_path / "two.bam", FusionBreakpoint("1", 1500, "+"), FusionBreakpoint("1", 2600, "+"))
    assert unannotated["status"] == "splice_ambiguous_candidates"
    assert unannotated["paths"][0]["event_linkage"]["status"] == "splice_ambiguous_event_junction"


def test_inputs_and_parameters_are_validated(tmp_path):
    s = Scenario()
    bam = write_bam(tmp_path / "rna.bam", [])
    with pytest.raises(ValueError, match="provenance"):
        reconstruct_sv_rna(None, event_id="e", reference_name="test", donor=s.donor, acceptor=s.acceptor,
                           regions=[], references=[s.donor_ref], sample_id="s", source="x", event_provenance={})
    for options in [dict(max_paths=0), dict(min_anchor_bases=2), dict(min_alternative_fraction=1.5),
                    dict(peptide_lengths=[]), dict(min_orf_amino_acids=0), dict(max_orf_candidates=0)]:
        with pytest.raises(ValueError):
            s.run(bam, **options)
    with pytest.raises(ValueError, match="interval"):
        s.run(bam, regions=[("3", 0, 10)])
    with pytest.raises(ValueError, match="versioned"):
        s.run(bam, references=[replace(s.donor_ref, reference_name="other")])


def test_cli_writes_the_api_result(tmp_path):
    s = Scenario()
    bam = write_bam(tmp_path / "rna.bam", s.tile("f", s.sequence, s.positions))
    data = dict(event_id="D--A", reference_name="test", sample_id="sample",
                donor=vars(s.donor), acceptor=vars(s.acceptor), regions=[["1", 2400, 2600]],
                references=[vars(r) for r in (s.donor_ref, s.acceptor_ref)],
                event_provenance=dict(caller="synthetic DNA call"))
    (tmp_path / "event.json").write_text(json.dumps(data))
    output = tmp_path / "result.json"
    commands.run(["sv-rna", "--bam", str(bam), "--input", str(tmp_path / "event.json"),
                  "--output", str(output), "--orf-output-prefix", str(tmp_path / "orfs"),
                  "--no-assembly", "--min-orf-amino-acids", "5", "--max-orf-candidates", "2"])
    result = json.loads(output.read_text())
    with pysam.AlignmentFile(str(bam)) as alignments:
        expected = reconstruct_sv_rna(alignments, source=str(bam), assemble=False,
                                      min_orf_amino_acids=5, max_orf_candidates=2,
                                      **sv_rna_input_from_dict(json.loads(json.dumps(data))))
    assert result == json.loads(json.dumps(expected))
    assert result["parameters"]["assemble"] is False and result["status"] == "event_linked_candidates"
    from isovar import export_sv_rna_orfs
    assert json.loads((tmp_path / "orfs.json").read_text()) == export_sv_rna_orfs(result)
    assert (tmp_path / "orfs.tsv").is_file()
    assert (tmp_path / "orfs.protein.fasta").is_file()
    assert (tmp_path / "orfs.nucleotide.fasta").is_file()


def corpus_bam(tmp_path, name, directory="coding-corpus"):
    data = json.loads(gzip.decompress((FUSIONS / directory / (name + ".input.json.gz")).read_bytes()))
    contigs = sorted({r["contig"] for r in data["references"]} | {data["fusion"]["donor"]["contig"],
                                                                   data["fusion"]["acceptor"]["contig"]})
    header = pysam.AlignmentHeader.from_dict(dict(SQ=[dict(SN=c, LN=300000000) for c in contigs]))
    sams = sorted({s for r in data["original_records"] for s in (r["sam"], r.get("partner_sam")) if s})
    unsorted, path = tmp_path / (name + ".unsorted"), tmp_path / (name + ".bam")
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for sam in sams:
            out.write(pysam.AlignedSegment.fromstring(sam, header))
    pysam.sort("-o", str(path), str(unsorted))
    pysam.index(str(path))
    fusion = data["fusion"]
    inputs = sv_rna_input_from_dict(dict(
        event_id=name, reference_name=fusion["reference_name"], sample_id=fusion["provenance"]["sample_id"],
        donor=fusion["donor"], acceptor=fusion["acceptor"], references=data["references"],
        event_provenance=dict(fixture=name, caller=fusion["provenance"].get("caller", "osteosarc catalogue"))))
    return data, path, inputs


def run_corpus(path, inputs, **kwargs):
    with pysam.AlignmentFile(str(path)) as bam:
        return reconstruct_sv_rna(bam, source=path.name, **inputs, **kwargs)


def test_original_short_reads_reproduce_the_validated_atp5mg_kmt2a_translation(tmp_path):
    data, bam, inputs = corpus_bam(tmp_path, "ATP5MG--KMT2A")
    supplied, = reconstruct_fusion(*fusion_from_dict(data))["translations"]
    result = run_corpus(bam, inputs)
    path, = result["paths"]
    junction, = path["junctions"]
    # Exactly at the nominated breakpoint, but also an annotated ATP5MG exon
    # end spliced to an annotated KMT2A exon start downstream on the same
    # strand: read-through splicing makes the same RNA.
    assert junction["relation"] == "splice_ambiguous_event_junction" and junction["breakpoint_assignment"] == [0, 0]
    assert junction["forward_splice_geometry"] and result["status"] == "splice_ambiguous_candidates"
    assert junction["kinds"] == ["N"] and junction["direct_fragments"] == 2  # STAR's CIGAR N, not SA.
    assert data["fusion"]["sequence"] in path["sequence"]
    translation, = path["translations"]
    assert translation["amino_acids"] == supplied["amino_acids"] == "MAQFVRNLVEKTPALVNG"
    assert translation["ends_with_stop_codon"] and translation["transcript_ids"] == supplied["donor_transcript_ids"]
    assert all(e["complete_5prime"] for e in translation["frame_evidence"])
    # Noncoding donor isoforms share the anchor, as in the supplied result.
    assert path["frame_status"] == "ambiguous" and translation["competing_noncoding_models"]
    peptides = {p["sequence"] for p in translation["candidate_peptides"]}
    assert {p["sequence"] for p in supplied["junction_peptides"]} <= peptides
    assert ({sam.split("\t")[0] for sam in result["original_records"].values()}
            == {r["sam"].split("\t")[0] for r in data["original_records"]})

    # Soft-clipped read-through adapter makes a second, technical path until
    # the explicit TruSeq profile trims it.
    clipped = run_corpus(bam, inputs, read_collector=ReadCollector(use_soft_clipped_bases=True))
    assert len(clipped["paths"]) == 2
    profile = ReadEndProfile("TruSeq", "1", (Adapter("TruSeq", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"),))
    trimmed = run_corpus(bam, inputs, read_collector=ReadCollector(
        use_soft_clipped_bases=True, read_end_profile=profile, trim_adapters=True))
    assert [p["sequence"] for p in trimmed["paths"]] == [path["sequence"]]


def test_sa_only_long_reads_do_not_become_a_placed_bcr_abl1_junction(tmp_path):
    data, bam, inputs = corpus_bam(tmp_path, "BCR--ABL1")
    supplied, = reconstruct_fusion(*fusion_from_dict(data))["translations"]
    result = run_corpus(bam, inputs)
    # The fixture retains primaries whose SA partners were not acquired.
    assert result["status"] == "no_candidate_paths"
    assert result["segment_path_notes"]["SA_declared_piece_not_linked"] == 2
    clipped = run_corpus(bam, inputs, read_collector=ReadCollector(use_soft_clipped_bases=True))
    assert clipped["status"] == "event_linked_candidates" and len(clipped["paths"]) == 2  # One per CCS read.
    assert {p["event_linkage"]["status"] for p in clipped["paths"]} == {"breakpoint_clip_partner_unplaced"}
    peptides = {p["sequence"] for path in clipped["paths"] for t in path["translations"]
                for p in t["candidate_peptides"] if "breakpoint_clip_partner_unplaced" in t["departure_relations"]}
    assert {p["sequence"] for p in supplied["junction_peptides"]} <= peptides
    assert any(supplied["amino_acids"] in t["amino_acids"] for p in clipped["paths"] for t in p["translations"])


@pytest.mark.parametrize("name,fragments,frame_statuses", [
    ("PARD3B--CDKN2B-AS1-CDKN2B-T1", 2, {"no_gene_anchor"}),
    ("FOXO3--STRADA-CCDC47-T1", 5, {"no_gene_anchor"}),
    ("TPST1--CRCP-T1", 61, {"reference_protein_only", "ambiguous"}),
])
def test_original_ont_split_reads_recover_the_junction_without_a_coding_frame(tmp_path, name, fragments,
                                                                                frame_statuses):
    data, bam, inputs = corpus_bam(tmp_path, name, "corpus")
    started = time.perf_counter()
    result = run_corpus(bam, inputs)
    assert time.perf_counter() - started < 30
    window = data["fusion"]
    event = [p for p in result["paths"] if p["event_linkage"]["status"] == "breakpoint_junction"]
    assert event and all(window["sequence"] in p["sequence"] for p in event[:1])
    # Unplaced junction bases are the supplied untemplated insertion plus
    # homology the supplied window assigned to one partner.
    junction = next(j for j in event[0]["junctions"] if j["relation"] == "breakpoint_junction")
    assert (len(junction["unplaced_bases"]) - sum(junction["breakpoint_assignment"])
            == window["junction_end"] - window["junction_start"])
    seeds = [s for s in result["seeds"] if s["relation"] == "breakpoint_junction"]
    assert sum(s["fragments"] for s in seeds) >= fragments
    # The supplied analysis found no coding donor frame; neither does this.
    assert {p["frame_status"] for p in event} <= frame_statuses
    assert not any("breakpoint_junction" in t["departure_relations"] for p in event for t in p["translations"])


LONG_READ = FUSIONS / "long-read"
LONG_READ_ENTRIES = {(e["event"], e["source"]): e for e in json.loads((LONG_READ / "manifest.json").read_text())}


def long_read_run(tmp_path, event, source):
    entry = LONG_READ_ENTRIES[event, source]
    raw = (LONG_READ / entry["file"]).read_bytes()
    assert sha256(raw).hexdigest() == entry["sha256"]
    lines = gzip.decompress(raw).decode().splitlines()
    header = pysam.AlignmentHeader.from_text("".join(line + "\n" for line in lines if line.startswith("@")))
    unsorted, path = tmp_path / "unsorted.bam", tmp_path / "long-read.bam"
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for line in lines:
            if not line.startswith("@"):
                out.write(pysam.AlignedSegment.fromstring(line, header))
    pysam.sort("-o", str(path), str(unsorted))
    pysam.index(str(path))
    data = json.loads(gzip.decompress((FUSIONS / entry["input"]).read_bytes()))
    fusion = data["fusion"]
    inputs = sv_rna_input_from_dict(dict(
        event_id=event, reference_name="GRCh38", sample_id="Sid-T1", donor=fusion["donor"],
        acceptor=fusion["acceptor"], references=data["references"], event_provenance=dict(fixture=entry["file"])))
    with pysam.AlignmentFile(str(path)) as bam:
        return reconstruct_sv_rna(bam, source=entry["url"], **inputs)


def adjacency_junction(result):
    return max((j for p in result["paths"] for j in p["junctions"] if j["relation"] == "breakpoint_junction"),
               key=lambda j: j["direct_fragments"])


def test_pacbio_places_homologous_junction_bases_past_the_catalogue_breakpoint(tmp_path):
    result = long_read_run(tmp_path, "TPST1--CRCP", "PacBio-T1")
    junction = adjacency_junction(result)
    # pbmm2 assigns all 8 homologous bases to CRCP; ONT reads leave them unplaced.
    assert junction["breakpoint_assignment"] == [3, -3] and junction["unplaced_bases"] == ""
    assert junction["direct_fragments"] == 20
    assert junction["direct_molecules"] is None
    assert junction["direct_cell_umi_support"]["status_counts"] == {"unresolved_xm": 20}
    assert result["status"] == "event_linked_candidates"


def test_noisy_ont_reads_count_as_direct_junction_support(tmp_path):
    result = long_read_run(tmp_path, "FOXO3--STRADA-CCDC47", "ONT-T1-tagged")
    junction = adjacency_junction(result)
    # v1.20.0 counted 0: no read matched the assembled consensus end to end.
    assert junction["direct_fragments"] == 12 and junction["direct_molecules"] is None
    assert junction["direct_cell_umi_support"]["observed_labels"] == 8
    assert junction["direct_cell_umi_support"]["unknown_library_segments"] == 12
    assert junction["direct_junction_sequences"][0] == ["GGA", 12]
    assert junction["direct_read_lineage"]["status_counts"] == {"unknown_producer": 12}


def test_regional_extensions_preserve_directly_witnessed_junction_orf(tmp_path):
    result = long_read_run(tmp_path, "FOXO3--STRADA-CCDC47", "ONT-T1-tagged")
    candidates = [c for p in result["paths"] for c in p["exploratory_orfs"]["candidates"]
                  if c["amino_acids"] == "MPLLYGYSVIEIYRRSNGTQP"]
    assert candidates
    assert {c["full_interval_support"]["fragments"] for c in candidates} == {9}
    # Alternative paths must not inflate support for the same original RNA.
    witnesses = {w["observation"] for c in candidates for w in c["full_interval_support"]["witnesses"]}
    assert len({tuple(result["observations"][w]["identity"][:2]) for w in witnesses}) == 9
    assert adjacency_junction(result)["direct_fragments"] == 12
    assert all(not c["initiation_observed"] and not c["translation_observed"] for c in candidates)


def test_long_reads_show_the_atp5mg_kmt2a_join_is_read_through_ambiguous(tmp_path):
    result = long_read_run(tmp_path, "ATP5MG--KMT2A", "PacBio-T1")
    junction, = [j for p in result["paths"] for j in p["junctions"]
                 if j["relation"] == "splice_ambiguous_event_junction" and j["breakpoint_assignment"] == [0, 0]][:1]
    assert junction["forward_splice_geometry"] and junction["direct_fragments"] == 5
    assert junction["direct_molecules"] is None
    assert junction["direct_cell_umi_support"]["status_counts"] == {"unresolved_xm": 5}
    assert result["status"] == "splice_ambiguous_candidates"
