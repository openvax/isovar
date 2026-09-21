"""Candidate initiation is separate from annotation and full-span RNA support."""

from dataclasses import replace
from types import SimpleNamespace

import pysam

from isovar.sv_rna import RnaObservation, _Model
from isovar.sv_rna_orfs import exploratory_orfs, record_evidence
from tests.test_sv_rna import long_read_run


def observation(key, sequence, positions, **kwargs):
    return RnaObservation(key=key, identity=("library", key, 0), records=(key,), sequence=sequence,
                          positions=tuple(positions), breaks=((5, 6, "split"),), reverse=False,
                          query_interval=(10, 10 + len(sequence)), missing_qualities=False,
                          secondary=False, **kwargs)


def inputs(sequence="ATGAAAATGAAATAA"):
    positions = tuple(("1" if q < 6 else "2", q, "+") for q in range(len(sequence)))
    reads = {"full": observation("full", sequence, positions)}
    junction = dict(query_interval=[5, 6], left=positions[5], right=positions[6], annotated=False,
                    relation="breakpoint_junction", direct_observations=list(reads))
    return sequence, positions, reads, junction


def run(sequence, positions, reads, junction, models=(), minimum=1, limit=100, molecule=lambda identity: None):
    return exploratory_orfs(sequence, positions, [junction], reads, models, molecule, minimum, limit)


def test_only_orfs_crossing_the_junction_are_reported_with_observed_stop():
    seq, pos, reads, junction = inputs()
    result = run(seq, pos, reads, junction)
    candidate, = result["candidates"]
    assert candidate["amino_acids"] == "MKMK"
    assert candidate["ends_with_stop_codon"] and candidate["query_interval"] == [0, 15]
    assert not candidate["initiation_observed"] and not candidate["translation_observed"]
    assert not candidate["annotated_start"] and candidate["reference_comparisons"] == []
    assert candidate["start_context"]["minus_three"] is None
    assert candidate["full_interval_support"]["fragments"] == 1
    assert not run(seq, pos, reads, junction, minimum=5)["candidates"]


def test_competing_starts_and_partial_terminal_codon_are_preserved_and_limit_is_explicit():
    seq, pos, reads, junction = inputs("ATGATGAAAC")
    result = run(seq, pos, reads, junction)
    assert [c["amino_acids"] for c in result["candidates"]] == ["MMK", "MK"]
    assert all(not c["ends_with_stop_codon"] for c in result["candidates"])
    assert [c["query_interval"] for c in result["candidates"]] == [[0, 9], [3, 9]]
    limited = run(seq, pos, reads, junction, limit=1)
    assert len(limited["candidates"]) == 1 and limited["candidate_limit_reached"]


def test_partial_reads_conflicting_placements_and_different_unplaced_bases_are_not_full_witnesses():
    seq, pos, reads, junction = inputs()
    reads["left"] = observation("left", seq[:9], pos[:9])
    reads["right"] = observation("right", seq[3:], pos[3:])
    reads["elsewhere"] = observation("elsewhere", seq, [("3", q, "+") for q in range(len(seq))])
    reads["missing"] = observation("missing", seq[:7] + seq[8:], pos[:7] + pos[8:])
    changed = list(pos)
    changed[7] = None
    reads["unplaced_error"] = observation("unplaced_error", seq[:7] + "C" + seq[8:], changed)
    junction["direct_observations"] = list(reads)
    candidate, = run(seq, pos, reads, junction)["candidates"]
    support = candidate["full_interval_support"]
    assert support["fragments"] == 1
    assert [w["observation"] for w in support["witnesses"]] == ["full"]
    # A stop in the path is not itself evidence for a complete RNA molecule.
    del reads["full"]
    junction["direct_observations"] = list(reads)
    candidate, = run(seq, pos, reads, junction)["candidates"]
    assert candidate["ends_with_stop_codon"] and candidate["full_interval_support"]["fragments"] == 0


def test_reverse_observation_offsets_and_read_group_scoped_molecule_labels():
    seq, pos, reads, junction = inputs()
    for key, rg in (("same", "library"), ("other", "another-library")):
        reads[key] = replace(reads["full"], key=key, identity=(rg, key, 0), reverse=True,
                             missing_qualities=True, sequence="CCC" + seq, positions=(None,) * 3 + pos,
                             query_interval=(10, 28))
    junction["direct_observations"] = list(reads)
    candidate, = run(seq, pos, reads, junction, molecule=lambda identity: ("cell", "umi"))["candidates"]
    support = candidate["full_interval_support"]
    assert support["segments"] == support["fragments"] == 3
    assert support["molecule_labels"] == 2 and support["missing_quality_segments"] == 2
    assert all(w["original_query_interval"] == [10, 25] for w in support["witnesses"])


def test_reference_start_classification_compares_same_atg_without_claiming_initiation():
    seq, pos, reads, junction = inputs()
    ref = SimpleNamespace(contig="1", strand="+", exons=((0, 15),), sequence="ATGAAACCCTAATAA",
                          cds_start=0, cds_end=9, transcript_id="coding", annotation="test")
    models = [_Model(ref), _Model(SimpleNamespace(**dict(vars(ref), cds_start=9, cds_end=12, transcript_id="utr")))]
    candidate, = run(seq, pos, reads, junction, models=models)["candidates"]
    assert candidate["annotated_start"] and not candidate["initiation_observed"]
    first, second = candidate["reference_comparisons"]
    assert first["start_kind"] == "annotated_start" and second["start_kind"] == "five_prime_UTR"
    assert first["amino_acids"] == "MKP" and first["shared_prefix_amino_acids"] == 2
    assert first["differs_from_reference_orf"]


def test_record_evidence_preserves_zero_missing_qual_and_mapq_unavailable():
    read = pysam.AlignedSegment()
    read.query_sequence = "A" * 12
    read.cigarstring = "2S4M2I4M"
    read.mapping_quality = 255
    read.set_tag("NM", 0)
    read.set_tag("rc", 0)
    read.set_tag("rq", 0.99)
    read.set_tag("np", 12)
    read.set_tag("ec", 11.5)
    evidence = record_evidence(read)
    assert evidence["mapping_quality"] is None and not evidence["base_qualities_available"]
    assert evidence["aligned_query_bases"] == 10
    assert evidence["tags"]["NM"] == evidence["tags"]["rc"] == 0
    assert evidence["tags"]["mg"] is None and evidence["tags"]["ic"] is None
    assert evidence["tags"]["rq"] > 0.98 and evidence["tags"]["np"] == 12 and evidence["tags"]["ec"] == 11.5
    read.mapping_quality = 0
    read.query_qualities = [0] * 12
    assert record_evidence(read)["mapping_quality"] == 0
    assert record_evidence(read)["base_qualities_available"]


def test_original_pacbio_upstream_orf_has_full_span_support_and_retains_tags(tmp_path):
    result = long_read_run(tmp_path, "TPST1--CRCP", "PacBio-T1")
    candidates = [c for p in result["paths"] for c in p["exploratory_orfs"]["candidates"]
                  if c["amino_acids"].startswith("MPSRRRGGSSLIP")]
    assert {len(c["amino_acids"]): c["full_interval_support"]["fragments"] for c in candidates} == {30: 15, 61: 1, 72: 1}
    candidate, = [c for c in candidates if len(c["amino_acids"]) == 30]
    assert candidate["amino_acids"] == "MPSRRRGGSSLIPGRMPILRFSVTTRYFSY"
    assert candidate["ends_with_stop_codon"] and not candidate["annotated_start"]
    assert candidate["start_context"]["sequence"] == "GCCAGGATGCCGTCC"
    reference, = [r for r in candidate["reference_comparisons"] if r["transcript_id"] == "ENST00000304842"]
    assert reference["start_kind"] == "five_prime_UTR"
    assert reference["transcript_offset"] == 284 and reference["annotated_cds_start"] == 425
    assert reference["amino_acids"] == "MPSRRRGGSSLIPDVGYLSEVDCPWPEHFPKIILSKISV"
    assert reference["shared_prefix_amino_acids"] == 13
    support = candidate["full_interval_support"]
    assert support["molecule_labels"] == support["missing_quality_segments"] == 15
    record_ids = {rid for w in support["witnesses"] for rid in result["observations"][w["observation"]]["records"]}
    evidence = [result["record_evidence"][rid] for rid in record_ids]
    assert all(not r["base_qualities_available"] and r["mapping_quality"] == 60 for r in evidence)
    assert all(r["tags"]["mg"] > 97 and r["tags"]["rq"] is None and r["tags"]["np"] is None for r in evidence)
    assert all(r["tags"]["CB"] and r["tags"]["XM"] and r["tags"]["im"] for r in evidence)
    assert {r["tags"]["ic"] for r in evidence} == {1, 2}
    assert {r["tags"]["rc"] for r in evidence} == {0, 1}
    assert {r["tags"]["rm"] for r in evidence} == {1}
