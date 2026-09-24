"""Dorado relationships are signal ancestry, not alignment or molecule counts."""

from collections import defaultdict
from copy import deepcopy
from dataclasses import replace

import pysam
import pytest

from isovar.cell_umi import CellUmiEvidence
from isovar.read_lineage import ReadLineage
from isovar.sv_rna import segment_identity
from isovar.sv_rna_orfs import exploratory_orfs
from tests.test_sv_rna import HEADER, Scenario, aligned, record, spanning
from tests.test_sv_rna_orfs import inputs, observation


def header():
    result = HEADER.to_dict()
    result["RG"] = [dict(ID=rg, PL="ONT", LB="same-library", SM="sample") for rg in ("a", "b")]
    result["PG"] = [dict(ID="basecaller", PN="dorado", CL="dorado basecaller model reads/"),
                    dict(ID="aligner", PN="minimap2", PP="basecaller")]
    return result


def read(name, group="a", flag=0, **tags):
    result = record(name, "1", 100, "15M", "ATGAAAATGAAATAA", group=group, flag=flag)
    for tag, value in tags.items():
        result.set_tag(tag, value)
    return result


def lineage(records, metadata=None):
    groups = defaultdict(list)
    for r in records:
        groups[segment_identity(r)].append(r)
    return ReadLineage(groups, header() if metadata is None else metadata)


def test_parent_and_split_children_share_signal_but_not_segment_identity():
    records = [read("parent", dx=0), read("child1", pi="parent", sp=0, dx=0),
               read("child2", pi="parent", sp=500, dx=0), read("child1", flag=256, pi="parent", sp=0, dx=0),
               read("child1", group="b", pi="parent", sp=0, dx=0)]
    evidence = lineage(records)
    result = evidence.support(evidence.groups)
    assert len(result["segment_ids"]) == 4
    assert result["resolved_signal_groups"] == 2  # Identical LB does not merge RGs.
    assert result["unresolved_segments"] == 0
    assert result["status_counts"] == {"simplex_read": 1, "split_read": 3}
    assert evidence.segment(("a", "child2", 0))["signal_group"] == ["a", "parent"]


def test_duplex_status_does_not_invent_parent_consensus_links_from_names():
    evidence = lineage([read("template", dx=-1), read("complement", dx=-1),
                        read("template;complement", dx=1), read("unpaired", dx=0)])
    result = evidence.support(evidence.groups)
    assert result["resolved_signal_groups"] == 1 and result["unresolved_segments"] == 3
    assert not result["all_segments_resolved"]
    assert result["status_counts"] == {"duplex_parent": 2, "duplex_consensus": 1, "simplex_read": 1}
    assert all(row["signal_group"] is None for row in evidence.evidence() if row["duplex_status"] != 0)


@pytest.mark.parametrize("change", ["no_pg", "no_rg", "other_platform", "unknown_chain", "broken_chain",
                                     "cycle", "duplicate_pg", "duplicate_rg", "wrong_pointer",
                                     "aligner_only", "missing_command", "malformed_command"])
def test_tags_and_unrelated_header_programs_do_not_identify_a_producer(change):
    metadata = header()
    r = read("child", pi="parent", dx=0)
    if change == "no_pg":
        metadata.pop("PG")
    elif change == "no_rg":
        metadata.pop("RG")
    elif change == "other_platform":
        metadata["RG"][0]["PL"] = "PACBIO"
    elif change == "unknown_chain":
        metadata["PG"].append(dict(ID="other", PN="unknown"))
    elif change == "broken_chain":
        metadata["PG"][0]["PP"] = "missing"
    elif change == "cycle":
        metadata["PG"][0]["PP"] = "aligner"
    elif change == "duplicate_pg":
        metadata["PG"] += deepcopy(metadata["PG"])
    elif change == "duplicate_rg":
        metadata["RG"] += deepcopy(metadata["RG"])
    elif change in ("aligner_only", "missing_command", "malformed_command"):
        metadata["PG"][0]["CL"] = {"aligner_only": "dorado aligner reference reads",
                                   "missing_command": "", "malformed_command": 'dorado basecaller "'}[change]
    else:
        r.set_tag("PG", "missing")
    row = lineage([r], metadata).segment(("a", "child", 0))
    assert row["status"] == "unknown_producer" and row["signal_group"] is None


@pytest.mark.parametrize("pointer", ["record", "group"])
def test_explicit_program_pointer_disambiguates_a_merged_header(pointer):
    metadata = header()
    metadata["PG"].append(dict(ID="other", PN="unknown"))
    r = read("child", pi="parent")
    if pointer == "record":
        r.set_tag("PG", "aligner")
    else:
        metadata["RG"][0]["PG"] = "aligner"
    row = lineage([r], metadata).segment(("a", "child", 0))
    assert row["signal_group"] == ["a", "parent"]
    assert row["duplex_status"] is None  # Simplex basecalling need not supply dx.


@pytest.mark.parametrize("tags", [dict(pi=""), dict(pi="child"), dict(pi=3), dict(dx=2), dict(dx=1.0),
                                  dict(sp=0), dict(pi="parent", sp=-1)])
def test_invalid_lineage_tags_remain_unresolved(tags):
    row = lineage([read("child", **tags)]).segment(("a", "child", 0))
    assert row["status"] == "invalid_tags" and row["signal_group"] is None


@pytest.mark.parametrize("other", [dict(pi="different", dx=0), dict(dx=0), dict(pi="parent", dx=1),
                                   dict(pi="parent", dx=0.0), dict(pi="parent", dx=0, sp=5)])
def test_conflicting_placements_are_order_independent(other):
    records = [read("child", pi="parent", dx=0), read("child", flag=2048, **other)]
    for ordered in (records, records[::-1]):
        row = lineage(ordered).segment(("a", "child", 0))
        assert row["status"] == "conflicting_tags" and row["signal_group"] is None


def test_visible_parent_chains_resolve_without_promoting_parent_to_support():
    evidence = lineage([read("child", pi="middle"), read("middle", pi="root"), read("root")])
    support = evidence.support([("a", "child", 0)])
    assert support["segment_ids"] == [["a", "child", 0]]
    assert support["resolved_signal_groups"] == 1
    assert len(evidence.evidence()) == 3
    assert evidence.segment(("a", "child", 0))["signal_group"] == ["a", "root"]


@pytest.mark.parametrize("parent_tags", [dict(dx=1), dict(pi="child")])
def test_unresolved_or_cyclic_parent_cannot_resolve_a_child(parent_tags):
    evidence = lineage([read("child", pi="parent"), read("parent", **parent_tags)])
    assert evidence.support(evidence.groups)["unresolved_segments"] == 2
    assert evidence.segment(("a", "child", 0))["status"] == "unresolved_parent"


def test_pair_bits_do_not_turn_ont_records_into_independent_signals():
    evidence = lineage([read("paired", flag=65)])
    assert evidence.segment(("a", "paired", 64))["status"] == "unsupported_segment"


def test_lineage_metadata_is_cached_without_accessing_sequences_or_bulk_tags():
    class MetadataOnly:
        def __init__(self, r):
            self.r, self.calls = r, 0

        def __getattr__(self, name):
            assert name not in ("query_sequence", "query_qualities", "cigartuples", "get_tags")
            return getattr(self.r, name)

        def get_tag(self, tag):
            self.calls += 1
            return self.r.get_tag(tag)

    r = MetadataOnly(read("child", pi="parent", dx=0, sp=0))
    evidence = ReadLineage({("a", "child", 0): [r]}, header())
    first = evidence.support(evidence.groups)
    calls = r.calls
    assert calls > 0
    assert evidence.support(evidence.groups) == first and r.calls == calls


def test_full_orf_lineage_uses_only_complete_compatible_witnesses():
    sequence, positions, observations, junction = inputs()
    observations["full"] = replace(observations["full"], identity=("a", "full", 0))
    for name, lo, hi in (("full2", 0, 15), ("left", 0, 9), ("right", 3, 15)):
        observations[name] = replace(observation(name, sequence[lo:hi], positions[lo:hi]), identity=("a", name, 0))
    observations["elsewhere"] = replace(observations["full"], key="elsewhere",
                                        positions=tuple(("3", q, "+") for q in range(15)))
    evidence = lineage([read(name, pi="parent" if name.startswith("full") else "partial-parent")
                        for name in ("full", "full2", "left", "right")])
    junction["direct_observations"] = list(observations)
    labels = CellUmiEvidence(evidence.groups, {}, "sample", "source")
    result = exploratory_orfs(sequence, positions, [junction], observations, [], labels.support,
                              1, 100, lineage=evidence.support)
    candidate, = result["candidates"]
    support = candidate["full_interval_support"]
    assert support["segments"] == 2 and len(support["witnesses"]) == 2
    assert support["read_lineage"]["resolved_signal_groups"] == 1
    assert support["read_lineage"]["segment_ids"] == [["a", "full", 0], ["a", "full2", 0]]
    # The two partial siblings still cannot establish the complete ORF together.
    junction["direct_observations"] = ["left", "right"]
    candidate, = exploratory_orfs(sequence, positions, [junction], observations, [], labels.support,
                                  1, 100, lineage=evidence.support)["candidates"]
    assert candidate["full_interval_support"]["read_lineage"]["resolved_signal_groups"] == 0


def test_sv_junctions_export_signal_ancestry_without_changing_reconstruction(tmp_path):
    scenario = Scenario()
    records = []
    for name in ("child1", "child2"):
        for r in aligned(name, scenario.sequence, scenario.positions):
            r.set_tag("pi", "parent")
            r.set_tag("dx", 0)
            records.append(r)
    unsorted, bam = tmp_path / "unsorted.bam", tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(unsorted), "wb", header=header()) as out:
        for r in records:
            out.write(r)
    pysam.sort("-o", str(bam), str(unsorted))
    pysam.index(str(bam))
    result = scenario.run(bam)
    assert result["paths"][0]["sequence"] == scenario.sequence
    junction, = spanning(result, "event_compatible_junction")
    assert junction["direct_segments"] == junction["direct_fragments"] == 2
    assert junction["direct_read_lineage"]["resolved_signal_groups"] == 1
    assert len(result["read_lineage"]["segments"]) == 2
    assert all(row["signal_group"] == ["a", "parent"] for row in result["read_lineage"]["segments"])


def test_unbuilt_direct_junction_reads_still_contribute_lineage(tmp_path):
    scenario = Scenario()
    sequence = scenario.donor_ref.sequence[40:100] + "TTT" + scenario.donor_ref.sequence[250:350]
    unsorted, bam = tmp_path / "unsorted.bam", tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(unsorted), "wb", header=header()) as out:
        for i in range(8):
            r = record("child%d" % i, "1", 1040, "60M3I1900N100M", sequence)
            r.set_tag("pi", "parent%d" % (i // 2))
            out.write(r)
    pysam.sort("-o", str(bam), str(unsorted))
    pysam.index(str(bam))
    result = scenario.run(bam, references=[scenario.donor_ref], max_extension_segments=2)
    junction, = spanning(result, "regional_novel_junction")
    assert junction["direct_segments"] == 8
    assert result["observation_counts"]["built_segments"] == 2
    assert junction["direct_read_lineage"]["resolved_signal_groups"] == 4
    assert len(result["read_lineage"]["segments"]) == 8


def test_header_lines_without_ids_are_ignored_not_fatal():
    metadata = header()
    metadata["RG"].append(dict(PL="ONT", SM="sample"))
    metadata["PG"].append(dict(PN="samtools"))
    evidence = lineage([read("x")], metadata)
    assert set(evidence.read_groups) == {"a", "b"}
    assert set(evidence.programs) == {"basecaller", "aligner"}
