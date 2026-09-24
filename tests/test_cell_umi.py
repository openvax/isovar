"""#334: library-scoped labels and honest junction/ORF evidence denominators."""

from collections import defaultdict
import json

import pysam
import pytest

from isovar.cell_umi import CellUmiEvidence
from isovar.read_identity import segment_identity
from tests.test_sv_rna import HEADER, Scenario, aligned, record, write_bam


def tagged(name="read", group="a", flag=0, **tags):
    read = record(name, "1", 1000, "6M", "ATGAAA", flag=flag, group=group)
    for tag, value in tags.items():
        read.set_tag(tag, value)
    return read


def evidence(reads, read_groups=None, programs=(), sample="sample", source="source"):
    groups = defaultdict(list)
    for read in reads:
        groups[segment_identity(read)].append(read)
    if read_groups is None:
        read_groups = [dict(ID=rg, SM="sample", LB="library") for rg in ("a", "b")]
    return CellUmiEvidence(groups, dict(RG=read_groups, PG=programs), sample, source)


def support(reads, **kwargs):
    labels = evidence(reads, **kwargs)
    return labels.support(labels.groups)


@pytest.mark.parametrize("read_groups,count,complete", [
    ([dict(ID="a", SM="s", LB="A"), dict(ID="b", SM="s", LB="B")], 2, 2),
    ([dict(ID="a", SM="s", LB="A"), dict(ID="b", SM="s", LB="A")], 1, 1),
    ([dict(ID="a", SM="s1", LB="A"), dict(ID="b", SM="s2", LB="A")], 2, 2),
    ([dict(ID="a", LB="A"), dict(ID="b", LB="A")], 1, 1),
    ([dict(ID="a"), dict(ID="b")], 2, None),
    ([dict(ID="a", LB="A"), dict(ID="b")], 2, None),
    ([], 2, None),
])
def test_scope_uses_declared_library_and_sample_not_just_rg_or_barcode(read_groups, count, complete):
    reads = [tagged("same-name", group=rg, CB="cell-1", UB="umi") for rg in ("a", "b")]
    result = support(reads, read_groups=read_groups)
    assert result["observed_labels"] == count
    assert result["complete_label_count"] == complete
    assert result["unresolved_segments"] == 0
    assert result["independent_molecules"] is None
    assert len(result["segment_ids"]) == 2


def test_missing_rg_is_input_scoped_with_unknown_library_and_no_portable_identity_inference():
    read = tagged(CB="cell", UB="umi")
    read.set_tag("RG", None)
    labels = [evidence([read], sample=sample, source=source)
              for sample, source in (("T1", "raw"), ("T2", "raw"), ("T1", "reprocessed"))]
    keys = []
    for data in labels:
        result = data.support(data.groups)
        assert result["observed_labels"] == result["unknown_library_segments"] == 1
        assert result["complete_label_count"] is None
        row, = data.evidence()["segments"]
        assert row["scope"]["basis"] == "input" and row["scope"]["library"] is None
        keys.append(tuple(row["label"]))
    assert len(set(keys)) == 3  # Separate namespaces, not proof of independent RNA.


def test_exact_reported_labels_preserve_cell_suffixes_and_raw_correction_disagreements():
    reads = [tagged(str(i), CB=cb, UB=ub, CR="raw-cell", UR="raw-umi") for i, (cb, ub) in enumerate(
        [("cell-1", "umi"), ("cell-1", "umi"), ("cell-2", "umi"), ("cell-1", "other")])]
    originals = [r.to_string() for r in reads]
    result = support(reads + [reads[0]])
    assert result["observed_labels"] == result["complete_label_count"] == 3
    assert len(result["segment_ids"]) == 4
    assert [r.to_string() for r in reads] == originals


@pytest.mark.parametrize("tags,status", [
    ({}, "missing_cell_barcode"),
    (dict(CR="raw-cell", UB="umi"), "raw_cell_barcode_only"),
    (dict(CB="cell"), "missing_umi"),
    (dict(CB="cell", UR="raw"), "raw_umi_only"),
    (dict(CB="cell", RX="possibly-raw", MI="molecule"), "raw_umi_only"),
    (dict(CB="cell", XM="methylation-or-umi"), "unresolved_xm"),
    (dict(CB="cell", UB=""), "invalid_tags"),
    (dict(CB="cell", UB=42), "invalid_tags"),
    (dict(CB=" ", UB="umi"), "invalid_tags"),
])
def test_partial_and_invalid_labels_never_become_a_complete_denominator(tags, status):
    reads = [tagged("known", CB="cell", UB="umi"), tagged("unknown", **tags)]
    result = support(reads)
    assert result["observed_labels"] == result["unresolved_segments"] == 1
    assert result["complete_label_count"] is None and not result["all_segments_labeled"]
    assert result["status_counts"] == {"resolved_label": 1, status: 1}


@pytest.mark.parametrize("tags", [dict(CB="different", UB="umi"), dict(CB="cell", UB="different")])
@pytest.mark.parametrize("flag", [256, 2048])
def test_conflicting_placements_are_order_independent(tags, flag):
    primary, other = tagged(CB="cell", UB="umi"), tagged(flag=flag, **tags)
    for reads in ([primary, other], [other, primary]):
        result = support(reads)
        assert result["observed_labels"] == 0 and result["unresolved_segments"] == 1
        assert result["status_counts"] == {"conflicting_tags": 1}


def test_never_complete_a_label_by_combining_partial_tags_from_different_records():
    reads = [tagged(CB="cell"), tagged(flag=2048, UB="umi")]
    assert support(reads)["status_counts"] == {"partial_label": 1}
    reads.append(tagged(flag=256, CB="cell", UB="umi"))
    result = support(reads)
    assert result["complete_label_count"] == 1
    labels = evidence(reads)
    assert labels.segment(segment_identity(reads[0]))["records_without_complete_label"] == 2


def test_conflicting_mate_is_checked_even_when_only_one_mate_is_a_witness():
    reads = [tagged(flag=65, CB="cell", UB="umi"), tagged(flag=129, CB="cell", UB="other")]
    labels = evidence(reads)
    result = labels.support([segment_identity(reads[0])])
    assert result["status_counts"] == {"conflicting_template_labels": 1}
    assert result["observed_labels"] == 0
    assert len(labels.evidence()["segments"]) == 2


def test_missing_mate_label_is_not_invented_or_counted_as_a_second_molecule():
    reads = [tagged(flag=65, CB="cell", UB="umi"), tagged(flag=129)]
    result = support(reads)
    assert result["observed_labels"] == result["unresolved_segments"] == 1
    assert result["complete_label_count"] is None


def test_duplicate_header_rg_is_unresolved_instead_of_selecting_the_first_library():
    result = support([tagged(CB="cell", UB="umi")], read_groups=[dict(ID="a", LB="A"), dict(ID="a", LB="B")])
    assert result["status_counts"] == {"ambiguous_read_group": 1}
    assert result["observed_labels"] == 0 and result["complete_label_count"] is None


CORRECT = dict(ID="correct", PN="isoseq", CL="isoseq correct input.bam output.bam")
ALIGN = dict(ID="align", PN="pbmm2", PP="correct", CL="pbmm2 align input.bam output.bam")


@pytest.mark.parametrize("programs,pg,expected", [
    ([CORRECT, ALIGN], None, 1),
    ([dict(ID="isoseq correct", CL="correct in.bam out.bam")], None, 1),
    ([dict(ID="tag", PN="isoseq", CL="isoseq tag in.bam out.bam")], None, 0),
    ([dict(ID="dedup", PN="isoseq", CL="isoseq groupdedup in.bam out.bam")], None, 0),
    ([dict(ID="bismark", PN="bismark", CL="bismark reads.fastq")], None, 0),
    ([], None, 0),
    ([CORRECT, dict(ALIGN, PP="missing")], None, 0),
    ([CORRECT, dict(ALIGN, PP="align")], None, 0),
    ([CORRECT, dict(ID="align", PN="pbmm2")], None, 0),
    ([CORRECT, dict(ID="bismark", PN="bismark")], "correct", 1),
    ([CORRECT, dict(ID="bismark", PN="bismark")], "bismark", 0),
    ([CORRECT, dict(ID="bismark", PN="bismark")], None, 0),
    ([CORRECT], "missing", 0),
    ([CORRECT, dict(CORRECT, CL="isoseq tag in out")], "correct", 0),
    ([CORRECT, dict(ID="retag", PN="isoseq", CL="isoseq tag in out", PP="correct")], None, 0),
])
def test_xm_requires_attributable_isoseq_correction_not_platform_or_tag_presence(programs, pg, expected):
    read = tagged(CB="cell", XM="umi")
    if pg is not None:
        read.set_tag("PG", pg)
    result = support([read], read_groups=[dict(ID="a", PL="PACBIO", LB="library")], programs=programs)
    assert result["observed_labels"] == expected
    assert result["status_counts"] == {"resolved_label" if expected else "unresolved_xm": 1}


def test_rg_program_pointer_and_ub_xm_disagreement():
    read = tagged(CB="cell", XM="isoseq-umi", UB="other-umi")
    groups = [dict(ID="a", LB="library", PG="correct")]
    programs = [CORRECT, dict(ID="unrelated", PN="bismark")]
    assert support([read], read_groups=groups, programs=programs)["status_counts"] == {"conflicting_tags": 1}
    read.set_tag("UB", "isoseq-umi")
    assert support([read], read_groups=groups, programs=programs)["complete_label_count"] == 1
    # A producer's unrelated XM must not override an explicit UB label.
    read.set_tag("PG", "unrelated")
    read.set_tag("XM", "..zzHH")
    assert support([read], read_groups=groups, programs=programs)["complete_label_count"] == 1


def test_conflicting_xm_provenance_across_placements_cannot_use_just_the_primary():
    reads = [tagged(CB="cell", XM="umi", PG="correct"),
             tagged(flag=2048, CB="cell", XM="unrelated", PG="bismark")]
    programs = [CORRECT, dict(ID="bismark", PN="bismark")]
    for records in (reads, reads[::-1]):
        result = support(records, programs=programs)
        assert result["observed_labels"] == 0 and result["complete_label_count"] is None
        assert result["status_counts"] == {"unresolved_xm": 1}


@pytest.mark.parametrize("libraries,missing,count", [(('A', 'B'), False, 2), (('A', 'A'), False, 1),
                                                    (('A', 'B'), True, None), ((None, None), False, None)])
def test_junction_and_full_orf_use_the_same_scoped_evidence(tmp_path, libraries, missing, count):
    scenario = Scenario()
    reads = []
    header = HEADER.to_dict()
    header["RG"] = [dict(ID=rg, SM="sample", **({"LB": library} if library else {}))
                    for rg, library in zip(("a", "b"), libraries)]
    for rg in ("a", "b"):
        for read in aligned("full-" + rg, scenario.sequence, scenario.positions, group=rg):
            if not (missing and rg == "b"):
                read.set_tag("CB", "cell")
                read.set_tag("UB", "umi")
            reads.append(read)
    result = scenario.run(write_bam(tmp_path / "rna.bam", reads, header=header))
    path, = result["paths"]
    junction, = path["junctions"]
    assert junction["direct_fragments"] == junction["direct_segments"] == 2
    assert junction["direct_molecules"] == count
    orfs = path["exploratory_orfs"]["candidates"]
    assert orfs
    for orf in orfs:
        full = orf["full_interval_support"]
        assert full["segments"] == 2 and full["molecule_labels"] == count
        assert full["cell_umi_support"] == junction["direct_cell_umi_support"]
    assert json.loads(json.dumps(result))["cell_umi_evidence"] == result["cell_umi_evidence"]


def test_partial_junction_witness_is_not_promoted_to_full_orf_support(tmp_path):
    scenario = Scenario()
    reads = aligned("full", scenario.sequence, scenario.positions)
    reads += aligned("partial", scenario.sequence[230:280], scenario.positions[230:280])
    for read in reads:
        if read.query_name == "full":
            read.set_tag("CB", "cell")
            read.set_tag("UB", "umi")
    header = HEADER.to_dict()
    header["RG"] = [dict(ID="a", LB="library", SM="sample")]
    result = scenario.run(write_bam(tmp_path / "rna.bam", reads, header=header))
    path, = result["paths"]
    junction, = path["junctions"]
    assert junction["direct_segments"] == 2 and junction["direct_molecules"] is None
    assert junction["direct_cell_umi_support"]["unresolved_segments"] == 1
    orf = next(c for c in path["exploratory_orfs"]["candidates"] if c["query_interval"][0] == 20)
    full = orf["full_interval_support"]
    assert full["segments"] == full["molecule_labels"] == 1
    assert full["cell_umi_support"]["unresolved_segments"] == 0
