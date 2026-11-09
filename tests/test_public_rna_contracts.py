"""Independent public consumers and shared-engine dependency boundaries."""

import ast
from pathlib import Path

import pysam
import pytest

from isovar.dna import reverse_complement_dna
from isovar.read_collector import ReadCollector
from isovar.rna_library import original_library_context
from isovar.rna_path_geometry import PathPositions, position_runs
from isovar.rna_reconstruction import discover_rna_paths, support_rna_paths

RNA = "ATGGCTACGTTGCAAGTCTAA"


def original_bam(tmp_path, strand, reads):
    """Write original records from independently specified RNA intervals."""
    path = tmp_path / "original.bam"
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "1", "LN": 1000}],
        "RG": [{"ID": "run", "SM": "sample", "LB": "library"}],
    })
    records = []
    for name, start, end, mate in reads:
        record = pysam.AlignedSegment(header)
        record.query_name = name
        sequence = RNA[start:end]
        record.query_sequence = sequence if strand == "+" else reverse_complement_dna(sequence)
        record.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
        record.flag = (16 if strand == "-" else 0) | (1 | mate if mate else 0)
        record.reference_id = 0
        record.reference_start = 10 + start if strand == "+" else 10 + len(RNA) - end
        record.mapping_quality = 60
        record.cigarstring = str(len(sequence)) + "M"
        record.set_tag("RG", "run")
        records.append(record)
    for record in records:
        if record.is_paired:
            partner = next(r for r in records if r.query_name == record.query_name and r is not record)
            record.next_reference_id = partner.reference_id
            record.next_reference_start = partner.reference_start
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for record in sorted(records, key=lambda r: (r.reference_start, r.query_name, r.flag)):
            out.write(record)
    pysam.index(str(path))
    return path


def candidate(strand, identifier="observed", sequence=RNA):
    positions = range(10, 10 + len(RNA))
    if strand == "-":
        positions = reversed(positions)
    return dict(id=identifier, sequence=sequence,
                placements=PathPositions(("1", p, strand) for p in positions),
                coding=[dict(status="complete_CDS", query_start=0, query_end=len(RNA))])


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("label", ["isoform", "sequence", "path"])
def test_public_discovery_budget_does_not_cap_original_support(tmp_path, strand, label):
    bam_path = original_bam(tmp_path, strand, [("a", 0, 21, 0), ("b", 0, 12, 0), ("c", 9, 21, 0)])
    collector = ReadCollector.for_rna_paths()
    regions = {"1": [(10, 31, strand)]}
    row = candidate(strand)
    with pysam.AlignmentFile(bam_path) as bam, original_library_context(
            bam, collector, tmp_path / "scratch", library_index=None,
            reference_name="test", read_filter_identity=None) as library:
        paths, discovery = discover_rna_paths(library, regions, collector, 1, 8, 3)
        assert paths == [(RNA, tuple(row["placements"]))]
        assert discovery["selected_segments"] == 1
        assert discovery["eligible_observed_segments"] == 3
        assert discovery["truncated"]
        support = support_rna_paths(library, [row], regions, collector, bam.header.to_dict(),
                                    "sample", "original-library", label=label)
    assert support["scope"] == "complete_declared_input"
    assert support["target_fragments"] == 3
    assert support["unassigned_target_fragments"] == 0
    assert row["rna_support"]["reads"] == row["rna_support"]["fragments"] == 3
    assert row["compatible_fragments"] == row["unique_fragments"] == 3
    assert row["coverage"] == [2] * 9 + [3] * 3 + [2] * 9
    assert len(row["whole_path_segments"]) == len(row["coding"][0]["whole_cds_segments"]) == 1
    assert row["coding"][0]["all_bases_observed"]
    assert position_runs(row["placements"]) == [(0, 21, "1", 10, 31, strand)]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_public_support_keeps_unsequenced_mate_gap_uncovered(tmp_path, strand):
    bam_path = original_bam(tmp_path, strand, [("paired", 0, 9, 64), ("paired", 12, 21, 128)])
    collector = ReadCollector.for_rna_paths()
    row = candidate(strand)
    with pysam.AlignmentFile(bam_path) as bam, original_library_context(
            bam, collector, tmp_path / "scratch", library_index=None,
            reference_name="test", read_filter_identity=None) as library:
        support = support_rna_paths(library, [row], {"1": [(10, 31, strand)]}, collector,
                                    bam.header.to_dict(), "sample", "original-library")
    assert support["target_fragments"] == 1
    assert row["rna_support"]["reads"] == 2
    assert row["rna_support"]["fragments"] == row["compatible_fragments"] == 1
    assert row["coverage"] == [1] * 9 + [0] * 3 + [1] * 9
    assert row["whole_path_segments"] == row["coding"][0]["whole_cds_segments"] == []
    assert not row["coding"][0]["all_bases_observed"]
    assert row["sequence_status"] == "partially_covered_reference"


def test_public_support_retains_ambiguous_isoform_membership(tmp_path):
    bam_path = original_bam(tmp_path, "+", [("shared", 0, 9, 0)])
    collector = ReadCollector.for_rna_paths()
    rows = [candidate("+", "first"), candidate("+", "second", RNA[:12] + "G" + RNA[13:])]
    with pysam.AlignmentFile(bam_path) as bam, original_library_context(
            bam, collector, tmp_path / "scratch", library_index=None,
            reference_name="test", read_filter_identity=None) as library:
        support = support_rna_paths(library, rows, {"1": [(10, 31, "+")]}, collector,
                                    bam.header.to_dict(), "sample", "original-library")
    assert support["target_fragments"] == 1
    assert support["equivalence_classes"] == [dict(isoform_ids=["first", "second"],
                                                    outside_targets=False, fragments=1)]
    assert [row["compatible_fragments"] for row in rows] == [1, 1]
    assert [row["unique_fragments"] for row in rows] == [0, 0]


def test_package_consumers_use_public_cross_module_contracts():
    package = Path(__file__).parents[1] / "isovar"
    violations = []
    paths = list(package.rglob("*.py")) + list((package.parent / "examples").rglob("*.py"))
    for path in paths:
        tree = ast.parse(path.read_text())
        aliases = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and ((node.level and path.is_relative_to(package)) or (node.module or "").startswith("isovar")):
                names = [part for part in (node.module or "").split(".") if part]
                names.extend(alias.name for alias in node.names)
                if any(name.startswith("_") and not name.startswith("__") for name in names):
                    violations.append((path.name, node.lineno, ast.unparse(node)))
                if node.module in (None, "isovar"):
                    aliases.update(alias.asname or alias.name for alias in node.names)
            elif isinstance(node, ast.Import):
                aliases.update(alias.asname or alias.name for alias in node.names
                               if alias.name.startswith("isovar."))
        for node in ast.walk(tree):
            if (path.is_relative_to(package) and isinstance(node, ast.Attribute) and ast.unparse(node.value) in aliases
                    and node.attr.startswith("_") and not node.attr.startswith("__")):
                violations.append((path.name, node.lineno, ast.unparse(node)))
    assert not violations, violations


def test_shared_engine_does_not_depend_on_antigen_adapters():
    package = Path(__file__).parents[1] / "isovar"
    for name in ("rna_path_geometry", "rna_observation", "rna_path_assembly", "rna_reconstruction"):
        for node in ast.walk(ast.parse((package / (name + ".py")).read_text())):
            if isinstance(node, ast.ImportFrom) and node.level:
                assert not (node.module or "").startswith(("sv_", "viral_", "antigen_")), (name, ast.unparse(node))
