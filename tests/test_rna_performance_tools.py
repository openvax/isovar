"""Performance evidence must fail closed on changed inputs or incomplete runs."""

from copy import deepcopy
import gzip
import json
import logging
from pathlib import Path
import shutil
from types import SimpleNamespace

import pytest

from isovar.allele_read import AlleleRead
from isovar.variant_sequence import VariantSequence
from tests.data.osteosarc.expansion import benchmark, collection_benchmark, collection_report, performance_report
from tests.data.osteosarc.expansion.inventory import write_json


def test_rna_fingerprint_preserves_identity_multiplicity_coverage_and_order():
    read = AlleleRead("AAA", "C", "GGG", "same", source_read_count=2)
    sequence = VariantSequence("AA", "C", "GG", [read])
    expected = benchmark.rna_fingerprint([sequence])
    assert benchmark.rna_fingerprint([deepcopy(sequence)]) == expected
    for name, value in [("prefix", "TAAA"), ("suffix", "GGGT"), ("allele", "T"),
                        ("name", "different"), ("source_read_count", 1)]:
        values = {k: getattr(read, k) for k in ("prefix", "allele", "suffix", "name", "source_read_count")}
        values[name] = value
        altered = VariantSequence("AA", "C", "GG", [AlleleRead(**values)])
        assert benchmark.rna_fingerprint([altered]) != expected
    longer = VariantSequence("AAA", "C", "GGG", [read])
    assert benchmark.rna_fingerprint([longer, sequence]) != benchmark.rna_fingerprint([sequence, longer])


@pytest.mark.parametrize("timeout", [0, -1, float("inf"), float("nan"), 1801])
def test_invalid_benchmark_budget_fails_before_creating_output(timeout):
    with pytest.raises(ValueError, match="finite timeout"):
        benchmark.run(SimpleNamespace(timeout=timeout))


def test_measurement_phases_are_separate_and_keep_failures(tmp_path):
    measurements = {}
    with benchmark.measure("production", tmp_path, measurements, 1, False):
        pass
    with pytest.raises(ValueError, match="test error"):
        with benchmark.measure("validation", tmp_path, measurements, 1, False):
            raise ValueError("test error")
    assert json.loads((tmp_path / "production.json").read_text())["status"] == "ok"
    assert json.loads((tmp_path / "validation.json").read_text())["status"] == "interrupted"
    assert set(measurements) == {"production", "validation"}


@pytest.fixture
def benchmark_inputs(tmp_path):
    """Stage original offline fixtures in the benchmark's acquisition layout."""
    corpus = Path(__file__).parent / "data/osteosarc/expansion/corpus"
    cases = json.loads((corpus / "manifest.json").read_text())["cases"]
    destination = tmp_path / "inputs"
    shutil.copytree(corpus / "references/GRCh38", destination / "reference-GRCh38")

    def prepare(variant_id, mode="defaults"):
        case = next(c for c in cases if c["variant"]["variant_id"] == variant_id and c["reference"] == "GRCh38")
        source = destination / "alignments" / case["source_id"]
        source.mkdir(parents=True)
        for suffix in ("", ".bai"):
            name = case["bam"] + suffix
            assert benchmark.digest(corpus / name) == case["files"][name]
            shutil.copyfile(corpus / name, source / ("regions-GRCh38.bam" + suffix))
        write_json(source / "regions-GRCh38.json", dict(bam_sha256=case["files"][case["bam"]]))
        write_json(destination / "inventory-GRCh38-validated.json", dict(variants=[case["variant"]]))
        return case, SimpleNamespace(destination=destination, source_id=case["source_id"],
                                     variant_id=variant_id, mode=mode, output=tmp_path / "run",
                                     timeout=60, profile=False)

    disabled = logging.root.manager.disable
    try:
        yield prepare
    finally:
        logging.disable(disabled)


@pytest.mark.parametrize("variant_id", [
    "ABCF2-chr7-151218156",  # Alternate evidence, but no RNA candidate.
    "ADGRF5-chr6-46856758",  # No alternate reads; ranking is never entered.
    "DYNC1H1-chr14-102030200",  # Positive control with multiple ranked proteins.
])
@pytest.mark.parametrize("mode", ["defaults", "primary_only"])
def test_benchmark_compares_real_empty_and_nonempty_results(tmp_path, benchmark_inputs, variant_id, mode):
    case, args = benchmark_inputs(variant_id, mode)
    entries = []
    for label in ("before", "after"):
        args.output = tmp_path / label
        benchmark.run(args)
        result = json.loads((args.output / "result.json").read_text())
        with gzip.open(args.output / "proteins.json.gz", "rt") as handle:
            proteins = json.load(handle)
        assert result["status"] == result["validation_status"] == "ok"
        assert result["counts"] == case[mode]["counts"]
        assert result["public_protein_count"] == len(case[mode]["proteins"])
        assert result["checked_protein_count"] == len(proteins)
        if case[mode]["proteins"]:
            assert result["checked_protein_count"] > result["public_protein_count"]
            assert all(p["validation_status"] == "ok" for p in proteins)
        else:
            assert result["checked_protein_count"] == 0
            assert result["rna_sha256"] == benchmark.rna_fingerprint([])
            assert result["trace"].get("rna_candidate_count", 0) == 0
        entries.append(f"{label}={args.output}")
    comparison = performance_report.collect_benchmarks(entries, ["before=after"])
    assert comparison["comparisons"] == [dict(before="before", after="after", exact_rna_and_ranked_proteins=True)]


@pytest.mark.parametrize("checked_count", [0, 1])
@pytest.mark.parametrize("error_type", [benchmark.AuditTimeout, ValueError])
def test_benchmark_retains_checked_count_on_validation_interruption(
        benchmark_inputs, monkeypatch, checked_count, error_type):
    _, args = benchmark_inputs("DYNC1H1-chr14-102030200")
    original_check = benchmark.protein_check
    completed = 0

    def interrupted_check(*args, **kwargs):
        nonlocal completed
        if completed == checked_count:
            raise error_type("injected validation interruption")
        result = original_check(*args, **kwargs)
        completed += 1
        return result

    monkeypatch.setattr(benchmark, "protein_check", interrupted_check)
    if error_type is ValueError:
        with pytest.raises(ValueError, match="injected validation interruption"):
            benchmark.run(args)
    else:
        benchmark.run(args)
    result = json.loads((args.output / "result.json").read_text())
    assert result["status"] == ("error" if error_type is ValueError else "resource_limit")
    assert result["trace"]["stage"] == "independent_protein_validation"
    assert result["checked_protein_count"] == completed == checked_count
    assert "validation_status" not in result
    assert not (args.output / "proteins.json.gz").exists()
    assert json.loads((args.output / "production.json").read_text())["status"] == "ok"
    assert json.loads((args.output / "validation.json").read_text())["status"] == "interrupted"
    with pytest.raises(ValueError, match="Incomplete comparison"):
        performance_report.collect_benchmarks([f"interrupted={args.output}"], ["interrupted=interrupted"])


def benchmark_pair(tmp_path):
    identity = dict(input_bam_sha256="bam", input_index_sha256="index", variant_id="variant",
                    mode="defaults", reference_manifest_sha256="reference")
    result = dict(status="ok", validation_status="ok", rna_sha256="rna", proteins_sha256="proteins",
                  counts={"alt": 2}, public_protein_count=1, checked_protein_count=20)
    directories = [tmp_path / label for label in ("before", "after")]
    for directory in directories:
        directory.mkdir()
        write_json(directory / "identity.json", identity)
        write_json(directory / "result.json", result)
        write_json(directory / "production.json", dict(wall_seconds=1))
    entries = [f"{p.name}={p}" for p in directories]
    return entries, directories


@pytest.mark.parametrize("variant_id", ["ABCF2-chr7-151218156", "DYNC1H1-chr14-102030200"])
@pytest.mark.parametrize("mode", ["defaults", "primary_only"])
@pytest.mark.parametrize("allocations", [False, True])
def test_collection_benchmark_preserves_original_evidence(benchmark_inputs, variant_id, mode, allocations):
    case, args = benchmark_inputs(variant_id, mode)
    args.allocations = allocations
    collection_benchmark.run(args)
    result = json.loads((args.output / "result.json").read_text())
    identity = json.loads((args.output / "identity.json").read_text())
    assert identity["allocation_instrumented"] is allocations
    assert result["status"] == "ok"
    assert result["counts"] == case[mode]["counts"]
    assert result["locus_read_count"] >= result["allele_read_count"]
    assert result["allele_read_count"] == sum(result["counts"]["read_objects"].values())
    for key in ("locus_reads_sha256", "allele_reads_sha256"):
        assert len(result[key]) == 64
    for key in ("locus_collection_wall_seconds", "locus_collection_cpu_seconds",
                "locus_collection_peak_rss_bytes", "post_collection_wall_seconds", "whole_process_peak_rss_bytes"):
        assert result[key] >= 0
    assert (args.output / "allocations.json").exists() is allocations
    if allocations:
        profile = json.loads((args.output / "allocations.json").read_text())
        assert profile["stage"] == "before_mate_merging"
        assert profile["locus_read_count"] >= result["locus_read_count"]
        assert profile["traced_peak_bytes"] >= profile["traced_current_bytes"] > 0
        assert profile["locations"]
    assert not collection_benchmark.tracemalloc.is_tracing()


def test_collection_benchmark_records_real_empty_bam_results(benchmark_inputs):
    import pysam

    _, args = benchmark_inputs("ABCF2-chr7-151218156")
    args.allocations = False
    source = args.destination / "alignments" / args.source_id
    bam = source / "regions-GRCh38.bam"
    with pysam.AlignmentFile(bam) as handle:
        header = handle.header
    # Only this test's staged copy becomes an empty, still-indexed BAM.
    with pysam.AlignmentFile(bam, "wb", header=header):
        pass
    pysam.index(str(bam))
    (source / "regions-GRCh38.json").write_text(json.dumps(dict(bam_sha256=benchmark.digest(bam))))
    collection_benchmark.run(args)
    result = json.loads((args.output / "result.json").read_text())
    assert result["status"] == "ok"
    assert result["locus_read_count"] == result["allele_read_count"] == 0
    assert result["locus_reads_sha256"] == result["allele_reads_sha256"] == collection_benchmark.read_fingerprint([])
    assert result["counts"]["reads"] == dict(ref=0, alt=0, other=0)


def test_collection_profiler_stops_tracing_before_snapshot_analysis(benchmark_inputs, monkeypatch):
    _, args = benchmark_inputs("ABCF2-chr7-151218156")
    args.allocations = True
    take_snapshot = collection_benchmark.take_snapshot
    analyzed = []

    def observed_snapshot():
        assert collection_benchmark.tracemalloc.is_tracing()
        # Only the benchmark's own name is patched, never the stdlib module.
        assert collection_benchmark.tracemalloc.take_snapshot is take_snapshot
        snapshot = take_snapshot()

        class ObservedSnapshot:
            def statistics(self, key):
                assert not collection_benchmark.tracemalloc.is_tracing()
                analyzed.append(key)
                return snapshot.statistics(key)

        return ObservedSnapshot()

    monkeypatch.setattr(collection_benchmark, "take_snapshot", observed_snapshot)
    collection_benchmark.run(args)
    assert analyzed == ["lineno"]


def test_collection_profiler_records_errors_and_restores_tracing(benchmark_inputs, monkeypatch):
    _, args = benchmark_inputs("ABCF2-chr7-151218156")
    args.allocations = True

    def fail(*args, **kwargs):
        assert collection_benchmark.tracemalloc.is_tracing()
        raise ValueError("injected collection failure")

    monkeypatch.setattr(collection_benchmark.ReadCollector, "get_locus_reads", fail)
    with pytest.raises(ValueError, match="injected collection failure"):
        collection_benchmark.run(args)
    assert not collection_benchmark.tracemalloc.is_tracing()
    result = json.loads((args.output / "result.json").read_text())
    assert result["status"] == "error"
    assert "injected collection failure" in result["error"]
    assert "counts" not in result  # An incomplete collection is not zero evidence.


def test_collection_error_after_counts_leaves_no_partial_evidence(benchmark_inputs, monkeypatch):
    _, args = benchmark_inputs("ABCF2-chr7-151218156")
    args.allocations = False
    peak_rss_bytes = collection_benchmark.peak_rss_bytes
    calls = []

    def fail_after_counts():
        calls.append(None)
        if len(calls) == 2:  # The final whole-process peak, after counts exist.
            raise KeyboardInterrupt
        return peak_rss_bytes()

    monkeypatch.setattr(collection_benchmark, "peak_rss_bytes", fail_after_counts)
    with pytest.raises(KeyboardInterrupt):
        collection_benchmark.run(args)
    assert len(calls) == 2
    result = json.loads((args.output / "result.json").read_text())
    assert result == dict(status="error", error="KeyboardInterrupt: ")


def test_collection_timeout_bounds_post_collection_work(benchmark_inputs, monkeypatch):
    import signal

    _, args = benchmark_inputs("ABCF2-chr7-151218156")
    args.allocations = False
    stages = []

    def bounded(name, function):
        def wrapped(*args, **kwargs):
            assert signal.getitimer(signal.ITIMER_REAL)[0] > 0, f"{name} ran without a time limit"
            stages.append(name)
            return function(*args, **kwargs)
        return wrapped

    for name in ("read_fingerprint", "allele_reads_from_locus_reads", "counts_from_evidence"):
        monkeypatch.setattr(collection_benchmark, name, bounded(name, getattr(collection_benchmark, name)))
    collection_benchmark.run(args)
    assert stages == ["read_fingerprint", "allele_reads_from_locus_reads", "read_fingerprint",
                      "counts_from_evidence"]
    assert signal.getitimer(signal.ITIMER_REAL) == (0.0, 0.0)


@pytest.mark.parametrize("collection_fails", [True, False])
def test_collection_record_write_failure_never_hides_the_outcome(
        benchmark_inputs, monkeypatch, capsys, collection_fails):
    _, args = benchmark_inputs("ABCF2-chr7-151218156")
    args.allocations = False
    write_json = collection_benchmark.write_json

    def disk_full(path, value):
        if path.name == "result.json":
            raise OSError(28, "No space left on device")
        return write_json(path, value)

    def fail(*args, **kwargs):
        raise ValueError("injected collection failure")

    monkeypatch.setattr(collection_benchmark, "write_json", disk_full)
    if collection_fails:
        monkeypatch.setattr(collection_benchmark.ReadCollector, "get_locus_reads", fail)
    with pytest.raises(ValueError if collection_fails else OSError):
        collection_benchmark.run(args)
    assert not (args.output / "result.json").exists()
    record = json.loads(capsys.readouterr().err.strip().splitlines()[-1])
    assert record["status"] == ("error" if collection_fails else "ok")
    assert ("counts" in record) is not collection_fails


@pytest.mark.parametrize("timeout", [0, -1, float("inf"), float("nan"), 1801])
def test_collection_benchmark_rejects_invalid_timeout(timeout):
    with pytest.raises(ValueError, match="finite timeout"):
        collection_benchmark.run(SimpleNamespace(timeout=timeout))


def test_collection_fingerprint_preserves_all_fields_and_order():
    from isovar.locus_read import LocusRead

    read = LocusRead("read", "ACG", [10000, None, 10001], [30, 0, 40], 10001, 10001, 1, 2)
    baseline = collection_benchmark.read_fingerprint([read])
    assert collection_benchmark.read_fingerprint([deepcopy(read)]) == baseline
    changes = dict(name="other", sequence="ATG", reference_positions=[10000, 10001, 10002],
                   quality_scores=[31, 0, 40], source_read_count=2, reference_base0_start_inclusive=10000,
                   reference_base0_end_exclusive=10002, read_base0_start_inclusive=0, read_base0_end_exclusive=3)
    assert set(changes) == set(read._fields)
    for field, value in changes.items():
        changed = deepcopy(read)
        setattr(changed, field, value)
        assert collection_benchmark.read_fingerprint([changed]) != baseline
        assert collection_benchmark.read_fingerprint([read, changed]) != collection_benchmark.read_fingerprint([changed, read])
    assert collection_benchmark.read_fingerprint([read, read]) != baseline


def collection_pair(tmp_path):
    identity = dict(source_bam_sha256="source", receipt_sha256="receipt", input_bam_sha256="bam",
                    input_index_sha256="index", inventory_sha256="inventory", benchmark_sha256="benchmark",
                    variant_id="variant", mode="defaults", merge_overlapping_fragments=True,
                    allocation_instrumented=False)
    result = dict(status="ok", locus_reads_sha256="locus", allele_reads_sha256="alleles",
                  locus_read_count=0, allele_read_count=0, counts={"alt": 0})
    entries, directories = [], []
    for label in ("before", "after"):
        directory = tmp_path / label
        directory.mkdir()
        write_json(directory / "identity.json", dict(identity, source_files={"isovar/read_collector.py": label}))
        write_json(directory / "result.json", result)
        entries.append(f"{label}={directory}")
        directories.append(directory)
    return entries, directories


def test_collection_comparison_accepts_complete_empty_results(tmp_path):
    entries, _ = collection_pair(tmp_path)
    result = collection_report.collect_runs(entries, ["before=after"])
    assert result["comparisons"] == [dict(before="before", after="after", exact_locus_and_allele_reads=True,
                                          allocation_instrumented=False)]


def test_collection_comparison_rejects_runs_of_identical_sources(tmp_path):
    entries, directories = collection_pair(tmp_path)
    for pair in ("after=after", "before=before"):
        with pytest.raises(ValueError, match="identical Isovar sources"):
            collection_report.collect_runs(entries, [pair])
    (directories[1] / "identity.json").write_text((directories[0] / "identity.json").read_text())
    with pytest.raises(ValueError, match="identical Isovar sources"):
        collection_report.collect_runs(entries, ["before=after"])


@pytest.mark.parametrize("filename,key,value", [
    ("identity", "source_bam_sha256", "different"), ("identity", "receipt_sha256", "different"),
    ("identity", "input_bam_sha256", "different"), ("identity", "input_index_sha256", "different"),
    ("identity", "inventory_sha256", "different"), ("identity", "benchmark_sha256", "different"),
    ("identity", "variant_id", "different"),
    ("identity", "mode", "primary_only"), ("identity", "merge_overlapping_fragments", False),
    ("identity", "allocation_instrumented", True), ("result", "status", "error"),
    ("result", "locus_reads_sha256", "different"), ("result", "allele_reads_sha256", "different"),
    ("result", "locus_read_count", 1), ("result", "allele_read_count", 1), ("result", "counts", {"alt": 1}),
])
def test_collection_comparison_rejects_changed_or_incomplete_results(tmp_path, filename, key, value):
    entries, directories = collection_pair(tmp_path)
    path = directories[1] / (filename + ".json")
    row = json.loads(path.read_text())
    row[key] = value
    path.write_text(json.dumps(row))
    with pytest.raises(ValueError):
        collection_report.collect_runs(entries, ["before=after"])


def test_benchmark_comparison_validates_complete_outputs(tmp_path):
    entries, _ = benchmark_pair(tmp_path)
    result = performance_report.collect_benchmarks(entries, ["before=after"])
    assert result["comparisons"] == [dict(before="before", after="after", exact_rna_and_ranked_proteins=True)]


@pytest.mark.parametrize("file,key,value", [
    ("identity", "mode", "primary_only"), ("identity", "input_bam_sha256", "different"),
    ("identity", "reference_manifest_sha256", "different"),
    ("result", "status", "resource_limit"), ("result", "validation_status", "error"),
    ("result", "rna_sha256", "different"), ("result", "proteins_sha256", "different"),
    ("result", "counts", {"alt": 1}), ("result", "checked_protein_count", 1),
])
def test_benchmark_comparison_rejects_changed_or_incomplete_evidence(tmp_path, file, key, value):
    entries, directories = benchmark_pair(tmp_path)
    path = directories[1] / (file + ".json")
    data = json.loads(path.read_text())
    data[key] = value
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        performance_report.collect_benchmarks(entries, ["before=after"])


@pytest.mark.parametrize("n_names", [0, 1, 7, 8, 9, 31, 100])
def test_support_bitsets_round_trip_without_modifying_input(n_names):
    names = [f"read-{i:03d}" for i in range(n_names)]
    proteins = [dict(supporting_read_names=names[::step], checks=[2, 1]) for step in (1, 2, 3)]
    row = dict(defaults=dict(proteins=proteins[:1], uncapped_ranked_proteins=proteins),
               primary_only=dict(proteins=proteins[1:]))
    untouched = deepcopy(row)
    packed = performance_report.pack_support_names(row)
    assert row == untouched
    for mode in ("defaults", "primary_only"):
        for field in row[mode]:
            for original, protein in zip(row[mode][field], packed[mode][field]):
                assert performance_report.supporting_read_names(packed, protein) == original["supporting_read_names"]
                assert protein["checks"] == [1, 2]


@pytest.mark.parametrize("changed,value", [("bits", ""), ("bits", "Ag=="), ("bits", "!invalid!"),
                                          ("names", ["different"]), ("encoding", "unknown")])
def test_support_bitsets_reject_corruption(changed, value):
    row = dict(defaults=dict(proteins=[dict(supporting_read_names=["original"], checks=[])]), primary_only={})
    packed = performance_report.pack_support_names(row)
    protein = packed["defaults"]["proteins"][0]
    if changed == "bits":
        protein["supporting_read_name_bits"] = value
    elif changed == "names":
        packed["supporting_read_name_table"] = value
    else:
        packed["support_encoding"] = value
    with pytest.raises(ValueError):
        performance_report.supporting_read_names(packed, protein)


def rerun_inputs(tmp_path, monkeypatch):
    monkeypatch.setattr(performance_report, "SOURCE_IDS", ["source"])
    audit, baseline = tmp_path / "audit", tmp_path / "baseline"
    for directory in (audit, baseline):
        (directory / "source").mkdir(parents=True)
    write_json(audit / "source/run.json", {})
    protein = dict(supporting_read_names=["original"], checks=[dict(status="ok")])
    mode = dict(status="ok", outcome="expected_top_protein", counts={"alt": 2},
                uncapped_validation_status="ok", proteins=[protein], uncapped_ranked_proteins=[protein])
    for variant_id in [performance_report.MT_VARIANT, *(f"nuclear-{i}" for i in range(43))]:
        current = dict(variant_id=variant_id, source_id="source", software="new",
                       defaults=deepcopy(mode), primary_only=deepcopy(mode), source_bam_sha256="bam")
        previous = deepcopy(current)
        previous["software"] = "old"
        if variant_id == performance_report.MT_VARIANT:
            for m in ("defaults", "primary_only"):
                previous[m].update(status="resource_limit", outcome="audit_timeout")
        write_json(audit / f"source/{variant_id}.json", current)
        write_json(baseline / f"source/{variant_id}.json", previous)
    return audit, baseline


def test_rerun_packaging_keeps_full_support_and_checks_all_controls(tmp_path, monkeypatch):
    result = performance_report.collect_reruns(*rerun_inputs(tmp_path, monkeypatch))
    assert result["unchanged_nuclear_rows"] == 43
    assert result["unchanged_nuclear_modes"] == 86
    assert result["recovered_mitochondrial_modes"] == 2
    assert len(result["checkpoints"]) == 44
    row, = result["rows"]
    assert performance_report.supporting_read_names(row, row["defaults"]["proteins"][0]) == ["original"]


@pytest.mark.parametrize("change", ["counts", "timeout", "validation", "empty", "input", "nuclear"])
def test_rerun_packaging_rejects_missing_or_changed_evidence(tmp_path, monkeypatch, change):
    audit, baseline = rerun_inputs(tmp_path, monkeypatch)
    variant_id = "nuclear-0" if change == "nuclear" else performance_report.MT_VARIANT
    path = audit / f"source/{variant_id}.json"
    row = json.loads(path.read_text())
    if change == "counts":
        row["defaults"]["counts"] = {"alt": 1}
    elif change == "timeout":
        row["defaults"]["status"] = "resource_limit"
    elif change == "validation":
        row["defaults"]["uncapped_validation_status"] = "error"
    elif change == "empty":
        row["defaults"]["proteins"] = []
    elif change == "input":
        row["source_bam_sha256"] = "different"
    else:
        row["defaults"]["proteins"][0]["supporting_read_names"] = ["different"]
    path.write_text(json.dumps(row))
    with pytest.raises(ValueError):
        performance_report.collect_reruns(audit, baseline)
