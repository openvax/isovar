from copy import deepcopy
from pathlib import Path

import pytest

from examples.recover_sid_sv_pilot import clone_cache, coverage, failed_pairs, prepare
from examples.sid_sv_audit.inventory import digest, identity, read_json, write_json
from examples.sid_sv_audit.reconstruct import implementation_identity, run_batch
from tests.test_sid_sv_audit import empty_audit as empty_audit


def test_successor_producer_requires_unchanged_inputs_and_preserved_failure(tmp_path):
    from examples.check_sid_sv_recovery import view_producer
    from examples.sid_sv_audit.reconstruct import result_path

    request = dict(source_id='source', geometry_id='geometry', orientation='forward',
                   inventory_sha256='inventory', references_sha256='references', acquisition_sha256='reads',
                   parameters=dict(dense_support=True, max_records=50000),
                   engine_id='old-engine', isovar_version='1.0.0', max_seconds=14400)
    old = dict(request=request, status='reconstruction_timeout')
    successor = dict(request=dict(request, engine_id='streaming-engine', isovar_version='2.0.0', max_seconds=3600),
                     status='no_candidate_paths')
    run = {field: successor['request'][field] for field in ('engine_id', 'isovar_version', 'parameters', 'max_seconds')}
    old_path, run_path = tmp_path / 'history/failed.json', tmp_path / 'controls/producer.json'
    write_json(old_path, old)
    write_json(run_path, run)
    outcome = result_path(tmp_path, 'source', 'geometry', 'forward')
    write_json(outcome, successor)
    pin = lambda p: dict(path=str(p.relative_to(tmp_path)), sha256=digest(p))
    key = 'source/geometry/forward'
    producers = {key: dict(previous_outcome=pin(old_path), run=pin(run_path), result_sha256=digest(outcome))}
    pins = {}
    assert view_producer(tmp_path, successor, {}, producers, pins) == run
    assert len(pins) == 2
    assert view_producer(tmp_path, successor, {'default': True}, {}, {}) == {'default': True}
    changed = deepcopy(successor)
    changed['request']['parameters']['max_records'] = 100
    with pytest.raises(ValueError, match='scientific inputs'):
        view_producer(tmp_path, changed, {}, producers, {})
    changed = deepcopy(successor)
    changed['request']['engine_id'] = 'undeclared-engine'
    with pytest.raises(ValueError, match='declared producer'):
        view_producer(tmp_path, changed, {}, producers, {})
    outcome.write_text('{}')
    with pytest.raises(ValueError, match='outcome changed'):
        view_producer(tmp_path, successor, {}, producers, {})
    old_path.write_text('{}')
    with pytest.raises(ValueError, match='producer control'):
        view_producer(tmp_path, successor, {}, producers, {})


def pilot(directory):
    rows = [dict(source_id="source", geometry_id="g", orientation=o,
                 status="not_assessable", acquisition_status="acquisition_error")
            for o in ("forward", "reverse")]
    data = dict(outcomes=rows, expected_geometry_views=400,
                inventory_sha256=digest(directory / "inventory.json.gz"),
                references_sha256=digest(directory / "references.json.gz"))
    path = directory / "pilot.json"
    write_json(path, data)
    return path


def test_frozen_original_pilot_accounts_for_all_failed_pairs():
    path = Path(__file__).parents[1] / "notes/data/sid-sv-pilot/coverage.json.gz"
    pairs = failed_pairs(read_json(path))
    assert len(pairs) == 22
    assert sum(len(row["original_outcomes"]) for row in pairs) == 44
    assert len({row["pair"] for row in pairs}) == 22
    assert {row["source_id"][:12] for row in pairs} == {"58f81602ec2f", "6e9147690525"}


def test_recovery_keeps_pilot_and_catalogue_denominators_separate(tmp_path, empty_audit):
    path = pilot(tmp_path)
    _, _, selection, _ = prepare(tmp_path / "recovery", tmp_path, path)
    assert selection["original_pilot_views"] == 2
    assert selection["original_catalogue_views"] == 400


def test_a_single_failed_orientation_cannot_shrink_the_scope():
    data = dict(outcomes=[dict(source_id="source", geometry_id="g", orientation="forward",
                              status="not_assessable")])
    with pytest.raises(ValueError, match="both failed orientations"):
        failed_pairs(data)


def test_prepare_keeps_pending_views_and_rejects_drift(tmp_path, empty_audit):
    directory = tmp_path / "recovery"
    path = pilot(tmp_path)
    _, _, selection, batches = prepare(directory, tmp_path, path)
    result = coverage(directory, selection, batches)
    assert result["expected_pairs"] == 1 and result["expected_views"] == 2
    assert result["acquisition_counts"] == {"pending": 2}
    assert result["outcome_counts"] == {"pending_reconstruction": 2}
    assert not result["all_inputs_acquired"] and not result["all_views_reconstructed"]
    assert not result["full_catalogue_complete"]
    assert prepare(directory, tmp_path, path)[2] == selection
    # Changing a source's original failure must not silently resume the old audit.
    changed = deepcopy(read_json(path))
    changed["outcomes"][0]["acquisition_status"] = "another_failure"
    other = tmp_path / "other-pilot.json"
    write_json(other, changed)
    with pytest.raises(ValueError, match="Recovery selection drift"):
        prepare(directory, tmp_path, other)


def test_real_completed_empty_bam_and_timeout_are_distinct(tmp_path, empty_audit):
    source, _, bam = empty_audit
    directory = tmp_path / "recovery"
    _, _, selection, batches = prepare(directory, tmp_path, pilot(tmp_path))
    batch = batches["g"]
    write_json(directory / "sources/source" / (batch["id"] + ".json"),
               dict(source_id="source", batch_id=batch["id"],
                    request=dict(source_identity=identity(source), batch=batch), status="acquired",
                    receipt=dict(records=0), unavailable_contigs=[],
                    bam=dict(path=str(bam), sha256=digest(bam)),
                    index=dict(path=str(bam) + ".bai", sha256=digest(str(bam) + ".bai"))))
    run_batch(directory, "source", batch["id"], dict(assemble=False), 30, implementation_identity())
    complete = coverage(directory, selection, batches)
    assert complete["all_inputs_acquired"] and complete["all_views_reconstructed"]
    assert complete["outcome_counts"] == {"no_candidate_paths": 2}
    path = directory / "results/source/g.forward.json.gz"
    result = read_json(path)
    details = directory / result["reconstruction"]["path"]
    details.write_bytes(b"interrupted or changed output")
    with pytest.raises(ValueError, match="changed completed reconstruction"):
        coverage(directory, selection, batches)
    # A recorded wall-time failure has no completed reconstruction or zero-support claim.
    path.unlink()
    write_json(path, dict(request=result["request"], status="reconstruction_timeout",
                          acquisition_status="acquired", input_limitations=[]))
    timed = coverage(directory, selection, batches)
    assert not timed["all_views_reconstructed"]
    forward, = [row for row in timed["outcomes"] if row["orientation"] == "forward"]
    assert forward["status"] == "reconstruction_timeout" and not forward["reconstructed"]
    assert forward["candidate_sequences"] is None


def test_cache_reuse_shares_read_bytes_but_isolates_metadata(tmp_path):
    source, destination = tmp_path / "original", tmp_path / "new-cache"
    source.mkdir()
    bam = source / "reads.bam"
    bam.write_bytes(b"original immutable BAM bytes")
    bam.chmod(0o444)
    metadata = source / "receipt.json"
    metadata.write_text('{"original":true}\n')
    clone_cache(source, destination)
    assert bam.stat().st_ino == (destination / "reads.bam").stat().st_ino
    (destination / "receipt.json").write_text('{"new":true}\n')
    assert metadata.read_text() == '{"original":true}\n'
    assert bam.read_bytes() == (destination / "reads.bam").read_bytes()
    with pytest.raises(ValueError, match="separate directory"):
        clone_cache(source, source / "nested")


def test_staged_run_preserves_pending_and_completed_coverage(tmp_path, monkeypatch, empty_audit):
    import examples.recover_sid_sv_pilot as runner

    source, _, bam = empty_audit
    path = pilot(tmp_path)
    directory = tmp_path / "recovery"

    def local_acquisition(directory, source, batch, *args, **kwargs):
        result = dict(source_id="source", batch_id=batch["id"],
                      request=dict(source_identity=identity(source), batch=batch), status="acquired",
                      receipt=dict(records=0), unavailable_contigs=[],
                      bam=dict(path=str(bam), sha256=digest(bam)),
                      index=dict(path=str(bam) + ".bai", sha256=digest(str(bam) + ".bai")))
        write_json(directory / "sources/source" / (batch["id"] + ".json"), result)
        return result

    monkeypatch.setattr(runner.acquisition, "acquire_batch", local_acquisition)
    args = ["recover", str(directory), "--base-audit", str(tmp_path), "--pilot-coverage", str(path),
            "--cache", str(tmp_path / "cache")]
    for stage in ("acquire", "reconstruct"):
        monkeypatch.setattr("sys.argv", [*args, "--stage", stage])
        runner.main()
    snapshots = [read_json(p) for p in directory.glob("recovery-coverage.*.json.gz")]
    assert len(snapshots) == 2
    assert {snapshot["all_views_reconstructed"] for snapshot in snapshots} == {True, False}
    assert all(snapshot["expected_views"] == 2 for snapshot in snapshots)
    assert all(row["original_outcomes"] for snapshot in snapshots for row in snapshot["outcomes"])


def test_independent_recovery_checker_requires_completed_views(tmp_path, empty_audit):
    from examples.check_sid_sv_recovery import verify

    source, _, bam = empty_audit
    directory = tmp_path / "recovery"
    _, _, _, batches = prepare(directory, tmp_path, pilot(tmp_path))
    batch = batches["g"]
    engine = implementation_identity()
    write_json(directory / "recovery-run.json", dict(engine_id=engine))
    with pytest.raises(ValueError, match="every selected input and reconstruction"):
        verify(directory, tmp_path / "checks")
    write_json(directory / "sources/source" / (batch["id"] + ".json"),
               dict(source_id="source", batch_id=batch["id"],
                    request=dict(source_identity=identity(source), batch=batch), status="acquired",
                    receipt=dict(records=0), unavailable_contigs=[],
                    bam=dict(path=str(bam), sha256=digest(bam)),
                    index=dict(path=str(bam) + ".bai", sha256=digest(str(bam) + ".bai"))))
    run_batch(directory, "source", batch["id"], dict(assemble=False, dense_support=True), 30, engine)
    ledger = verify(directory, tmp_path / "checks")
    assert ledger["intended_views"] == ledger["accounted_views"] == 2
    assert ledger["outcome_counts"] == {"no_candidate_paths": 2}
    assert ledger["rows"][0]["acquired_records"] == ledger["rows"][0]["regional_alignment_records"] == 0
    assert ledger["independent_orf_translations"] == 0
    assert ledger["independent_checks"] == []
    assert ledger["verified_reconstructed_views"] == 2 and ledger["unresolved_views"] == 0

    forward = directory / "results/source/g.forward.json.gz"
    completed = read_json(forward)
    failed = dict(request=completed["request"], acquisition_status="acquired",
                  input_limitations=[], status="reconstruction_timeout", error="wall-time budget exceeded")
    forward.unlink()
    write_json(forward, failed)
    with pytest.raises(ValueError, match="every selected input and reconstruction"):
        verify(directory, tmp_path / "strict-checks")
    partial = verify(directory, tmp_path / "partial-checks", allow_unresolved=True)
    assert partial["accounted_views"] == 2 and partial["verified_reconstructed_views"] == 1
    assert partial["unresolved_views"] == 1 and not partial["all_views_reconstructed"]
    assert partial["rows"][0]["views"]["forward"]["candidates"] is None
    assert partial["rows"][0]["views"]["reverse"]["candidates"] == 0
    assert partial["pins"]["results/source/g.forward.json.gz"] == digest(forward)
    forward.unlink()
    write_json(forward, dict(failed, candidates=[]))
    with pytest.raises(ValueError, match="unknown reconstruction/support counts"):
        verify(directory, tmp_path / "invalid-checks", allow_unresolved=True)
    forward.unlink()
    with pytest.raises(ValueError, match="recorded outcome for every selected view"):
        verify(directory, tmp_path / "missing-checks", allow_unresolved=True)


def test_mixed_checkout_is_rejected_before_preparation(tmp_path, monkeypatch, capsys):
    import examples.recover_sid_sv_pilot as runner

    monkeypatch.setattr(runner.isovar, "__file__", str(tmp_path / "other/isovar/__init__.py"))
    def unexpected_preparation(*args):
        raise AssertionError("Mixed package must be rejected before input preparation")
    monkeypatch.setattr(runner, "prepare", unexpected_preparation)
    monkeypatch.setattr("sys.argv", ["recover", str(tmp_path / "output"),
                                   "--base-audit", str(tmp_path / "base"),
                                   "--pilot-coverage", str(tmp_path / "pilot.json"),
                                   "--cache", str(tmp_path / "cache")])
    with pytest.raises(SystemExit) as error:
        runner.main()
    assert error.value.code == 2
    assert "same checkout" in capsys.readouterr().err
    assert not (tmp_path / "output").exists() and not (tmp_path / "cache").exists()
