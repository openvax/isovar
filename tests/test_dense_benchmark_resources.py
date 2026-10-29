"""A completed reconstruction is not a successfully exported result."""

import json
from pathlib import Path
import signal
import socket
import time
from types import SimpleNamespace

import pytest

from examples import benchmark_dense_sv as driver
from isovar import rna_library, sv_rna_dense, sv_rna_stream
from isovar.sid_audit_bundle import bundle_sid_audit
from tests.test_audit_bundle import original_sid_audit


@pytest.mark.parametrize('phase', ['fingerprint', 'serialization'])
def test_resource_guard_covers_export_and_marks_partial_results(tmp_path, request, monkeypatch, phase):
    bundle = bundle_sid_audit(request.getfixturevalue('original_sid_audit')[0], tmp_path / 'bundle')
    run_id = next(iter(bundle.manifest['metadata']['sv_rna_runs']))
    output = tmp_path / 'metrics.json'
    monkeypatch.setattr('sys.argv', ['benchmark_dense_sv', '--checkout', str(Path(driver.__file__).resolve().parents[1]),
                                   '--bundle', str(bundle.directory), '--run', run_id, '--max-records', '50000',
                                   '--seconds', '30', '--min-free-bytes', '1', '--output', str(output)])
    exporting = [False]
    monkeypatch.setattr(driver.shutil, 'disk_usage', lambda path: SimpleNamespace(free=0 if exporting[0] else 100))
    original_fingerprint = driver.fingerprint

    def fingerprint(value):
        if phase == 'fingerprint':
            exporting[0] = True
            time.sleep(2)
        return original_fingerprint(value)

    def write_partial(value, handle, **kwargs):
        handle.write('{"partial":')
        exporting[0] = True
        time.sleep(2)

    monkeypatch.setattr(driver, 'fingerprint', fingerprint)
    monkeypatch.setattr(driver.json, 'dump', write_partial)
    # main deliberately installs process-local instrumentation and offline
    # guards; restore these when exercising it inside the test process.
    for obj, names in ((socket, ('create_connection',)), (socket.socket, ('connect',)),
                       (rna_library, ('prepare_rna_library',)), (sv_rna_dense, ('_select', '_support')), (sv_rna_stream._SupportStream, ('consume',))):
        for name in names:
            monkeypatch.setattr(obj, name, getattr(obj, name))
    previous = signal.getsignal(signal.SIGALRM)
    try:
        assert driver.main() == 2
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)
    report = json.loads(output.read_text())
    assert report['reconstruction_status'] == 'completed'
    assert report['status'] == 'result_free_space_limit'
    assert not Path(str(output) + '.result.json.gz').exists()
    if phase == 'serialization':
        assert Path(report['partial_result']).is_file()
    else:
        assert 'partial_result' not in report
