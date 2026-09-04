# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import os
from pathlib import Path
import shutil
import subprocess

import pytest


SOURCE_ROOT = Path(__file__).resolve().parents[1]


def _fake_python(python_path):
    python_path.parent.mkdir(parents=True, exist_ok=True)
    python_path.write_text("""#!/usr/bin/env bash
printf '%s\\n' "$*" >> "$PYTHON_INVOCATION_LOG"
if [[ "${1:-}" == "-c" && "${2:-}" == "import xdist" ]]; then
    exit "${XDIST_IMPORT_STATUS:-0}"
fi
exit 0
""")
    python_path.chmod(0o755)
    return python_path


def _path_traps(tmp_path):
    trap_dir = tmp_path / "path-traps"
    trap_dir.mkdir()
    for command_name in ("python", "python3", "pytest", "ruff"):
        command_path = trap_dir / command_name
        command_path.write_text("""#!/bin/sh
printf '%s\\n' "$0 $*" >> "$PATH_TRAP_LOG"
exit 97
""")
        command_path.chmod(0o755)
    return trap_dir


def _run_gate(tmp_path, script_name, *args, env_updates=None):
    script_path = tmp_path / script_name
    shutil.copy2(SOURCE_ROOT / script_name, script_path)
    env = os.environ.copy()
    env.update({
        "PATH": "%s%s%s" % (
            _path_traps(tmp_path),
            os.pathsep,
            env["PATH"],
        ),
        "PATH_TRAP_LOG": str(tmp_path / "path-trap.log"),
        "PYTHON": str(_fake_python(tmp_path / "selected python")),
        "PYTHON_INVOCATION_LOG": str(tmp_path / "python-invocations.log"),
        "TEST_SH_MAX": "1",
        "TEST_SH_MIN": "1",
    })
    if env_updates:
        for name, value in env_updates.items():
            if value is None:
                env.pop(name, None)
            else:
                env[name] = value
    result = subprocess.run(
        ["bash", str(script_path), *args],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    invocation_log = tmp_path / "python-invocations.log"
    invocations = (
        invocation_log.read_text().splitlines()
        if invocation_log.exists()
        else [])
    return result, invocations, tmp_path / "path-trap.log"


def test_lint_uses_selected_python_module_not_path_ruff(tmp_path):
    result, invocations, trap_log = _run_gate(tmp_path, "lint.sh")

    assert result.returncode == 0, result.stderr
    assert invocations == ["-m ruff check isovar/ tests/"]
    assert not trap_log.exists()


@pytest.mark.parametrize(
    "script_name, expected_invocations",
    [
        ("lint.sh", ["-m ruff check isovar/ tests/"]),
        ("test.sh", [
            "-c import xdist",
            "-m pytest -n 1 --cov=isovar/ --cov-report=term-missing tests",
        ]),
    ],
)
@pytest.mark.parametrize("environment", ["active-venv", "repo-venv"])
def test_gate_resolves_virtualenv_python_when_python_is_unset(
        tmp_path,
        script_name,
        expected_invocations,
        environment):
    if environment == "active-venv":
        venv = tmp_path / "active venv"
        env_updates = {
            "PYTHON": None,
            "VIRTUAL_ENV": str(venv),
        }
    else:
        venv = tmp_path / ".venv"
        env_updates = {
            "PYTHON": None,
            "VIRTUAL_ENV": None,
        }
    selected_python = _fake_python(venv / "bin" / "python")

    result, invocations, trap_log = _run_gate(
        tmp_path,
        script_name,
        env_updates=env_updates)

    assert result.returncode == 0, result.stderr
    assert invocations == expected_invocations
    if script_name == "test.sh":
        assert "python=%s" % selected_python in result.stderr
    assert not trap_log.exists()


@pytest.mark.parametrize(
    "xdist_import_status, expected_pytest_invocation",
    [
        ("0", "-m pytest -n 1 --cov=isovar/ --cov-report=term-missing "
              "tests selected-test"),
        ("1", "-m pytest --cov=isovar/ --cov-report=term-missing "
              "tests selected-test"),
    ],
)
def test_test_gate_uses_same_python_for_plugin_probe_and_pytest(
        tmp_path,
        xdist_import_status,
        expected_pytest_invocation):
    result, invocations, trap_log = _run_gate(
        tmp_path,
        "test.sh",
        "selected-test",
        env_updates={"XDIST_IMPORT_STATUS": xdist_import_status})

    assert result.returncode == 0, result.stderr
    assert invocations == [
        "-c import xdist",
        expected_pytest_invocation,
    ]
    assert not trap_log.exists()
