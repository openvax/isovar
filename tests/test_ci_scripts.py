"""CI setup retries transient Ensembl release-asset failures."""

import os
from pathlib import Path
import subprocess


SOURCE_ROOT = Path(__file__).resolve().parents[1]
INSTALL_SCRIPT = SOURCE_ROOT / ".github/scripts/install_ensembl_data.sh"


def run_installer(tmp_path, always_fail_release=""):
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    pyensembl = fake_bin / "pyensembl"
    pyensembl.write_text("""#!/usr/bin/env bash
set -eu
release=""
while [[ $# -gt 0 ]]; do
  if [[ "$1" == "--release" ]]; then
    release="$2"
    break
  fi
  shift
done
count_file="$FAKE_STATE/$release"
count=0
if [[ -f "$count_file" ]]; then
  count="$(<"$count_file")"
fi
count=$((count + 1))
printf '%s' "$count" > "$count_file"
printf '%s:%s\n' "$release" "$count" >> "$INSTALL_LOG"
if [[ "$release" == "$ALWAYS_FAIL_RELEASE" || "$count" -eq 1 ]]; then
  exit 1
fi
""")
    pyensembl.chmod(0o755)
    sleep = fake_bin / "sleep"
    sleep.write_text("""#!/usr/bin/env bash
printf '%s\n' "$1" >> "$SLEEP_LOG"
""")
    sleep.chmod(0o755)
    state = tmp_path / "state"
    state.mkdir()
    install_log = tmp_path / "install.log"
    sleep_log = tmp_path / "sleep.log"
    env = os.environ.copy()
    env.update({
        "ALWAYS_FAIL_RELEASE": always_fail_release,
        "FAKE_STATE": str(state),
        "INSTALL_LOG": str(install_log),
        "PATH": f"{fake_bin}{os.pathsep}{env['PATH']}",
        "SLEEP_LOG": str(sleep_log),
    })
    result = subprocess.run(
        ["bash", str(INSTALL_SCRIPT)],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    return result, install_log.read_text().splitlines(), sleep_log.read_text().splitlines()


def test_ensembl_installer_retries_each_release(tmp_path):
    result, installs, sleeps = run_installer(tmp_path)

    assert result.returncode == 0, result.stderr
    assert installs == ["87:1", "87:2", "75:1", "75:2", "102:1", "102:2"]
    assert sleeps == ["15", "15", "15"]


def test_ensembl_installer_stops_after_bounded_attempts(tmp_path):
    result, installs, sleeps = run_installer(tmp_path, always_fail_release="87")

    assert result.returncode == 1
    assert installs == ["87:1", "87:2", "87:3"]
    assert sleeps == ["15", "30"]
    assert "Failed to install Ensembl release 87 after 3 attempts" in result.stderr
