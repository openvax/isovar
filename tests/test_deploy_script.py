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
import sys

import pytest


def _run(command, cwd, env=None):
    return subprocess.run(
        command,
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )


def _git(repo, *args):
    return subprocess.run(
        ["git"] + list(args),
        cwd=repo,
        capture_output=True,
        text=True,
        check=True,
    )


def _release_repo(tmp_path):
    source_root = Path(__file__).resolve().parents[1]
    repo = tmp_path / "release-repo"
    package_dir = repo / "isovar"
    package_dir.mkdir(parents=True)
    shutil.copy(source_root / "deploy.sh", repo / "deploy.sh")
    (package_dir / "__init__.py").write_text('__version__ = "1.7.2"\n')
    (repo / ".gitignore").write_text("__pycache__/\nbuild/\ndist/\n")
    for script_name, output in (("lint.sh", "lint-gate"), ("test.sh", "test-gate")):
        script = repo / script_name
        script.write_text("#!/bin/sh\necho %s\n" % output)
        script.chmod(0o755)

    _git(repo, "init")
    _git(repo, "checkout", "-b", "main")
    _git(repo, "config", "user.name", "Test Release")
    _git(repo, "config", "user.email", "release@example.com")
    _git(repo, "add", ".")
    _git(repo, "commit", "-m", "Initial release fixture")
    origin = tmp_path / "origin.git"
    subprocess.run(
        ["git", "init", "--bare", str(origin)],
        capture_output=True,
        text=True,
        check=True,
    )
    _git(repo, "remote", "add", "origin", str(origin))
    _git(repo, "push", "--set-upstream", "origin", "main")
    return repo


def _deploy(repo, *args, env_updates=None):
    env = os.environ.copy()
    for name in tuple(env):
        if name.startswith("COV_CORE_"):
            del env[name]
    env["PYTHON"] = sys.executable
    if env_updates:
        env.update(env_updates)
    return _run(["bash", "deploy.sh"] + list(args), cwd=repo, env=env)


def _fake_release_python(tmp_path):
    fake_python = tmp_path / "fake-python"
    fake_python.write_text("""#!/bin/bash
set -eu
if [[ "${1:-}" == "-" ]]; then
  exec "${REAL_PYTHON}" "$@"
elif [[ "${1:-}" == "-m" && "${2:-}" == "pip" ]]; then
  exit 0
elif [[ "${1:-}" == "-m" && "${2:-}" == "build" ]]; then
  version="$("${REAL_PYTHON}" -c \
    'import re; from pathlib import Path; print(re.search(r"__version__ = \\\"([^\\\"]+)\\\"", Path("isovar/__init__.py").read_text()).group(1))')"
  mkdir -p dist
  touch "dist/isovar-${version}-py3-none-any.whl"
  touch "dist/isovar-${version}.tar.gz"
  echo "built ${version}"
  exit 0
elif [[ "${1:-}" == "-m" && "${2:-}" == "twine" ]]; then
  echo "twine-check"
  exit 0
elif [[ "${1:-}" == "release_upload.py" ]]; then
  echo "release-upload $*"
  exit 0
fi
exec "${REAL_PYTHON}" "$@"
""")
    fake_python.chmod(0o755)
    return fake_python


def test_deploy_dry_run_validates_without_mutating(tmp_path):
    repo = _release_repo(tmp_path)

    result = _deploy(repo, "--dry-run", "v1.7.3")

    assert result.returncode == 0, result.stderr
    assert "lint-gate" in result.stdout
    assert "test-gate" in result.stdout
    assert "would bump 1.7.2 to 1.7.3" in result.stdout
    assert "no files or remote state changed" in result.stdout
    assert (repo / "isovar" / "__init__.py").read_text() == \
        '__version__ = "1.7.2"\n'
    assert _git(repo, "rev-list", "--count", "HEAD").stdout.strip() == "1"
    assert _git(repo, "status", "--porcelain").stdout == ""


def test_deploy_rejects_non_release_branch_before_gates(tmp_path):
    repo = _release_repo(tmp_path)
    _git(repo, "checkout", "-b", "feature")

    result = _deploy(repo, "--dry-run", "1.7.3")

    assert result.returncode == 1
    assert "only allowed from main or master" in result.stderr
    assert "lint-gate" not in result.stdout


def test_deploy_rejects_dirty_tree_before_gates(tmp_path):
    repo = _release_repo(tmp_path)
    (repo / "isovar" / "__init__.py").write_text('__version__ = "dirty"\n')

    result = _deploy(repo, "--dry-run", "1.7.3")

    assert result.returncode == 1
    assert "Working tree not clean" in result.stderr
    assert "lint-gate" not in result.stdout


def test_deploy_rejects_existing_local_tag_before_gates(tmp_path):
    repo = _release_repo(tmp_path)
    _git(repo, "tag", "v1.7.3")

    result = _deploy(repo, "--dry-run", "1.7.3")

    assert result.returncode == 1
    assert "Tag v1.7.3 already exists locally at a different release commit" in \
        result.stderr
    assert "lint-gate" not in result.stdout


@pytest.mark.parametrize("version", [
    "definitely-not-a-version",
    "01.7.3",
    "1.07.3",
    "1.7.03",
    "1.7.3rc1",
])
def test_deploy_rejects_invalid_version_before_gates(tmp_path, version):
    repo = _release_repo(tmp_path)

    result = _deploy(repo, "--dry-run", version)

    assert result.returncode == 1
    assert "Invalid release version" in result.stderr
    assert "lint-gate" not in result.stdout


def test_deploy_rejects_version_downgrade_before_gates(tmp_path):
    repo = _release_repo(tmp_path)

    result = _deploy(repo, "1.7.1")

    assert result.returncode == 1
    assert "Release version must increase: 1.7.2 -> 1.7.1" in result.stderr
    assert "lint-gate" not in result.stdout


@pytest.mark.parametrize("script_name", ["lint.sh", "test.sh"])
def test_deploy_stops_before_mutation_when_local_gate_fails(
        tmp_path,
        script_name):
    repo = _release_repo(tmp_path)
    gate_script = repo / script_name
    gate_script.write_text("#!/bin/sh\nexit 7\n")
    _git(repo, "add", script_name)
    _git(repo, "commit", "-m", "Make fixture gate fail")
    _git(repo, "push", "origin", "main")

    result = _deploy(repo, "1.7.3")

    assert result.returncode == 7
    assert (repo / "isovar" / "__init__.py").read_text() == \
        '__version__ = "1.7.2"\n'
    assert _git(repo, "rev-list", "--count", "HEAD").stdout.strip() == "2"
    assert _git(repo, "tag", "--list", "v1.7.3").stdout == ""


def test_deploy_rejects_existing_remote_tag_before_mutation(tmp_path):
    repo = _release_repo(tmp_path)
    _git(repo, "tag", "v1.7.3")
    _git(repo, "push", "origin", "v1.7.3")
    _git(repo, "tag", "--delete", "v1.7.3")

    result = _deploy(repo, "1.7.3")

    assert result.returncode == 1
    assert "Tag v1.7.3 already exists on origin" in result.stderr
    assert (repo / "isovar" / "__init__.py").read_text() == \
        '__version__ = "1.7.2"\n'
    assert _git(repo, "rev-list", "--count", "HEAD").stdout.strip() == "1"


def test_deploy_fails_closed_when_remote_tags_cannot_be_checked(tmp_path):
    repo = _release_repo(tmp_path)
    _git(repo, "remote", "set-url", "origin", str(tmp_path / "missing.git"))

    result = _deploy(repo, "1.7.3")

    assert result.returncode == 1
    assert "Could not determine whether v1.7.3 exists on origin" in result.stderr
    assert (repo / "isovar" / "__init__.py").read_text() == \
        '__version__ = "1.7.2"\n'
    assert _git(repo, "rev-list", "--count", "HEAD").stdout.strip() == "1"


def test_deploy_resumes_when_valid_local_tag_has_not_reached_origin(tmp_path):
    repo = _release_repo(tmp_path)
    version_path = repo / "isovar" / "__init__.py"
    version_path.write_text('__version__ = "1.7.3"\n')
    _git(repo, "add", "isovar/__init__.py")
    _git(repo, "commit", "-m", "Bump version to 1.7.3")
    _git(repo, "push", "origin", "main")
    _git(repo, "tag", "v1.7.3")

    fake_python = _fake_release_python(tmp_path)

    result = _deploy(
        repo,
        "1.7.3",
        env_updates={
            "PYTHON": str(fake_python),
            "REAL_PYTHON": sys.executable,
        })

    assert result.returncode == 0, result.stderr
    assert "Resuming release from existing local tag v1.7.3" in result.stdout
    assert "refs/tags/v1.7.3" in _git(
        repo, "ls-remote", "--tags", "origin", "refs/tags/v1.7.3").stdout
    assert _git(repo, "rev-list", "--count", "HEAD").stdout.strip() == "2"


def test_deploy_resumes_after_version_commit_push_fails(tmp_path):
    repo = _release_repo(tmp_path)
    origin = Path(_git(repo, "remote", "get-url", "origin").stdout.strip())
    rejecting_hook = origin / "hooks" / "pre-receive"
    rejecting_hook.write_text("#!/bin/sh\nexit 1\n")
    rejecting_hook.chmod(0o755)

    first_attempt = _deploy(repo, "1.7.3")

    assert first_attempt.returncode != 0
    assert (repo / "isovar" / "__init__.py").read_text() == \
        '__version__ = "1.7.3"\n'
    assert '__version__ = "1.7.2"' in _git(
        repo, "show", "origin/main:isovar/__init__.py").stdout
    assert _git(repo, "tag", "--list", "v1.7.3").stdout == ""

    rejecting_hook.unlink()
    fake_python = _fake_release_python(tmp_path)
    second_attempt = _deploy(
        repo,
        "1.7.3",
        env_updates={
            "PYTHON": str(fake_python),
            "REAL_PYTHON": sys.executable,
        })

    assert second_attempt.returncode == 0, second_attempt.stderr
    assert _git(repo, "rev-list", "--count", "HEAD").stdout.strip() == "2"
    assert '__version__ = "1.7.3"' in _git(
        repo, "show", "origin/main:isovar/__init__.py").stdout
    assert "refs/tags/v1.7.3" in _git(
        repo, "ls-remote", "--tags", "origin", "refs/tags/v1.7.3").stdout


def test_deploy_bumps_pushes_builds_uploads_and_tags(tmp_path):
    repo = _release_repo(tmp_path)
    fake_python = _fake_release_python(tmp_path)

    result = _deploy(
        repo,
        "1.7.3",
        env_updates={
            "PYTHON": str(fake_python),
            "REAL_PYTHON": sys.executable,
        })

    assert result.returncode == 0, result.stderr
    assert "built 1.7.3" in result.stdout
    assert "twine-check" in result.stdout
    assert "release-upload" in result.stdout
    assert "Deployed isovar 1.7.3" in result.stdout
    assert (repo / "isovar" / "__init__.py").read_text() == \
        '__version__ = "1.7.3"\n'
    assert _git(repo, "log", "-1", "--pretty=%s").stdout.strip() == \
        "Bump version to 1.7.3"
    assert '__version__ = "1.7.3"' in _git(
        repo, "show", "origin/main:isovar/__init__.py").stdout
    assert _git(repo, "tag", "--list", "v1.7.3").stdout.strip() == "v1.7.3"
    assert "refs/tags/v1.7.3" in _git(
        repo, "ls-remote", "--tags", "origin", "refs/tags/v1.7.3").stdout
    assert _git(repo, "status", "--porcelain").stdout == ""
