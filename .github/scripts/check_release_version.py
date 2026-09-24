"""Fail a pull request whose version would collide with a released or base version.

Every PR bumps isovar.__version__ (see AGENTS.md). Two PRs from the same base
can carry the same bump; once one is released, the other must not merge with
that version (#345). Run on the PR's merge commit:

    python .github/scripts/check_release_version.py BASE_VERSION_FILE

where BASE_VERSION_FILE is the base branch's isovar/__init__.py. Released
versions come from PyPI and the repository's v* tags, which are fetched with
``git ls-remote``; an unreachable index is an error, never a pass.
"""

import json
from pathlib import Path
import re
import subprocess
import sys
from urllib.request import urlopen

from packaging.version import Version

ROOT = Path(__file__).resolve().parents[2]


def version_of(init_text):
    match = re.search(r'^__version__ = "([^"]+)"$', init_text, re.M)
    if match is None:
        raise ValueError("No __version__ assignment found")
    return Version(match.group(1))


def released_versions():
    with urlopen("https://pypi.org/pypi/isovar/json", timeout=30) as response:
        versions = {Version(v) for v in json.load(response)["releases"]}
    tags = subprocess.run(["git", "ls-remote", "--tags", "origin", "refs/tags/v*"],
                          cwd=ROOT, check=True, capture_output=True, text=True).stdout
    versions.update(Version(ref.rsplit("/v", 1)[1]) for ref in re.findall(r"refs/tags/v[^\s^]+", tags))
    return versions


def problems(version, base_version, released):
    found = []
    if version <= base_version:
        found.append("version %s must be greater than the base branch's %s" % (version, base_version))
    if version in released:
        found.append("version %s is already released or tagged" % version)
    return found


def main(base_version_file):
    version = version_of((ROOT / "isovar/__init__.py").read_text())
    base_version = version_of(Path(base_version_file).read_text())
    found = problems(version, base_version, released_versions())
    for problem in found:
        print("::error::" + problem + "; bump isovar/__init__.py")
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
