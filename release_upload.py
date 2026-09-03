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

"""Bounded, idempotent PyPI uploads for isovar releases."""

import argparse
import json
import os
from pathlib import Path
from secrets import token_hex
import subprocess
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import urlopen


PYPI_JSON_BASE_URL = "https://pypi.org/pypi"
DEFAULT_UPLOAD_TIMEOUT_SECONDS = 60.0
DEFAULT_VERIFY_TIMEOUT_SECONDS = 30.0
DEFAULT_VERIFY_POLL_SECONDS = 1.0


class ReleaseUploadError(RuntimeError):
    """Raised when the complete release artifact set cannot be verified."""


def expected_release_filenames(project, version):
    """Return the expected pure-Python wheel and source-distribution names."""
    wheel_project = project.replace("-", "_")
    return frozenset({
        "%s-%s-py3-none-any.whl" % (wheel_project, version),
        "%s-%s.tar.gz" % (project, version),
    })


def pypi_release_filenames(
        project,
        version,
        json_base_url=PYPI_JSON_BASE_URL,
        request_timeout_seconds=10.0):
    """Return the filenames published for one release in PyPI metadata."""
    url = "%s/%s/json?cache_bust=%s" % (
        json_base_url.rstrip("/"),
        quote(project, safe=""),
        token_hex(8),
    )
    try:
        with urlopen(url, timeout=request_timeout_seconds) as response:
            payload = json.load(response)
    except HTTPError as error:
        if error.code == 404:
            return frozenset()
        raise
    releases = payload.get("releases", {})
    return frozenset(item["filename"] for item in releases.get(version, ()))


def upload_distribution(
        distribution_path,
        python_executable=sys.executable,
        timeout_seconds=DEFAULT_UPLOAD_TIMEOUT_SECONDS,
        repository_url=None):
    """Upload one distribution with a finite timeout and quiet progress."""
    command = [
        python_executable,
        "-m",
        "twine",
        "upload",
        "--disable-progress-bar",
    ]
    if repository_url:
        command.extend(["--repository-url", repository_url])
    command.append(str(distribution_path))
    subprocess.run(command, check=True, timeout=timeout_seconds)


def wait_for_release_file(
        filename,
        fetch_release_filenames,
        timeout_seconds=DEFAULT_VERIFY_TIMEOUT_SECONDS,
        poll_seconds=DEFAULT_VERIFY_POLL_SECONDS):
    """Poll release metadata until ``filename`` is visible or time expires."""
    deadline = time.monotonic() + timeout_seconds
    while True:
        try:
            if filename in set(fetch_release_filenames()):
                return True
        except (json.JSONDecodeError, OSError, URLError):
            # PyPI metadata can transiently fail while its release view is
            # converging. Keep the retry bounded by the same deadline.
            pass
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return False
        time.sleep(min(poll_seconds, remaining))


def publish_release(
        distribution_paths,
        expected_filenames,
        fetch_release_filenames,
        upload_file,
        verify_timeout_seconds=DEFAULT_VERIFY_TIMEOUT_SECONDS,
        verify_poll_seconds=DEFAULT_VERIFY_POLL_SECONDS):
    """Upload missing artifacts and reconcile every response with PyPI state.

    A timeout or nonzero Twine exit is ambiguous: the server may have accepted
    the immutable file before the client lost its response. Verify server state
    after every attempt and consider a file published only when its exact name
    appears in PyPI metadata.
    """
    distribution_paths = tuple(distribution_paths)
    paths_by_name = {
        Path(distribution_path).name: Path(distribution_path)
        for distribution_path in distribution_paths
    }
    expected = frozenset(expected_filenames)
    if (len(paths_by_name) != len(distribution_paths) or
            frozenset(paths_by_name) != expected):
        raise ReleaseUploadError(
            "Distribution files do not match the expected release set: "
            "expected=%s, observed=%s" % (
                sorted(expected),
                sorted(paths_by_name),
            )
        )

    published = frozenset(fetch_release_filenames())
    for filename in sorted(expected):
        if filename in published:
            print("Already published: %s" % filename)
            continue

        upload_error = None
        try:
            upload_file(paths_by_name[filename])
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
            upload_error = error

        if not wait_for_release_file(
                filename,
                fetch_release_filenames,
                timeout_seconds=verify_timeout_seconds,
                poll_seconds=verify_poll_seconds):
            message = "PyPI did not publish %s after the upload attempt" % filename
            if upload_error is not None:
                raise ReleaseUploadError(message) from upload_error
            raise ReleaseUploadError(message)
        if upload_error is not None:
            print("Reconciled ambiguous upload from PyPI state: %s" % filename)
        else:
            print("Published: %s" % filename)
        published = frozenset(fetch_release_filenames())

    published = frozenset(fetch_release_filenames())
    missing = expected - published
    if missing:
        raise ReleaseUploadError(
            "PyPI release is missing expected files: %s" % sorted(missing)
        )
    return published


def main(argv=None):
    """Upload and verify a complete isovar release."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--json-base-url", default=PYPI_JSON_BASE_URL)
    parser.add_argument("--repository-url")
    parser.add_argument(
        "--upload-timeout-seconds",
        type=float,
        default=float(os.environ.get(
            "ISOVAR_UPLOAD_TIMEOUT_SECONDS",
            DEFAULT_UPLOAD_TIMEOUT_SECONDS,
        )),
    )
    parser.add_argument(
        "--verify-timeout-seconds",
        type=float,
        default=float(os.environ.get(
            "ISOVAR_VERIFY_TIMEOUT_SECONDS",
            DEFAULT_VERIFY_TIMEOUT_SECONDS,
        )),
    )
    parser.add_argument("distributions", nargs="+")
    args = parser.parse_args(argv)

    expected = expected_release_filenames(args.project, args.version)

    def fetch_release_filenames():
        return pypi_release_filenames(
            args.project,
            args.version,
            json_base_url=args.json_base_url,
        )

    def upload_file(path):
        upload_distribution(
            path,
            python_executable=sys.executable,
            timeout_seconds=args.upload_timeout_seconds,
            repository_url=args.repository_url,
        )

    published = publish_release(
        args.distributions,
        expected_filenames=expected,
        fetch_release_filenames=fetch_release_filenames,
        upload_file=upload_file,
        verify_timeout_seconds=args.verify_timeout_seconds,
    )
    print(
        "Verified PyPI release files: %s" %
        ", ".join(sorted(expected & published))
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
