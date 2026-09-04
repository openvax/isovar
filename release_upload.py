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
import hashlib
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
SHA256_HEX_LENGTH = 64
HASH_READ_SIZE_BYTES = 1024 * 1024


class ReleaseUploadError(RuntimeError):
    """Raised when the complete release artifact set cannot be verified."""


def expected_release_filenames(project, version):
    """Return the expected pure-Python wheel and source-distribution names."""
    wheel_project = project.replace("-", "_")
    return frozenset({
        "%s-%s-py3-none-any.whl" % (wheel_project, version),
        "%s-%s.tar.gz" % (project, version),
    })


def sha256_file(path):
    """Return the hexadecimal SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as file_handle:
        while True:
            chunk = file_handle.read(HASH_READ_SIZE_BYTES)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _validated_sha256(value, context):
    """Normalize a SHA-256 value or fail closed on malformed metadata."""
    if not isinstance(value, str):
        raise ReleaseUploadError("Missing SHA-256 digest for %s" % context)
    normalized = value.lower()
    if (len(normalized) != SHA256_HEX_LENGTH or
            any(character not in "0123456789abcdef" for character in normalized)):
        raise ReleaseUploadError(
            "Invalid SHA-256 digest for %s: %r" % (context, value)
        )
    return normalized


def pypi_release_digests(
        project,
        version,
        json_base_url=PYPI_JSON_BASE_URL,
        request_timeout_seconds=10.0):
    """Return a filename-to-SHA-256 map from PyPI release metadata."""
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
            return {}
        raise
    if not isinstance(payload, dict):
        raise ReleaseUploadError(
            "Invalid PyPI release metadata for %s %s" % (project, version)
        )
    releases = payload.get("releases", {})
    if not isinstance(releases, dict):
        raise ReleaseUploadError(
            "Invalid PyPI release metadata for %s %s" % (project, version)
        )
    release_items = releases.get(version, ())
    if not isinstance(release_items, (list, tuple)):
        raise ReleaseUploadError(
            "Invalid PyPI release metadata for %s %s" % (project, version)
        )
    result = {}
    for item in release_items:
        try:
            filename = item["filename"]
            sha256 = item["digests"]["sha256"]
        except (KeyError, TypeError) as error:
            raise ReleaseUploadError(
                "Invalid PyPI artifact metadata for %s %s" % (project, version)
            ) from error
        if not isinstance(filename, str) or not filename:
            raise ReleaseUploadError(
                "Invalid PyPI artifact filename for %s %s: %r" % (
                    project,
                    version,
                    filename,
                )
            )
        sha256 = _validated_sha256(sha256, "PyPI artifact %s" % filename)
        previous = result.get(filename)
        if previous is not None and previous != sha256:
            raise ReleaseUploadError(
                "Conflicting PyPI SHA-256 digests for %s" % filename
            )
        result[filename] = sha256
    return result


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


def _published_file_matches(filename, expected_sha256, published_digests):
    """Return true for an exact file, false if absent, and fail if different."""
    expected_sha256 = _validated_sha256(
        expected_sha256,
        "local artifact %s" % filename,
    )
    if filename not in published_digests:
        return False
    actual_sha256 = _validated_sha256(
        published_digests[filename],
        "PyPI artifact %s" % filename,
    )
    if actual_sha256 != expected_sha256:
        raise ReleaseUploadError(
            "PyPI SHA-256 mismatch for %s: local=%s, pypi=%s" % (
                filename,
                expected_sha256,
                actual_sha256,
            )
        )
    return True


def wait_for_release_file(
        filename,
        expected_sha256,
        fetch_release_digests,
        timeout_seconds=DEFAULT_VERIFY_TIMEOUT_SECONDS,
        poll_seconds=DEFAULT_VERIFY_POLL_SECONDS):
    """Poll until ``filename`` is visible with the expected SHA-256 digest."""
    deadline = time.monotonic() + timeout_seconds
    while True:
        try:
            published = dict(fetch_release_digests())
            if _published_file_matches(
                    filename,
                    expected_sha256,
                    published):
                return True
        except (json.JSONDecodeError, OSError, URLError):
            # PyPI metadata can transiently fail while its release view is
            # converging. Keep the retry bounded by the same deadline.
            pass
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return False
        time.sleep(min(poll_seconds, remaining))


def _matching_release_filenames(expected_digests, published_digests):
    """Return exact matches and fail if any expected filename has new bytes."""
    matching = set()
    for filename, expected_sha256 in expected_digests.items():
        if _published_file_matches(
                filename,
                expected_sha256,
                published_digests):
            matching.add(filename)
    return frozenset(matching)


def publish_release(
        distribution_paths,
        expected_filenames,
        fetch_release_digests,
        upload_file,
        verify_timeout_seconds=DEFAULT_VERIFY_TIMEOUT_SECONDS,
        verify_poll_seconds=DEFAULT_VERIFY_POLL_SECONDS):
    """Upload missing artifacts and reconcile every response with PyPI state.

    A timeout or nonzero Twine exit is ambiguous: the server may have accepted
    the immutable file before the client lost its response. Verify server state
    after every attempt and consider a file published only when its exact name
    appears in PyPI metadata with the exact local SHA-256 digest.
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

    try:
        expected_digests = {
            filename: sha256_file(paths_by_name[filename])
            for filename in sorted(expected)
        }
    except OSError as error:
        raise ReleaseUploadError(
            "Could not hash release artifact: %s" % error
        ) from error

    for filename in sorted(expected):
        published_digests = dict(fetch_release_digests())
        matching = _matching_release_filenames(
            expected_digests,
            published_digests,
        )
        if filename in matching:
            print("Already published with matching SHA-256: %s" % filename)
            continue

        upload_error = None
        try:
            upload_file(paths_by_name[filename])
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
            upload_error = error

        if not wait_for_release_file(
                filename,
                expected_digests[filename],
                fetch_release_digests,
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

    published_digests = dict(fetch_release_digests())
    matching = _matching_release_filenames(
        expected_digests,
        published_digests,
    )
    missing = expected - matching
    if missing:
        raise ReleaseUploadError(
            "PyPI release is missing expected files: %s" % sorted(missing)
        )
    return published_digests


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

    def fetch_release_digests():
        return pypi_release_digests(
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
        fetch_release_digests=fetch_release_digests,
        upload_file=upload_file,
        verify_timeout_seconds=args.verify_timeout_seconds,
    )
    print(
        "Verified PyPI release files: %s" %
        ", ".join(sorted(expected & set(published)))
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
