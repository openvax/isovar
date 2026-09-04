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


import io
import json
from pathlib import Path
import subprocess
from urllib.error import HTTPError, URLError

import pytest

import release_upload
from release_upload import (
    ReleaseUploadError,
    expected_release_filenames,
    main,
    pypi_release_digests,
    publish_release,
    sha256_file,
    upload_distribution,
    wait_for_release_file,
)


EXPECTED = expected_release_filenames("isovar", "1.7.3")
WHEEL = "isovar-1.7.3-py3-none-any.whl"
SDIST = "isovar-1.7.3.tar.gz"
PATHS = tuple(Path("dist") / filename for filename in EXPECTED)
VALID_SHA256 = "a" * 64
OTHER_SHA256 = "b" * 64


@pytest.fixture
def release_files(tmp_path):
    dist_dir = tmp_path / "dist"
    dist_dir.mkdir()
    contents = {
        WHEEL: b"wheel bytes",
        SDIST: b"source distribution bytes",
    }
    paths = []
    for filename in sorted(EXPECTED):
        path = dist_dir / filename
        path.write_bytes(contents[filename])
        paths.append(path)
    return tuple(paths)


def _local_digests(paths):
    return {path.name: sha256_file(path) for path in paths}


def test_expected_release_filenames_are_exact_wheel_and_sdist():
    assert EXPECTED == {WHEEL, SDIST}


def test_sha256_file_reads_the_complete_file_in_chunks(tmp_path, monkeypatch):
    path = tmp_path / "artifact"
    path.write_bytes(b"abc")
    monkeypatch.setattr(release_upload, "HASH_READ_SIZE_BYTES", 2)

    assert sha256_file(path) == (
        "ba7816bf8f01cfea414140de5dae2223"
        "b00361a396177a9cb410ff61f20015ad"
    )


def test_pypi_release_digests_reads_project_release_map(monkeypatch):
    observed = {}

    def open_url(url, timeout):
        observed.update(url=url, timeout=timeout)
        return io.BytesIO(json.dumps({
            "releases": {
                "1.7.3": [{
                    "filename": SDIST,
                    "digests": {"sha256": VALID_SHA256.upper()},
                }],
            },
        }).encode("utf-8"))

    monkeypatch.setattr(release_upload, "urlopen", open_url)
    monkeypatch.setattr(release_upload, "token_hex", lambda length: "fresh-token")

    assert pypi_release_digests(
        "isovar",
        "1.7.3",
        json_base_url="https://packages.example/pypi/",
        request_timeout_seconds=7,
    ) == {SDIST: VALID_SHA256}
    assert observed == {
        "url": "https://packages.example/pypi/isovar/json?cache_bust=fresh-token",
        "timeout": 7,
    }


def test_pypi_release_digests_uses_fresh_url_for_each_poll(monkeypatch):
    urls = []
    tokens = iter(["first-token", "second-token"])

    def open_url(url, timeout):
        urls.append(url)
        return io.BytesIO(b'{"releases": {}}')

    monkeypatch.setattr(release_upload, "urlopen", open_url)
    monkeypatch.setattr(release_upload, "token_hex", lambda length: next(tokens))

    pypi_release_digests("isovar", "1.7.3")
    pypi_release_digests("isovar", "1.7.3")

    assert urls == [
        "https://pypi.org/pypi/isovar/json?cache_bust=first-token",
        "https://pypi.org/pypi/isovar/json?cache_bust=second-token",
    ]


def test_pypi_release_digests_treats_missing_project_as_empty(monkeypatch):
    def open_url(url, timeout):
        raise HTTPError(url, 404, "Not Found", {}, None)

    monkeypatch.setattr(release_upload, "urlopen", open_url)
    assert pypi_release_digests("isovar", "1.7.3") == {}


@pytest.mark.parametrize("payload", [
    None,
    [],
    {"releases": None},
    {"releases": []},
    {"releases": {"1.7.3": {}}},
])
def test_pypi_release_digests_rejects_invalid_release_metadata(
        monkeypatch,
        payload):
    monkeypatch.setattr(
        release_upload,
        "urlopen",
        lambda url, timeout: io.BytesIO(json.dumps(payload).encode("utf-8")),
    )

    with pytest.raises(ReleaseUploadError, match="Invalid PyPI release metadata"):
        pypi_release_digests("isovar", "1.7.3")


@pytest.mark.parametrize("item", [
    None,
    {},
    {"filename": SDIST},
    {"filename": SDIST, "digests": {}},
])
def test_pypi_release_digests_rejects_missing_metadata(monkeypatch, item):
    monkeypatch.setattr(
        release_upload,
        "urlopen",
        lambda url, timeout: io.BytesIO(json.dumps({
            "releases": {"1.7.3": [item]},
        }).encode("utf-8")),
    )

    with pytest.raises(ReleaseUploadError, match="Invalid PyPI artifact metadata"):
        pypi_release_digests("isovar", "1.7.3")


@pytest.mark.parametrize("digest", [None, "", "not-a-digest", "g" * 64])
def test_pypi_release_digests_rejects_invalid_sha256(monkeypatch, digest):
    monkeypatch.setattr(
        release_upload,
        "urlopen",
        lambda url, timeout: io.BytesIO(json.dumps({
            "releases": {"1.7.3": [{
                "filename": SDIST,
                "digests": {"sha256": digest},
            }]},
        }).encode("utf-8")),
    )

    with pytest.raises(ReleaseUploadError, match="SHA-256 digest"):
        pypi_release_digests("isovar", "1.7.3")


@pytest.mark.parametrize("filename", [None, ""])
def test_pypi_release_digests_rejects_invalid_filename(monkeypatch, filename):
    monkeypatch.setattr(
        release_upload,
        "urlopen",
        lambda url, timeout: io.BytesIO(json.dumps({
            "releases": {"1.7.3": [{
                "filename": filename,
                "digests": {"sha256": VALID_SHA256},
            }]},
        }).encode("utf-8")),
    )

    with pytest.raises(ReleaseUploadError, match="artifact filename"):
        pypi_release_digests("isovar", "1.7.3")


def test_pypi_release_digests_rejects_conflicting_duplicate(monkeypatch):
    monkeypatch.setattr(
        release_upload,
        "urlopen",
        lambda url, timeout: io.BytesIO(json.dumps({
            "releases": {"1.7.3": [
                {"filename": SDIST, "digests": {"sha256": VALID_SHA256}},
                {"filename": SDIST, "digests": {"sha256": OTHER_SHA256}},
            ]},
        }).encode("utf-8")),
    )

    with pytest.raises(ReleaseUploadError, match="Conflicting PyPI SHA-256"):
        pypi_release_digests("isovar", "1.7.3")


def test_upload_distribution_is_bounded_and_disables_progress(monkeypatch):
    observed = {}

    def run(command, check, timeout):
        observed.update(command=command, check=check, timeout=timeout)

    monkeypatch.setattr(subprocess, "run", run)
    upload_distribution(
        "dist/isovar-1.7.3.tar.gz",
        python_executable="release-python",
        timeout_seconds=17,
    )

    assert observed == {
        "command": [
            "release-python",
            "-m",
            "twine",
            "upload",
            "--disable-progress-bar",
            "dist/isovar-1.7.3.tar.gz",
        ],
        "check": True,
        "timeout": 17,
    }


def test_wait_for_release_file_observes_eventual_matching_digest():
    responses = iter([{}, {SDIST: VALID_SHA256}])

    assert wait_for_release_file(
        SDIST,
        VALID_SHA256,
        lambda: next(responses),
        timeout_seconds=1,
        poll_seconds=0,
    )


def test_wait_for_release_file_retries_transient_metadata_error():
    responses = iter([
        URLError("temporary metadata failure"),
        {SDIST: VALID_SHA256},
    ])

    def fetch_release_digests():
        response = next(responses)
        if isinstance(response, Exception):
            raise response
        return response

    assert wait_for_release_file(
        SDIST,
        VALID_SHA256,
        fetch_release_digests,
        timeout_seconds=1,
        poll_seconds=0,
    )


def test_wait_for_release_file_rejects_same_name_with_different_digest():
    with pytest.raises(ReleaseUploadError, match="SHA-256 mismatch"):
        wait_for_release_file(
            SDIST,
            VALID_SHA256,
            lambda: {SDIST: OTHER_SHA256},
            timeout_seconds=1,
            poll_seconds=0,
        )


@pytest.mark.parametrize(
    "upload_error",
    [
        pytest.param(
            subprocess.TimeoutExpired(["twine", "upload"], timeout=60),
            id="timeout",
        ),
        pytest.param(
            subprocess.CalledProcessError(1, ["twine", "upload"]),
            id="nonzero-exit",
        ),
    ],
)
def test_ambiguous_failure_is_reconciled_only_by_matching_digest(
        release_files,
        upload_error):
    local_digests = _local_digests(release_files)
    published = {}

    def upload_file(path):
        published[path.name] = local_digests[path.name]
        raise upload_error

    result = publish_release(
        release_files,
        expected_filenames=EXPECTED,
        fetch_release_digests=lambda: published,
        upload_file=upload_file,
        verify_timeout_seconds=0,
    )

    assert result == local_digests


def test_ambiguous_failure_with_different_server_bytes_fails_hard(release_files):
    published = {}

    def upload_file(path):
        published[path.name] = OTHER_SHA256
        raise subprocess.TimeoutExpired(["twine", "upload", path], timeout=60)

    with pytest.raises(ReleaseUploadError, match="SHA-256 mismatch"):
        publish_release(
            release_files,
            expected_filenames=EXPECTED,
            fetch_release_digests=lambda: published,
            upload_file=upload_file,
            verify_timeout_seconds=0,
        )


def test_partial_release_uploads_only_missing_file_and_retry_is_safe(
        release_files):
    local_digests = _local_digests(release_files)
    published = {WHEEL: local_digests[WHEEL]}
    uploaded = []

    def upload_file(path):
        uploaded.append(path.name)
        published[path.name] = local_digests[path.name]

    for _ in range(2):
        publish_release(
            release_files,
            expected_filenames=EXPECTED,
            fetch_release_digests=lambda: published,
            upload_file=upload_file,
            verify_timeout_seconds=0,
        )

    assert uploaded == [SDIST]


@pytest.mark.parametrize("conflicting_filename", [WHEEL, SDIST])
def test_existing_digest_mismatch_aborts_before_any_upload(
        release_files,
        conflicting_filename):
    published = {conflicting_filename: OTHER_SHA256}
    uploaded = []

    with pytest.raises(ReleaseUploadError, match="SHA-256 mismatch"):
        publish_release(
            release_files,
            expected_filenames=EXPECTED,
            fetch_release_digests=lambda: published,
            upload_file=lambda path: uploaded.append(path),
            verify_timeout_seconds=0,
        )

    assert uploaded == []


def test_failed_upload_without_server_artifact_fails_hard(release_files):
    def upload_file(path):
        raise subprocess.TimeoutExpired(["twine", "upload", path], timeout=60)

    with pytest.raises(ReleaseUploadError, match="did not publish"):
        publish_release(
            release_files,
            expected_filenames=EXPECTED,
            fetch_release_digests=dict,
            upload_file=upload_file,
            verify_timeout_seconds=0,
        )


def test_successful_upload_with_different_server_bytes_fails_hard(release_files):
    published = {}

    def upload_file(path):
        published[path.name] = OTHER_SHA256

    with pytest.raises(ReleaseUploadError, match="SHA-256 mismatch"):
        publish_release(
            release_files,
            expected_filenames=EXPECTED,
            fetch_release_digests=lambda: published,
            upload_file=upload_file,
            verify_timeout_seconds=0,
        )


def test_final_verification_rejects_digest_changed_after_exact_matches(
        release_files):
    local_digests = _local_digests(release_files)
    changed_digests = dict(local_digests)
    changed_digests[SDIST] = OTHER_SHA256
    responses = iter([
        local_digests,
        local_digests,
        changed_digests,
    ])

    with pytest.raises(ReleaseUploadError, match="SHA-256 mismatch"):
        publish_release(
            release_files,
            expected_filenames=EXPECTED,
            fetch_release_digests=lambda: next(responses),
            upload_file=lambda path: pytest.fail("unexpected upload"),
        )


def test_final_verification_rejects_missing_artifact(release_files):
    local_digests = _local_digests(release_files)
    responses = iter([
        local_digests,
        local_digests,
        {WHEEL: local_digests[WHEEL]},
    ])

    with pytest.raises(ReleaseUploadError, match="missing expected files"):
        publish_release(
            release_files,
            expected_filenames=EXPECTED,
            fetch_release_digests=lambda: next(responses),
            upload_file=lambda path: pytest.fail("unexpected upload"),
        )


def test_unexpected_distribution_set_is_rejected_before_hash_or_upload():
    with pytest.raises(ReleaseUploadError, match="do not match"):
        publish_release(
            ["dist/isovar-1.7.3.tar.gz"],
            expected_filenames=EXPECTED,
            fetch_release_digests=dict,
            upload_file=lambda path: None,
        )


def test_duplicate_distribution_basename_is_rejected_before_hash_or_upload():
    paths = list(PATHS) + [Path("other") / next(iter(EXPECTED))]

    with pytest.raises(ReleaseUploadError, match="do not match"):
        publish_release(
            paths,
            expected_filenames=EXPECTED,
            fetch_release_digests=dict,
            upload_file=lambda path: None,
        )


def test_missing_distribution_file_is_rejected_before_network_or_upload():
    fetched = []
    uploaded = []

    with pytest.raises(ReleaseUploadError, match="Could not hash"):
        publish_release(
            PATHS,
            expected_filenames=EXPECTED,
            fetch_release_digests=lambda: fetched.append(True),
            upload_file=lambda path: uploaded.append(path),
        )

    assert fetched == []
    assert uploaded == []


def test_main_passes_exact_artifact_contract_to_publisher(monkeypatch, capsys):
    observed = {}

    def publish_release_call(distribution_paths, **kwargs):
        observed.update(
            distribution_paths=tuple(distribution_paths),
            expected_filenames=kwargs["expected_filenames"],
            fetch_release_digests=kwargs["fetch_release_digests"],
        )
        return {filename: VALID_SHA256 for filename in EXPECTED}

    monkeypatch.setattr(release_upload, "publish_release", publish_release_call)

    assert main([
        "--project", "isovar",
        "--version", "1.7.3",
        *[str(path) for path in PATHS],
    ]) == 0
    assert observed["distribution_paths"] == tuple(str(path) for path in PATHS)
    assert observed["expected_filenames"] == EXPECTED
    assert callable(observed["fetch_release_digests"])
    assert "Verified PyPI release files" in capsys.readouterr().out
