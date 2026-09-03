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
    pypi_release_filenames,
    publish_release,
    upload_distribution,
    wait_for_release_file,
)


EXPECTED = expected_release_filenames("isovar", "1.7.3")
PATHS = tuple(Path("dist") / filename for filename in EXPECTED)


def test_expected_release_filenames_are_exact_wheel_and_sdist():
    assert EXPECTED == {
        "isovar-1.7.3-py3-none-any.whl",
        "isovar-1.7.3.tar.gz",
    }


def test_pypi_release_filenames_reads_project_release_map(monkeypatch):
    observed = {}

    def open_url(url, timeout):
        observed.update(url=url, timeout=timeout)
        return io.BytesIO(json.dumps({
            "releases": {
                "1.7.3": [{"filename": "isovar-1.7.3.tar.gz"}],
            },
        }).encode("utf-8"))

    monkeypatch.setattr(release_upload, "urlopen", open_url)
    monkeypatch.setattr(release_upload, "token_hex", lambda length: "fresh-token")

    assert pypi_release_filenames(
        "isovar",
        "1.7.3",
        json_base_url="https://packages.example/pypi/",
        request_timeout_seconds=7,
    ) == {"isovar-1.7.3.tar.gz"}
    assert observed == {
        "url": "https://packages.example/pypi/isovar/json?cache_bust=fresh-token",
        "timeout": 7,
    }


def test_pypi_release_filenames_uses_fresh_url_for_each_poll(monkeypatch):
    urls = []
    tokens = iter(["first-token", "second-token"])

    def open_url(url, timeout):
        urls.append(url)
        return io.BytesIO(b'{"releases": {}}')

    monkeypatch.setattr(release_upload, "urlopen", open_url)
    monkeypatch.setattr(release_upload, "token_hex", lambda length: next(tokens))

    pypi_release_filenames("isovar", "1.7.3")
    pypi_release_filenames("isovar", "1.7.3")

    assert urls == [
        "https://pypi.org/pypi/isovar/json?cache_bust=first-token",
        "https://pypi.org/pypi/isovar/json?cache_bust=second-token",
    ]


def test_pypi_release_filenames_treats_missing_project_as_empty(monkeypatch):
    def open_url(url, timeout):
        raise HTTPError(url, 404, "Not Found", {}, None)

    monkeypatch.setattr(release_upload, "urlopen", open_url)
    assert pypi_release_filenames("isovar", "1.7.3") == frozenset()


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


def test_wait_for_release_file_observes_eventual_metadata():
    responses = iter([set(), {"isovar-1.7.3.tar.gz"}])

    assert wait_for_release_file(
        "isovar-1.7.3.tar.gz",
        lambda: next(responses),
        timeout_seconds=1,
        poll_seconds=0,
    )


def test_wait_for_release_file_retries_transient_metadata_error():
    responses = iter([
        URLError("temporary metadata failure"),
        {"isovar-1.7.3.tar.gz"},
    ])

    def fetch_release_filenames():
        response = next(responses)
        if isinstance(response, Exception):
            raise response
        return response

    assert wait_for_release_file(
        "isovar-1.7.3.tar.gz",
        fetch_release_filenames,
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
def test_ambiguous_failure_after_server_acceptance_is_reconciled(upload_error):
    published = set()

    def fetch_release_filenames():
        return published

    def upload_file(path):
        published.add(path.name)
        raise upload_error

    result = publish_release(
        PATHS,
        expected_filenames=EXPECTED,
        fetch_release_filenames=fetch_release_filenames,
        upload_file=upload_file,
        verify_timeout_seconds=0,
    )

    assert EXPECTED <= result


def test_partial_release_uploads_only_missing_file_and_retry_is_safe():
    wheel = "isovar-1.7.3-py3-none-any.whl"
    sdist = "isovar-1.7.3.tar.gz"
    published = {wheel}
    uploaded = []

    def fetch_release_filenames():
        return published

    def upload_file(path):
        uploaded.append(path.name)
        published.add(path.name)

    for _ in range(2):
        publish_release(
            PATHS,
            expected_filenames=EXPECTED,
            fetch_release_filenames=fetch_release_filenames,
            upload_file=upload_file,
            verify_timeout_seconds=0,
        )

    assert uploaded == [sdist]


def test_failed_upload_without_server_artifact_fails_hard():
    def upload_file(path):
        raise subprocess.TimeoutExpired(["twine", "upload", path], timeout=60)

    with pytest.raises(ReleaseUploadError, match="did not publish"):
        publish_release(
            PATHS,
            expected_filenames=EXPECTED,
            fetch_release_filenames=frozenset,
            upload_file=upload_file,
            verify_timeout_seconds=0,
        )


def test_unexpected_distribution_set_is_rejected_before_upload():
    with pytest.raises(ReleaseUploadError, match="do not match"):
        publish_release(
            ["dist/isovar-1.7.3.tar.gz"],
            expected_filenames=EXPECTED,
            fetch_release_filenames=frozenset,
            upload_file=lambda path: None,
        )


def test_duplicate_distribution_basename_is_rejected_before_upload():
    paths = list(PATHS) + [Path("other") / next(iter(EXPECTED))]

    with pytest.raises(ReleaseUploadError, match="do not match"):
        publish_release(
            paths,
            expected_filenames=EXPECTED,
            fetch_release_filenames=frozenset,
            upload_file=lambda path: None,
        )


def test_main_passes_exact_artifact_contract_to_publisher(monkeypatch, capsys):
    observed = {}

    def publish_release_call(distribution_paths, **kwargs):
        observed.update(
            distribution_paths=tuple(distribution_paths),
            expected_filenames=kwargs["expected_filenames"],
        )
        return EXPECTED

    monkeypatch.setattr(release_upload, "publish_release", publish_release_call)

    assert main([
        "--project", "isovar",
        "--version", "1.7.3",
        *[str(path) for path in PATHS],
    ]) == 0
    assert observed == {
        "distribution_paths": tuple(str(path) for path in PATHS),
        "expected_filenames": EXPECTED,
    }
    assert "Verified PyPI release files" in capsys.readouterr().out
