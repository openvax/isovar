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
import stat
import subprocess
import sys
import tarfile


ROOT_FILES = (
    "LICENSE",
    "MANIFEST.in",
    "README.md",
    "deploy.sh",
    "lint.sh",
    "pyproject.toml",
    "release_upload.py",
    "requirements.txt",
    "test.sh",
)
RELEASE_TESTS = (
    "test_deploy_script.py",
    "test_gate_scripts.py",
    "test_release_upload.py",
)

SV_SUPPORT = (
    "examples/osteosarc_sv_validation.py", "docs/sv-rna.md", "tests/__init__.py",
    "tests/osteosarc_protein_helpers.py", "tests/reference_identity.py",
    "tests/data/osteosarc/protein_references.py", "tests/data/osteosarc/expansion/references.py",
    "tests/data/osteosarc/expansion/inventory.py",
)


def test_sdist_contains_release_tooling_required_by_release_tests(tmp_path):
    source_root = Path(__file__).resolve().parents[1]
    checkout = tmp_path / "checkout"
    checkout.mkdir()
    for filename in ROOT_FILES:
        shutil.copy2(source_root / filename, checkout / filename)
    shutil.copytree(source_root / "isovar", checkout / "isovar")
    packaged_tests = checkout / "tests"
    packaged_tests.mkdir()
    for filename in RELEASE_TESTS:
        shutil.copy2(source_root / "tests" / filename, packaged_tests / filename)

    for name in SV_SUPPORT:
        target = checkout / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_root / name, target)
    shutil.copytree(source_root / "tests/data/fusions", checkout / "tests/data/fusions",
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))

    dist_dir = checkout / "dist"
    dist_dir.mkdir()
    env = os.environ.copy()
    for name in tuple(env):
        if name.startswith("COV_CORE_"):
            del env[name]
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "from setuptools.build_meta import build_sdist; "
            "build_sdist('dist')",
        ],
        cwd=checkout,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr

    archives = tuple(dist_dir.glob("isovar-*.tar.gz"))
    assert len(archives) == 1
    with tarfile.open(archives[0], "r:gz") as archive:
        root = archive.getnames()[0].split("/", 1)[0]
        names = set(archive.getnames())
        release_scripts = ("deploy.sh", "lint.sh", "test.sh")
        for script_name in release_scripts:
            assert "%s/%s" % (root, script_name) in names
        assert "%s/release_upload.py" % root in names
        assert "%s/tests/test_deploy_script.py" % root in names
        assert "%s/tests/test_gate_scripts.py" % root in names
        assert "%s/tests/test_release_upload.py" % root in names
        for script_name in release_scripts:
            script_info = archive.getmember("%s/%s" % (root, script_name))
            assert script_info.mode & stat.S_IXUSR

        for name in SV_SUPPORT + ("tests/data/fusions/build_validation_events.py",
                                  "tests/data/fusions/validation-events.json"):
            assert "%s/%s" % (root, name) in names
        options = {"filter": "data"} if hasattr(tarfile, "data_filter") else {}
        archive.extractall(tmp_path / "unpacked", **options)
    # Exercise original-read reconstruction from the actual archive, with no
    # repository-relative imports, network, or external reference acquisition.
    smoke = subprocess.run([sys.executable, "-c", """
from pathlib import Path
from tempfile import TemporaryDirectory
from examples.osteosarc_sv_validation import DATA
from tests.data.fusions.audit_three_fusions import run_case
from isovar.sv_rna_comparison import compare_sv_rna_predictions
with TemporaryDirectory() as scratch:
    result, _ = run_case(DATA / 'three-fusions', 'PARD3B--CDKN2B-AS1-CDKN2B',
                         'ONT-T1', 'reverse', Path(scratch))
comparison = compare_sv_rna_predictions(result, dict(event_id=result['event_id'],
    reference_name=result['reference_name'], source='packaged-fixture-independent-oracle',
    predictions=[dict(prediction_id='known', amino_acids='MNISNIHISTQKKKKKKSRF', complete=True)]))
matched, = [h for h in comparison['hypotheses']
            if h['comparisons'][0]['status'] == 'same_amino_acid_sequence']
assert matched['full_interval_fragments'] == 4
assert matched['orf']['rna_support']['fragment_query_orientations']['original_query'] == 4
"""], cwd=tmp_path / "unpacked" / root, env=env, capture_output=True, text=True)
    assert smoke.returncode == 0, smoke.stderr
