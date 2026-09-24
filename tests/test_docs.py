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

import inspect
from pathlib import Path
import re

from isovar import IsovarResult, run_isovar
from isovar.default_parameters import DEFAULT_FILTER_FLAGS

ROOT = Path(__file__).resolve().parent.parent


def test_run_isovar_docstring_documents_every_parameter():
    doc = run_isovar.__doc__ or ""
    for name in inspect.signature(run_isovar).parameters:
        assert re.search(r"^\s*%s : " % name, doc, re.M), "%s undocumented" % name


def test_readme_filter_names_are_isovar_result_properties():
    readme = (ROOT / "README.md").read_text()
    thresholds = re.findall(r"""['"](?:min|max)_([a-z_]+)['"]""", readme)
    negated_flags = re.findall(r"`not_([a-z_]+)`", readme)
    assert thresholds and negated_flags
    for name in thresholds + negated_flags + list(DEFAULT_FILTER_FLAGS):
        assert hasattr(IsovarResult, name), name
    for name in DEFAULT_FILTER_FLAGS:
        assert "`%s`" % name in readme, name


REPOSITORY_URL = "https://github.com/openvax/isovar/blob/master/"


PAGES = [ROOT / "README.md", ROOT / "CHANGELOG.md", ROOT / "RELEASING.md", *(ROOT / "docs").glob("*.md")]


def heading_anchors(text):
    """GitHub's anchors for a Markdown page's headings, outside code blocks."""
    anchors, fenced = set(), False
    for line in text.splitlines():
        if line.startswith("```"):
            fenced = not fenced
        elif not fenced and re.match(r"#{1,6} ", line):
            slug = re.sub(r"[^\w\- ]", "", line.lstrip("#").strip().lower()).replace(" ", "-")
            anchor, n = slug, 0
            while anchor in anchors:
                n += 1
                anchor = "%s-%d" % (slug, n)
            anchors.add(anchor)
    return anchors


def test_markdown_links_to_this_repository_resolve():
    for page in PAGES:
        for target, fragment in re.findall(r"\]\(([^)#\s]*)(?:#([^)\s]*))?\)", page.read_text()):
            if target.startswith(REPOSITORY_URL):
                path = ROOT / target[len(REPOSITORY_URL):]
            elif "://" in target or target.startswith("mailto:"):
                continue
            else:
                path = page.parent / target if target else page
            assert path.exists(), "%s links to missing %s" % (page.name, target)
            if fragment and path.suffix == ".md":
                assert fragment in heading_anchors(path.read_text()), (
                    "%s links to missing section %s#%s" % (page.name, target, fragment))
