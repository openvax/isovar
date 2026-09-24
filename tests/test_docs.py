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

from isovar import run_isovar

ROOT = Path(__file__).resolve().parent.parent


def test_run_isovar_docstring_documents_every_parameter():
    doc = run_isovar.__doc__ or ""
    for name in inspect.signature(run_isovar).parameters:
        assert re.search(r"^\s*%s : " % name, doc, re.M), "%s undocumented" % name

