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

from pathlib import Path

from isovar import run_isovar
from isovar.dataframe_builder import DataFrameBuilder


def test_run_isovar_docstring_matches_signature_and_return_type():
    doc = run_isovar.__doc__
    assert doc is not None
    assert "decompression_threads : int" in doc
    assert "decompress_threads" not in doc
    assert "list of IsovarResult" in doc
    assert "Generator of IsovarResult" not in doc


def test_readme_filter_flags_use_existing_isovar_result_property():
    readme = Path("README.md").read_text()
    assert "protein_sequence_matches_predicted_mutation_effect" in readme
    assert "not_protein_sequence_matches_predicted_mutation_effect" in readme
    assert "protein_sequence_matches_predicted_effect" not in readme


def test_dataframe_builder_exclude_docstring_describes_omitted_fields():
    doc = DataFrameBuilder.__init__.__doc__
    assert doc is not None
    assert "Field names from element_class which should be omitted" in doc
    assert "should be used as columns" not in doc
