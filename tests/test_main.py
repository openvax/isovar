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

from isovar import run_isovar, isovar_results_to_dataframe
from .common import eq_
from .testing_helpers import data_path

def test_isovar_main_to_dataframe():
    results = run_isovar(
        variants=data_path("data/b16.f10/b16.vcf"),
        alignment_file=data_path("data/b16.f10/b16.combined.sorted.bam"))
    df = isovar_results_to_dataframe(results)
    print(df)
    eq_(len(df), 4)
    # B16 test data has 2/4 variants with enough coverage
    # to translate protein sequences
    eq_(df["passes_all_filters"].sum(), 2)

