# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from isovar.common import normalize_base0_range_indices
from .common import eq_

def normalize_base0_range_indices_valid_range():
    eq_(normalize_base0_range_indices(3, 7, 10),
        (3, 7))

def normalize_base0_range_indices_negative_indices():
    eq_(normalize_base0_range_indices(-3, -2, 10),
        (7, 8))

def normalize_base0_range_indices_start_None():
    eq_(normalize_base0_range_indices(None, -2, 10),
        (0, 8))
def normalize_base0_range_indices_end_None():
    eq_(normalize_base0_range_indices(-3, None, 10),
        (7, 10))
