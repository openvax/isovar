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

"""
Cooperative initializer base for adapters backed by IsovarResult iterables.
"""


class IsovarResultProvider(object):
    """
    Terminal cooperative base for providers initialized from IsovarResult
    collections.

    Subclasses should forward ``isovar_results`` with ``super()`` so ordinary
    multiple inheritance composes provider indexes safely.
    """

    def __init__(self, isovar_results, *args, **kwargs):
        super().__init__(*args, **kwargs)
