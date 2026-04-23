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
Transcript-anchored edit derived from a local Isovar assembly.
"""

from .value_object import ValueObject


class TranscriptAssemblyEdit(ValueObject):
    """
    Wrap a TranscriptEdit with the transcript it is relative to.

    TranscriptEdit coordinates only make sense for a particular transcript
    sequence, so PhaseGroup stores this transcript-anchored form rather than
    bare edits.
    """

    __slots__ = [
        "transcript_id",
        "transcript_name",
        "edit",
    ]

    @property
    def cdna_start(self):
        return self.edit.cdna_start

    @property
    def cdna_end(self):
        return self.edit.cdna_end

    @property
    def alt_bases(self):
        return self.edit.alt_bases

    @property
    def source_variant(self):
        return self.edit.source_variant
