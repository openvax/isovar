# Copyright (c) 2016-2019. Mount Sinai School of Medicine
#
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
Helper functions for extracting LocusRead objects from read alignment files at
specified locations.
"""

from __future__ import print_function, division, absolute_import

from .default_parameters import (
    MIN_READ_MAPPING_QUALITY,
    USE_DUPLICATE_READS,
    USE_SECONDARY_ALIGNMENTS,
)
from .locus_read import  LocusRead
from .logging import get_logger

logger = get_logger(__name__)


