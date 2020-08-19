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

from __future__ import print_function, division, absolute_import

import sys

from ..logging import get_logger
from .reference_context_args import (
    make_reference_context_arg_parser,
    reference_contexts_dataframe_from_args
)
from .output_args import add_output_args, write_dataframe


logger = get_logger(__name__)

parser = make_reference_context_arg_parser()
parser = add_output_args(
    parser,
    filename="isovar-reference-contexts-result.csv")


def run(args=None):
    if args is None:
        args = sys.argv[1:]
    args = parser.parse_args(args)
    reference_contexts_df = reference_contexts_dataframe_from_args(args)
    logger.info(reference_contexts_df)
    write_dataframe(reference_contexts_df, args)
