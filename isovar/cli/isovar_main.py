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
Primary Isovar command, used to collect information about variants,
the RNA reads which overlap and protein sequences which can be constructed
from reads that support the variant.
"""
from __future__ import print_function, division, absolute_import

import sys



from ..logging import get_logger
from ..dataframe_helpers import isovar_results_to_dataframe

from .main_args import run_isovar_from_parsed_args, make_isovar_arg_parser

from .output_args import add_output_args, write_dataframe

logger = get_logger(__name__)

def run(args=None):
    if args is None:
        args = sys.argv[1:]
    parser = make_isovar_arg_parser()
    parser = add_output_args(
        parser,
        filename="isovar-results.csv")
    args = parser.parse_args(args)
    logger.info(args)
    isovar_results = run_isovar_from_parsed_args(args)
    df = isovar_results_to_dataframe(isovar_results)
    logger.info(df)
    write_dataframe(df, args)
