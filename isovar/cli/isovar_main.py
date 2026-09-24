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
Run the complete pipeline: collect RNA reads overlapping each variant, assemble
and translate mutant protein sequences, and evaluate result filters.
"""

from ..dataframe_helpers import isovar_results_to_dataframe
from .commands import run_dataframe_command
from .main_args import run_isovar_from_parsed_args, make_isovar_arg_parser
from .output_args import add_output_args


def isovar_results_dataframe_from_args(args):
    return isovar_results_to_dataframe(run_isovar_from_parsed_args(args))


def run(args=None, *, prog=None):
    parser = add_output_args(
        make_isovar_arg_parser(prog=prog, description=__doc__),
        filename="isovar-results.csv",
        description="CSV with one row per variant")
    run_dataframe_command(parser, isovar_results_dataframe_from_args, args)
