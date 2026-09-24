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
Export the names and sequences of reads supporting each variant's alt allele.
"""

from .commands import run_dataframe_command
from .output_args import add_output_args
from .rna_args import make_rna_reads_arg_parser, variant_reads_dataframe_from_args

parser = add_output_args(
    make_rna_reads_arg_parser(description=__doc__),
    filename="isovar-variant-reads-result.csv",
    description="CSV of alt-allele read sequences")


def run(args=None, *, prog=None):
    run_dataframe_command(parser, variant_reads_dataframe_from_args, args, prog)
