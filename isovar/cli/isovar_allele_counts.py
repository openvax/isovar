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
Count reads and fragments supporting the ref, alt and other alleles at each
variant locus.
"""

from .commands import run_dataframe_command
from .output_args import add_output_args
from .rna_args import make_rna_reads_arg_parser, allele_counts_dataframe_from_args
from .validation import CommandInputError

parser = add_output_args(
    make_rna_reads_arg_parser(description=__doc__),
    filename="isovar-allele-counts-result.csv",
    description="CSV of read and fragment counts")
_cells = parser.add_argument_group("Cell/UMI labels (single-cell data)")
_cells.add_argument(
    "--cell-umi-labels", action="store_true",
    help="Add num_{ref,alt,other}_{cells,cell_umi_labels,unlabeled_segments} and "
         "num_cells_with_ref_and_alt columns from CB/UB tags; needs --sample-id")
_cells.add_argument("--sample-id", help="Sample the reads came from; scopes the labels")
_cells.add_argument("--source", help="Identity of the read set (default: --bam)")


def _counts(args):
    if args.cell_umi_labels and not args.sample_id:
        raise CommandInputError("--cell-umi-labels needs --sample-id")
    return allele_counts_dataframe_from_args(args)


def run(args=None, *, prog=None):
    run_dataframe_command(parser, _counts, args, prog)
