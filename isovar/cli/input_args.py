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
Shared input-flag helpers defining Isovar's unified CLI input vocabulary (#656).

Canonical spellings:
- ``--bam``: indexed BAM/CRAM RNA alignments (defined in ``rna_args.py``)
- ``--input``: structured JSON / JSON.gz specification or upstream Isovar export
- ``--manifest``: multi-source or cohort manifest JSON
- ``--genome``: Varcode/PyEnsembl reference genome name or Ensembl release
- ``--genome-fasta``: local reference genome FASTA file (hidden alias: ``--reference-fasta``)
- ``--gtf``: local versioned transcript/CDS GTF annotation
- ``--sample-id``: sample identifier scoping read/event records
- ``--source``: read-set / library identifier (hidden alias: ``--source-id``)
- ``--library-cache``: reusable original-RNA index directory (hidden alias: ``--library-cache-dir``)
- ``--scratch-dir``: directory for temporary disk indexes (hidden alias: ``--sv-scratch-dir``)
"""

from .output_args import add_hidden_alias


def add_input_arg(
    parser,
    *,
    required=True,
    help="Input JSON (or .json.gz) file",
    action="store",
    type=str,
    metavar="FILE",
    dest="input",
    aliases=(),
    group=None,
):
    """Add canonical ``--input FILE`` option with optional hidden aliases."""
    target = group if group is not None else parser
    kwargs = {
        "dest": dest,
        "required": required,
        "action": action,
        "type": type,
        "metavar": metavar,
        "help": help,
    }
    action_obj = target.add_argument("--input", **kwargs)
    for alias in aliases:
        add_hidden_alias(target, alias, dest=dest, type=type, action=action)
    return action_obj


def add_manifest_arg(
    parser,
    *,
    required=True,
    help="Input manifest JSON file",
    type=str,
    metavar="FILE",
    dest="manifest",
    group=None,
):
    """Add canonical ``--manifest FILE`` option."""
    target = group if group is not None else parser
    return target.add_argument(
        "--manifest",
        dest=dest,
        required=required,
        type=type,
        metavar=metavar,
        help=help,
    )


def add_sample_id_arg(
    parser,
    *,
    required=False,
    help="Sample the reads came from",
    dest="sample_id",
    aliases=(),
    group=None,
):
    """Add canonical ``--sample-id`` option with optional hidden aliases."""
    target = group if group is not None else parser
    action_obj = target.add_argument(
        "--sample-id",
        dest=dest,
        required=required,
        help=help,
    )
    for alias in aliases:
        add_hidden_alias(target, alias, dest=dest)
    return action_obj


def add_source_arg(
    parser,
    *,
    required=False,
    help="Identity of the original RNA read set (default: --bam)",
    dest="source",
    metavar="SOURCE",
    aliases=(),
    group=None,
):
    """Add canonical ``--source`` option with optional hidden aliases (e.g. ``--source-id``)."""
    target = group if group is not None else parser
    if aliases and required:
        # When a required flag has a legacy alias, argparse cannot enforce
        # `required=True` on only the primary action without rejecting the alias.
        # A required mutually exclusive group accepts either spelling cleanly.
        req_group = target.add_mutually_exclusive_group(required=True)
        action_obj = req_group.add_argument(
            "--source", dest=dest, metavar=metavar, help=help
        )
        for alias in aliases:
            add_hidden_alias(req_group, alias, dest=dest)
        return action_obj
    action_obj = target.add_argument(
        "--source",
        dest=dest,
        required=required,
        metavar=metavar,
        help=help,
    )
    for alias in aliases:
        add_hidden_alias(target, alias, dest=dest)
    return action_obj


def add_genome_fasta_arg(
    parser,
    *,
    required=False,
    help="Reference genome FASTA file",
    metavar="FILE",
    dest="genome_fasta",
    aliases=(),
    group=None,
):
    """Add canonical ``--genome-fasta FILE`` with optional hidden aliases (e.g. ``--reference-fasta``)."""
    target = group if group is not None else parser
    if aliases and required:
        req_group = target.add_mutually_exclusive_group(required=True)
        action_obj = req_group.add_argument(
            "--genome-fasta", dest=dest, metavar=metavar, help=help
        )
        for alias in aliases:
            add_hidden_alias(req_group, alias, dest=dest)
        return action_obj
    action_obj = target.add_argument(
        "--genome-fasta",
        dest=dest,
        required=required,
        metavar=metavar,
        help=help,
    )
    for alias in aliases:
        add_hidden_alias(target, alias, dest=dest)
    return action_obj


def add_gtf_arg(
    parser,
    *,
    required=False,
    help="Local, versioned transcript/CDS annotation GTF",
    metavar="FILE",
    dest="gtf",
    group=None,
):
    """Add canonical ``--gtf FILE`` option."""
    target = group if group is not None else parser
    return target.add_argument(
        "--gtf",
        dest=dest,
        required=required,
        metavar=metavar,
        help=help,
    )


def add_library_cache_arg(
    parser,
    *,
    help="Reusable original-RNA index directory",
    metavar="DIR",
    dest="library_cache",
    aliases=(),
    group=None,
):
    """Add canonical ``--library-cache DIR`` with optional hidden aliases."""
    target = group if group is not None else parser
    action_obj = target.add_argument(
        "--library-cache",
        dest=dest,
        metavar=metavar,
        help=help,
    )
    for alias in aliases:
        add_hidden_alias(target, alias, dest=dest)
    return action_obj


def add_scratch_dir_arg(
    parser,
    *,
    help="Directory for temporary disk indexes",
    metavar="DIR",
    dest="scratch_dir",
    type=str,
    aliases=(),
    group=None,
):
    """Add canonical ``--scratch-dir DIR`` with optional hidden aliases."""
    target = group if group is not None else parser
    action_obj = target.add_argument(
        "--scratch-dir",
        dest=dest,
        type=type,
        metavar=metavar,
        help=help,
    )
    for alias in aliases:
        add_hidden_alias(target, alias, dest=dest, type=type)
    return action_obj
