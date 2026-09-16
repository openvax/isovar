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
Check parsed commandline arguments before any reads are loaded, so that
common mistakes end with a one-line argparse error instead of a traceback
(or a traceback after the whole analysis has run).
"""

import re
from os.path import dirname, expanduser, isdir, isfile

from pysam import AlignmentFile
from varcode.reference import infer_genome

# htslib opens http/https/s3/gs/ftp URLs directly and varcode's load_vcf
# downloads HTTP(S) URLs to a temporary file, so a remote path is legal
# input even though it isn't on the local filesystem. Match a scheme rather
# than a fixed list so a Windows drive letter ("C:\...") isn't mistaken for
# one.
URL_SCHEME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]+://")


def is_url(path):
    """
    True if the path names a remote resource rather than a local file.
    """
    return bool(URL_SCHEME_RE.match(path))


def check_parsed_args(parser, args):
    """
    Call parser.error (exit status 2) for missing input files, a missing
    or unindexed alignment file, an unknown --genome, no variants, or an
    --output path whose directory doesn't exist.

    Remote (URL) inputs are left alone: they can't be checked without
    network access, and whichever loader consumes them reports its own
    error.
    """
    for flag, paths in (
            ("--vcf", args.vcf),
            ("--maf", args.maf),
            ("--json-variants", args.json_variants)):
        for path in paths:
            if not is_url(path) and not isfile(path):
                parser.error("%s file not found: %s" % (flag, path))

    if not (args.vcf or args.maf or args.variant or args.json_variants):
        parser.error(
            "no variants given; use --vcf, --maf, --variant or --json-variants")

    if args.variant and not args.genome:
        parser.error("--genome is required when using --variant")

    # --genome is ignored for MAF files (each row names its own reference)
    # and for serialized VariantCollections, so only check it when it will
    # actually be used to interpret coordinates.
    if args.genome and (args.vcf or args.variant):
        try:
            infer_genome(args.genome)
        except ValueError as e:
            parser.error("--genome %s: %s" % (args.genome, e))

    bam_path = getattr(args, "bam", None)
    if bam_path is not None and not is_url(bam_path):
        if not isfile(bam_path):
            parser.error("--bam file not found: %s" % bam_path)
        try:
            with AlignmentFile(bam_path) as alignment_file:
                is_indexed = (
                    (alignment_file.is_bam or alignment_file.is_cram)
                    and alignment_file.has_index())
        except (ValueError, OSError) as e:
            parser.error("--bam %s: %s" % (bam_path, e))
        else:
            if not is_indexed:
                parser.error(
                    "--bam must be a coordinate-sorted, indexed BAM or CRAM file "
                    "(try 'samtools sort' and 'samtools index'): %s" % bam_path)

    # pandas expands "~" when writing, so expand it here too before asking
    # whether the destination directory exists
    output_dir = dirname(expanduser(args.output))
    if output_dir and not isdir(output_dir):
        parser.error("--output directory does not exist: %s" % output_dir)
