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

from os.path import dirname, isdir, isfile

from pysam import AlignmentFile
from varcode.reference import infer_genome


def check_parsed_args(parser, args):
    """
    Call parser.error (exit status 2) for missing input files, a missing
    or unindexed alignment file, an unknown --genome, no variants, or an
    --output path whose directory doesn't exist.
    """
    for flag, paths in (
            ("--vcf", args.vcf),
            ("--maf", args.maf),
            ("--json-variants", args.json_variants)):
        for path in paths:
            if not isfile(path):
                parser.error("%s file not found: %s" % (flag, path))

    if not (args.vcf or args.maf or args.variant or args.json_variants):
        parser.error(
            "no variants given; use --vcf, --maf, --variant or --json-variants")

    if args.variant and not args.genome:
        parser.error("--genome is required when using --variant")

    if args.genome:
        try:
            infer_genome(args.genome)
        except ValueError as e:
            parser.error("--genome %s: %s" % (args.genome, e))

    bam_path = getattr(args, "bam", None)
    if bam_path is not None:
        if not isfile(bam_path):
            parser.error("--bam file not found: %s" % bam_path)
        try:
            with AlignmentFile(bam_path) as alignment_file:
                is_indexed = (
                    (alignment_file.is_bam or alignment_file.is_cram)
                    and alignment_file.has_index())
        except (ValueError, OSError) as e:
            parser.error("--bam %s: %s" % (bam_path, e))
        if not is_indexed:
            parser.error(
                "--bam must be a coordinate-sorted, indexed BAM or CRAM file "
                "(try 'samtools sort' and 'samtools index'): %s" % bam_path)

    output_dir = dirname(args.output)
    if output_dir and not isdir(output_dir):
        parser.error("--output directory does not exist: %s" % output_dir)
