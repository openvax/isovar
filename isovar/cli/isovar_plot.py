"""Publication-ready figures for one mutation, using the standard RNA pipeline."""

import sys

from varcode.cli import variant_collection_from_args

from ..default_parameters import (
    PLOT_COMPARE_ASSEMBLY, PLOT_DPI, PLOT_MAX_ROWS, PLOT_OUTPUT_DIRECTORY, PLOT_VIEW,
)
from ..visualization import (
    _plot_imports, collect_visualization_data, save_variant_figures, timestamped_run_directory,
)
from .rna_args import alignment_file_from_args, read_collector_from_args
from .translation_args import make_translation_arg_parser, protein_sequence_creator_kwargs_from_args


parser = make_translation_arg_parser(description=(
    "Plot the protein context, read overlap and local transcript models around one mutation. "
    "Writes white-background SVG/PNG figures and evidence JSON into a new UTC-stamped directory."))
parser.add_argument("--view", choices=("all", "protein", "assembly", "transcripts"), default=PLOT_VIEW,
                    help="Figure panels (default %(default)s).")
parser.add_argument("--compare-assembly", action="store_true", default=PLOT_COMPARE_ASSEMBLY,
                    help="Compare assembly on/off with the same collected reads and all other settings unchanged.")
parser.add_argument("--output-dir", default=PLOT_OUTPUT_DIRECTORY,
                    help="Parent of new UTC-stamped run directories (default %(default)s).")
parser.add_argument("--max-rows", type=int, default=PLOT_MAX_ROWS,
                    help="Maximum displayed span/model rows; never limits analysis (default %(default)s).")
parser.add_argument("--dpi", type=int, default=PLOT_DPI, help="PNG resolution (default %(default)s); SVG is vector.")


def run(args=None):
    args = parser.parse_args(sys.argv[1:] if args is None else args)
    if args.max_rows < 2 or args.dpi < 72:
        parser.error("--max-rows must be >= 2 and --dpi must be >= 72")
    try:
        _plot_imports()
    except ImportError as error:
        parser.error(str(error))
    variants = variant_collection_from_args(args)
    if len(variants) != 1:
        parser.error("Select exactly one mutation (use --variant or a single-record variant file).")
    variant = next(iter(variants))
    collector = read_collector_from_args(args)
    with alignment_file_from_args(args) as alignment:
        evidence = collector.read_evidence_for_variant(variant, alignment)
    data = collect_visualization_data(
        variant, evidence, creator_kwargs=protein_sequence_creator_kwargs_from_args(args),
        compare_assembly=args.compare_assembly)
    data["inputs"] = dict(bam=args.bam, vcf=args.vcf, maf=args.maf, json_variants=args.json_variants)
    data["collection_settings"] = dict(
        min_mapping_quality=collector.min_mapping_quality,
        use_duplicate_reads=collector.use_duplicate_reads,
        use_secondary_alignments=collector.use_secondary_alignments,
        use_soft_clipped_bases=collector.use_soft_clipped_bases,
        merge_overlapping_fragments=collector.merge_overlapping_fragments)
    directory = save_variant_figures(data, timestamped_run_directory(args.output_dir),
                                    view=args.view, max_rows=args.max_rows, dpi=args.dpi)
    print(directory)


if __name__ == "__main__":
    run()
