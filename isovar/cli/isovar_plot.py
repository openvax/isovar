"""Publication-ready figures for one mutation, using the standard RNA pipeline."""

import sys
from dataclasses import asdict
from collections.abc import Mapping
from pathlib import Path
from urllib.parse import urlsplit

from varcode.cli import variant_collection_from_args

from ..default_parameters import (
    PLOT_ALL_PROTEINS, PLOT_COMPARE_ASSEMBLY, PLOT_DPI, PLOT_MAX_ROWS, PLOT_OUTPUT_DIRECTORY, PLOT_VIEW, PLOT_VIEWS,
)
from ..logging import configure_cli_logging
from ..visualization import (
    _plot_imports, collect_visualization_data, save_variant_figures, timestamped_run_directory,
)
from .rna_args import alignment_file_from_args, read_collector_from_args
from .commands import parser_for_program
from .translation_args import (
    add_protein_selection_args,
    make_translation_arg_parser,
    protein_sequence_creator_kwargs_from_args,
)


parser = make_translation_arg_parser(description=(
    "Plot the protein context, read overlap and local transcript models around one mutation. "
    "Writes white-background SVG/PNG figures and evidence JSON into a new UTC-stamped directory."))
parser.add_argument("--view", choices=PLOT_VIEWS, default=PLOT_VIEW,
                    help="Figure panels (default %(default)s).")
parser.add_argument("--compare-assembly", action="store_true", default=PLOT_COMPARE_ASSEMBLY,
                    help="Compare assembly on/off with the same collected reads and all other settings unchanged.")
parser.add_argument("--output-dir", default=PLOT_OUTPUT_DIRECTORY,
                    help="Parent of new UTC-stamped run directories (default %(default)s).")
parser.add_argument("--max-rows", type=int, default=PLOT_MAX_ROWS,
                    help="Maximum displayed span/model rows; never limits analysis (default %(default)s).")
parser.add_argument("--dpi", type=int, default=PLOT_DPI, help="PNG resolution (default %(default)s); SVG is vector.")
parser.add_argument("--all-proteins", action="store_true", default=PLOT_ALL_PROTEINS,
                    help="Also write paginated protein/frame alternatives; removes only the protein result cap.")
parser.add_argument("--sample-label", help="Explicit sample/technology label for figures (default: BAM filename).")
add_protein_selection_args(parser)


def run(args=None, *, prog=None):
    command_parser = parser_for_program(parser, prog)
    args = command_parser.parse_args(sys.argv[1:] if args is None else args)
    if args.max_rows < 2 or args.dpi < 72:
        command_parser.error("--max-rows must be >= 2 and --dpi must be >= 72")
    try:
        _plot_imports()
    except ImportError as error:
        command_parser.error(str(error))
    configure_cli_logging()
    variants = variant_collection_from_args(args)
    if len(variants) != 1:
        command_parser.error("Select exactly one mutation (use --variant or a single-record variant file).")
    variant = next(iter(variants))
    collector = read_collector_from_args(args)
    with alignment_file_from_args(args) as alignment:
        evidence = collector.read_evidence_for_variant(variant, alignment)
    creator_kwargs = protein_sequence_creator_kwargs_from_args(args)
    if args.all_proteins:
        creator_kwargs["max_protein_sequences_per_variant"] = None
    data = collect_visualization_data(
        variant, evidence, creator_kwargs=creator_kwargs,
        compare_assembly=args.compare_assembly)
    label = args.sample_label or Path(urlsplit(args.bam).path).name
    data["provenance"] = dict(sample_label=label)
    data["inputs"] = dict(bam=args.bam, vcf=args.vcf, maf=args.maf, json_variants=args.json_variants)
    data["collection_settings"] = dict(
        min_mapping_quality=collector.min_mapping_quality,
        use_duplicate_reads=collector.use_duplicate_reads,
        use_secondary_alignments=collector.use_secondary_alignments,
        use_soft_clipped_bases=collector.use_soft_clipped_bases,
        use_reads_without_base_qualities=collector.use_reads_without_base_qualities,
        merge_overlapping_fragments=collector.merge_overlapping_fragments,
        infer_read_ends=collector.infer_read_ends,
        trim_adapters=collector.trim_adapters,
        trim_poly_a=collector.trim_poly_a,
        read_end_profile=(
            {group: asdict(profile) for group, profile in collector.read_end_profile.items()}
            if isinstance(collector.read_end_profile, Mapping) else
            asdict(collector.read_end_profile) if collector.read_end_profile is not None else None))
    directory = save_variant_figures(data, timestamped_run_directory(args.output_dir),
                                    view=args.view, max_rows=args.max_rows, dpi=args.dpi)
    if args.all_proteins:
        from ..protein_comparison import save_protein_comparison
        save_protein_comparison([dict(source=args.bam, label=label, visualization=data)],
                                directory / "protein-alternatives", dpi=args.dpi)
    print(directory)
