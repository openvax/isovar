"""Reconstruct original ALT/REF RNA paths before assessing supplied CDS models."""

import argparse
from contextlib import ExitStack
from pathlib import Path

import pysam

from ..contig_resolution import variant_reference_name
from ..rna_library import prepare_rna_library
from ..sv_rna_audit import read_json, write_json
from ..variant_rna_paths import reconstruct_variant_rna_paths, variant_rna_path_input_from_dict
from .rna_args import add_rna_args, read_collector_from_args


def make_parser(prog='isovar variant-rna-paths'):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    add_rna_args(parser)
    parser.set_defaults(min_mapping_quality=0)
    parser.add_argument('--input', required=True, help='Explicit literal variant, assembly and versioned models JSON')
    parser.add_argument('--source', required=True, help='Original RNA library identity; keep cultures/libraries separate')
    parser.add_argument('--output', required=True, help='New JSON or JSON.gz result path')
    parser.add_argument('--library-cache', help='Reusable original-RNA index directory')
    parser.add_argument('--scratch-dir')
    parser.add_argument('--reference-fasta', help='Reference for decoding CRAM input')
    for name, default in [('max-discovery-reads', 10000), ('max-seed-segments', 10000), ('max-paths', 256),
                          ('min-overlap', 30), ('min-anchor-bases', 12), ('flank', 1000)]:
        parser.add_argument('--' + name, type=int, default=default, help='Default: %(default)s')
    return parser


def run(args=None, *, prog=None):
    parser = make_parser(prog or 'isovar variant-rna-paths')
    options = parser.parse_args(args)
    try:
        output = Path(options.output)
        if output.exists():
            raise ValueError('Output exists; choose a new result path')
        supplied = variant_rna_path_input_from_dict(read_json(options.input))
        collector = read_collector_from_args(options)
        parameters = {name: getattr(options, name) for name in
                      ('max_discovery_reads', 'max_seed_segments', 'max_paths', 'min_overlap', 'min_anchor_bases', 'flank')}
        with ExitStack() as contexts:
            bam = contexts.enter_context(pysam.AlignmentFile(options.bam, reference_filename=options.reference_fasta,
                                         threads=options.num_rna_decompression_threads))
            library = None
            if options.library_cache:
                prepared = prepare_rna_library(options.bam, options.library_cache, read_collector=collector,
                                               reference_name=variant_reference_name(supplied['variant']))
                library = contexts.enter_context(prepared.open())
            result = reconstruct_variant_rna_paths(bam, **supplied, **parameters, source=options.source,
                read_collector=collector, library_index=library, scratch_dir=options.scratch_dir)
        output.parent.mkdir(parents=True, exist_ok=True)
        write_json(output, result)
    except (OSError, ValueError, TypeError) as error:
        parser.error(str(error))
