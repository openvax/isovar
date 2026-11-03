"""Convert recorded RNA hypotheses to a source-qualified ORFdata Proteome.

Inputs are the JSON written by `isovar protein-hypotheses`
(isovar.protein_hypotheses.v2), `isovar sv-rna` (isovar.sv_rna_orfs.v4 or
isovar.sv_rna_candidates.v5), `isovar fusion` (isovar.fusion_rna.v3) or the
variant ORF export (isovar.variant_rna_orfs.v1, library API only). The
--reference-context build, and the Ensembl release when the input records it,
must match what produced the input.
"""

import argparse
from pathlib import Path

from orfdata import ReferenceContext

from .output_args import add_log_level_arg
from ..json_files import read_json
from ..logging import configure_cli_logging, get_logger
from ..rna_proteome import from_rna_result
from .validation import check_output_path, input_error

logger = get_logger(__name__)


def run(args=None, *, prog=None):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', action='append', required=True,
                        help='Recorded Isovar RNA JSON/JSON.gz result; repeat to combine sources')
    parser.add_argument('--reference-context', required=True,
                        help='JSON object of ORFdata ReferenceContext fields; must agree with recorded inputs')
    parser.add_argument('--output', required=True, help='ORFdata Proteome JSON or JSON.gz output')
    parser.add_argument('--fasta', help='Also write all retained protein hypotheses as FASTA')
    add_log_level_arg(parser)
    options = parser.parse_args(args)
    configure_cli_logging(options.log_level)
    try:
        inputs = {Path(path).resolve() for path in options.input + [options.reference_context]}
        outputs = [Path(options.output).resolve()]
        if options.fasta:
            outputs.append(Path(options.fasta).resolve())
        if len(set(outputs)) != len(outputs) or inputs.intersection(outputs):
            raise ValueError('Output paths must be distinct from each other and the input/reference paths')
        for output in outputs:
            check_output_path(output)
        reference = ReferenceContext(**read_json(options.reference_context))
        proteome = None
        for path in options.input:
            converted = from_rna_result(read_json(path), reference=reference)
            try:
                proteome = converted if proteome is None else proteome.merge(converted)
            except ValueError as error:
                if 'Conflicting producer record identity' not in str(error):
                    raise
                raise ValueError(f'{path} repeats a sample and source already merged; give each '
                                 'input a distinct --source when it is produced') from error
        proteome.save(options.output)
        if options.fasta:
            proteome.write_fasta(options.fasta)
        logger.info("Wrote %s%s", options.output, " and " + str(options.fasta) if options.fasta else "")
    except (ValueError, TypeError, OSError) as error:
        input_error(parser, error)
