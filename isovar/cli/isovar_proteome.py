"""Convert recorded RNA hypotheses to a source-qualified ORFdata Proteome."""

import argparse
import gzip
import json
from pathlib import Path

from orfdata import ReferenceContext

from ..rna_proteome import from_rna_result
from .validation import check_output_path


def _read(path):
    path = Path(path)
    with (gzip.open if path.suffix == '.gz' else open)(path, 'rt') as handle:
        return json.load(handle)


def run(args=None, *, prog=None):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    parser.add_argument('--input', action='append', required=True,
                        help='Recorded Isovar RNA JSON/JSON.gz result; repeat to combine sources')
    parser.add_argument('--reference-context', required=True,
                        help='JSON object of ORFdata ReferenceContext fields; must agree with recorded inputs')
    parser.add_argument('--output', required=True, help='ORFdata Proteome JSON or JSON.gz output')
    parser.add_argument('--fasta', help='Also write all retained protein hypotheses as FASTA')
    options = parser.parse_args(args)
    try:
        inputs = {Path(path).resolve() for path in options.input + [options.reference_context]}
        outputs = [Path(options.output).resolve()]
        if options.fasta:
            outputs.append(Path(options.fasta).resolve())
        if len(set(outputs)) != len(outputs) or inputs.intersection(outputs):
            raise ValueError('Output paths must be distinct from each other and the input/reference paths')
        for output in outputs:
            check_output_path(output)
        reference = ReferenceContext(**_read(options.reference_context))
        proteome = None
        for path in options.input:
            converted = from_rna_result(_read(path), reference=reference)
            proteome = converted if proteome is None else proteome.merge(converted)
        proteome.save(options.output)
        if options.fasta:
            proteome.write_fasta(options.fasta)
    except (ValueError, TypeError, OSError) as error:
        parser.error(str(error))
