"""Recompute the bounded DDX17 audit from pinned original reads and placements."""

from collections import Counter, defaultdict
import gzip
import hashlib
import json
from pathlib import Path

import pysam

from isovar.dna import reverse_complement_dna
from isovar.genetic_code import translate_cdna


HERE = Path(__file__).parent / 'ddx17_audit'
OFFSET = 38470000
SITE = 38489857


def build(directory=HERE):
    """Keep alternative placements and distinguish local uniqueness from phasing."""
    manifest = json.loads((directory / 'manifest.json').read_text())
    for name, digest in manifest['files'].items():
        if hashlib.sha256((directory / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Input digest mismatch: {name}')
    alleles = json.loads((directory / 'allele-fragments.json').read_text())
    allele_by_name = defaultdict(set)
    for allele, names in alleles.items():
        for name in names:
            allele_by_name[name].add(allele)
    models = json.loads((directory / 'models.json').read_text())
    reference = ''.join((directory / 'reference.fa').read_text().splitlines()[1:])
    by_fragment = defaultdict(list)
    with pysam.AlignmentFile(directory / 'expanded.star.bam') as bam:
        for read in bam:
            by_fragment[read.query_name].append(read)
    junctions = defaultdict(lambda: defaultdict(set))
    counts = Counter()
    coverage = [0] * len(reference)
    alt_coverage = [0] * len(reference)
    alt_blocks = []
    for name, reads in by_fragment.items():
        labels = allele_by_name.get(name, {'not_in_original_allele_set'})
        unique = reads[0].get_tag('NH') == 1
        for allele in sorted(labels):
            counts[f'{allele}_mapped_fragments'] += 1
        if unique:
            for allele in sorted(labels):
                counts[f'{allele}_locally_unique_fragments'] += 1
            positions = {position for read in reads for position in read.get_reference_positions()}
            for position in positions:
                coverage[position] += 1
                if 'alt' in labels:
                    alt_coverage[position] += 1
            if 'alt' in labels:
                alt_blocks.extend((left + OFFSET, right + OFFSET) for read in reads
                                  for left, right in read.get_blocks())
        for read in reads:
            position = read.reference_start + OFFSET
            for index, (operation, length) in enumerate(read.cigartuples):
                if operation == 3:
                    row = junctions[position, position + length]
                    row['fragments'].add(name)
                    for allele in sorted(labels):
                        row[f'{allele}_fragments'].add(name)
                    if unique:
                        row['locally_unique_fragments'].add(name)
                    row['anchor_pairs'].add((read.cigartuples[index - 1][1],
                                             read.cigartuples[index + 1][1]))
                    row['cigars'].add(read.cigarstring)
                if operation in (0, 2, 3, 7, 8):
                    position += length
    rows = []
    for (left, right), support in sorted(junctions.items()):
        row = dict(intron_start_0based=left, intron_end_0based=right,
                   genomic_motif=reference[left - OFFSET:left - OFFSET + 2] + '-' +
                   reference[right - OFFSET - 2:right - OFFSET],
                   **{key: sorted(value) for key, value in support.items()})
        row['anchor_pairs'] = [list(pair) for pair in row['anchor_pairs']]
        rows.append(row)
    # These are explicitly reference-path predictions, not ALT-phased proteins.
    hypotheses = []
    for model in models['models']:
        if not model['complete']:
            continue
        last = model['upstream_exon_last_offset']
        partial_codon = (last - model['cds_start'] + 1) % 3
        prefix = model['sequence'][last - partial_codon + 1:last + 1] if partial_codon else ''
        for endpoint, label in ((38492055, 'retained intron'),
                                (38490436, 'extended intronic acceptor'),
                                (38491154, 'second intronic acceptor')):
            sequence = prefix + reverse_complement_dna(reference[38488115 - OFFSET:endpoint - OFFSET])
            amino_acids, stopped = translate_cdna(sequence)
            assert stopped
            stop_high = endpoint - (len(amino_acids) * 3 - partial_codon)
            hypotheses.append(dict(transcript=model['id'], path=label,
                                   first_intronic_base_1based=endpoint,
                                   reference_frame_extension=amino_acids + '*',
                                   first_stop_1based=[stop_high - 2, stop_high],
                                   mutation_downstream_of_stop=SITE < stop_high - 2,
                                   status='reference-path prediction; not an ALT-phased protein'))
    intron = coverage[38488115 - OFFSET:38492055 - OFFSET]
    with gzip.open(directory / 'selected-original-alignments.json.gz', 'rt') as stream:
        original = json.load(stream)
    result = dict(
        call_id='chr22:38489856:GA>G', gene='DDX17', strand='-',
        coding_reference_contexts=models['coding_reference_contexts'],
        original_allele_fragments={key: len(value) for key, value in alleles.items()},
        fragments_in_multiple_allele_groups=sum(len(labels) > 1 for labels in allele_by_name.values()),
        original_selected_mapped_read_ends=len(original),
        original_selected_spliced_read_ends=sum('N' in row['CIGAR_LONG'] for row in original),
        counts=dict(counts), junctions=rows,
        locally_unique_alt_block_extent_1based=[min(a for a, b in alt_blocks) + 1,
                                               max(b for a, b in alt_blocks)],
        intron_min_locally_unique_fragment_depth=min(intron),
        intron_uncovered_bases=sum(value == 0 for value in intron),
        coverage_start_0based=38487800, coverage_bin_size=25,
        coverage_max_per_bin=[max(coverage[start - OFFSET:start - OFFSET + 25])
                              for start in range(38487800, 38492300, 25)],
        alt_coverage_max_per_bin=[max(alt_coverage[start - OFFSET:start - OFFSET + 25])
                                  for start in range(38487800, 38492300, 25)],
        reference_frame_hypotheses=hypotheses,
        conclusion='RNA supports intron-containing paths and alternative splicing; '
                   'no mutation-phased coding path is established. The original protein '
                   'pipeline stopped at its annotated coding-reference gate before RNA assembly.')
    (directory / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    result = build()
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('junctions', 'coverage_max_per_bin',
                                     'alt_coverage_max_per_bin')}, indent=2))
