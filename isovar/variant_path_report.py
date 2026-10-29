"""Cohort presentation of conditional literal-allele RNA path evidence."""

from html import escape
from pathlib import Path

from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator
from reportlab.platypus import Image, KeepTogether, Paragraph


ALLELE_COLOR = '#d26928'
RNA_COLOR = '#128681'
FRAME_COLOR = '#687780'


def path_figure(candidate: dict, directory: Path) -> Path:
    """Plot original-segment coverage, genomic blocks and conditional frames."""
    figure = Figure(figsize=(7.2, 2.8), layout='constrained')
    mapping, frames, coverage = figure.subplots(3, 1, sharex=True, gridspec_kw={'height_ratios': [1, 1, 2]})
    blocks = candidate['mapping']
    for index, (start, end, contig, low, high, strand) in enumerate(blocks):
        mapping.broken_barh([(start, end - start)], (.2, .3), facecolors=RNA_COLOR)
        if index in (0, len(blocks) - 1):
            mapping.text(start if index == 0 else end, .7, f'{contig}:{low + 1}-{high} ({strand})',
                         ha='left' if index == 0 else 'right', fontsize=6)
    mapping.set_ylim(0, 1)
    mapping.axis('off')
    spans = sorted({(row['translation_start'], row['translation_end']) for row in candidate['coding']})
    for index, (start, end) in enumerate(spans):
        frames.broken_barh([(start, end - start)], (index, .6), facecolors=FRAME_COLOR)
    frames.set_ylim(-.2, max(1, len(spans)))
    frames.set_yticks([])
    frames.set_ylabel('CDS frames' if spans else 'No CDS frame', fontsize=7)
    frames.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
    coverage.fill_between(range(len(candidate['sequence'])), candidate['coverage'], color=RNA_COLOR, alpha=.85)
    for start, end in {tuple(row['query_interval']) for row in candidate['focal_alleles'] if row['state'] == 'ALT'}:
        for axis in (mapping, frames, coverage):
            axis.axvline(start, color=ALLELE_COLOR, linewidth=1)
            if end > start:
                axis.axvspan(start, end, color=ALLELE_COLOR, alpha=.18)
    coverage.set_ylabel('Original segments', fontsize=7)
    coverage.set_xlabel('Position in observed RNA (nt)', fontsize=7)
    coverage.tick_params(labelsize=7)
    coverage.spines[['top', 'right']].set_visible(False)
    coverage.set_ylim(bottom=0)
    coverage.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
    destination = directory / candidate['id']
    for extension in ('png', 'svg', 'pdf'):
        figure.savefig(str(destination) + '.' + extension, dpi=180, facecolor='white')
    return Path(str(destination) + '.png')


def variant_path_sections(data, task, figures, p, table, styles):
    """Render a compact summary without treating RNA frames as mutant proteins."""
    title = Paragraph(escape(task.get('title', task['id'])), styles['Heading2'])
    variant = data['variant']
    position = variant['start']
    coordinate = f"{data['resolved_contig']}:{position}-{position + max(1, len(variant['ref'])) - 1}"
    counts = data['allele_counts']
    rows = [(state, counts[state]['fragments'], counts[state]['segments'],
             data['discovery'][state]['truncated'] if state in data['discovery'] else 'Not assembled')
            for state in ('ALT', 'REF', 'OTHER')]
    story = [title, p(escape(f"{data['reference_name']} IGV: {coordinate}; original library: {data['source']}")),
        table(('Allele', 'Fragments', 'Segments', 'Discovery truncated'), rows, (80, 100, 100, 230)),
        p(f"{len(data['paths'])} retained RNA paths, including alternative orientations and contexts. "
          'Counts can overlap across states and paths; do not sum them as molecules or proteins. '
          'Frames are conditional on supplied CDS models. RNA does not establish translation or variant causation.')]
    if data.get('provenance', {}).get('scope'):
        story.append(p(escape(data['provenance']['scope']), 'Small'))
    candidates = sorted(data['paths'], key=lambda row: (
        -any(f['state'] == 'ALT' for f in row['focal_alleles']),
        -bool(row['coding']), -row['linkage']['single_read_whole_path_fragments'],
        -row['supporting_fragments'], -len(row['sequence']), row['id']))[:4]
    directory = figures / 'variant_rna_paths' / task['id']
    directory.mkdir(parents=True, exist_ok=True)
    for candidate in candidates:
        linkage = candidate['linkage']
        states = ', '.join(sorted({f['state'] for f in candidate['focal_alleles']}))
        description = (f"{states} context: {len(candidate['sequence'])} nt; "
            f"{candidate['supporting_fragments']} contributing fragments; "
            f"{linkage['single_read_whole_path_fragments']} whole-path single-read fragments; "
            f"{linkage['paired_template_only_whole_path_fragments']} paired-only whole-path templates.")
        frames = candidate['coding']
        outcomes = ', '.join(sorted({frame['focal_coverage'].replace('_', ' ') for frame in frames}))
        frame_text = ('Conditional frame / allele relation: ' + outcomes + '.' if frames else
                      'No CDS frame was assigned to this RNA path. ' + candidate['protein_interpretation'].replace('_', ' ') + '.')
        story.append(KeepTogether([p('<b>' + escape(description) + '</b>'),
            Image(str(path_figure(candidate, directory)), width=510, height=198), p(escape(frame_text), 'Small')]))
    return story
