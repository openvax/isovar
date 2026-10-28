"""Source-separated landscapes and sequence exemplars from cohort results."""

from collections import Counter, defaultdict
import gzip
from html import escape
import json
import math
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import (
    Flowable, Image, KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle,
)

from .sv_rna_audit import file_sha256, read_json, write_json
from .rna_evidence import content_identifier


INK = "#223743"
MUTANT = "#d26928"
GERMLINE = "#773f9c"
COSOMATIC = "#187eaa"
UNKNOWN = "#687780"
TEAL = "#128681"
PALE = "#edf4f4"


def _counts(row):
    event = row["protein_hypotheses"]["events"][0]
    return event, event["allele_support"]["alt"]["fragments"]


def _best_translation(protein):
    return max(protein["translations"], key=lambda t: (
        t["reference_comparison"]["focal_amino_acid_change"],
        t["coding_support"]["allele_observation_span_support"]["fragments"],
        t["coding_support"]["minimum_fragment_depth"], t["translation_id"]))


def _collect(result, highlights=()):
    counts, categories, depths = Counter(), Counter(), defaultdict(list)
    examples = defaultdict(list)
    if result.get("status") not in ("completed", "completed_with_errors", "in_progress"):
        return counts, categories, depths, examples
    path = Path(result["output"]) / result["results"]
    if file_sha256(path) != result["result_sha256"]:
        raise ValueError("Report input digest mismatch: " + str(path))
    counts["nominated"] = result["nominated_alleles"]
    with gzip.open(path, "rt") as stream:
        for line in stream:
            row = json.loads(line)
            counts["processed"] += 1
            if row["status"] != "assessed":
                counts["unresolved_chromosomes" if row["status"] == "unresolved_chromosome" else "errors"] += 1
                continue
            event, alt = _counts(row)
            category = row["category"]
            categories[category] += 1
            depths[category].append(alt)
            counts["alt"] += alt > 0
            counts["alt20"] += alt >= 20
            proteins = event["protein_hypotheses"]
            counts["protein"] += bool(proteins)
            counts["exon"] += row["exon_overlap"]
            changed = [p for p in proteins if any(t["reference_comparison"]["focal_amino_acid_change"]
                                                for t in p["translations"])]
            counts["origin_conflicts"] += bool(row["nomination"].get("origin_conflict"))
            counts["changed_protein"] += bool(changed) and not row["nomination"].get("origin_conflict", False)
            eligible = not row["nomination"].get("origin_conflict", False)
            counts["changed_protein_span20"] += eligible and any(
                t["coding_support"]["allele_observation_span_support"]["fragments"] >= 20
                for protein in changed for t in protein["translations"]
                if t["reference_comparison"]["focal_amino_acid_change"])
            counts["changed_protein_single20"] += eligible and any(
                t["coding_support"]["single_read_span_support"]["fragments"] >= 20
                for protein in changed for t in protein["translations"]
                if t["reference_comparison"]["focal_amino_acid_change"])
            nc = row["nc_orfs"].get("result")
            counts["nc_overlap"] += bool(row["nc_orfs"].get("overlapping_models"))
            counts["nc_unresolved_models"] += len(row["nc_orfs"].get("unresolved_models", []))
            changed_nc = [p for p in nc["protein_candidates"] if p["focal_amino_acid_change"]] if nc else []
            counts["changed_nc"] += bool(changed_nc) and not row["nomination"].get("origin_conflict", False)
            counts["changed_nc_span20"] += eligible and any(
                p["coding_support"]["allele_observation_span_support"]["fragments"] >= 20 for p in changed_nc)
            # The full export keeps witnesses and every alternative. Exemplar
            # selection needs only the selected protein/translation, not its
            # potentially large evidence-set ledger.
            chosen = next(iter(changed), proteins[0] if proteins else None)
            compact_proteins = []
            if chosen:
                best = _best_translation(chosen)
                compact_proteins = [dict(chosen, translations=[best])]
            displayed = dict(row, protein_hypotheses=dict(events=[dict(
                event, protein_hypotheses=compact_proteins)]), nc_orfs={})
            v = row['nomination']
            key = f"{str(v['contig']).removeprefix('chr')}:{v['start']}:{v['ref']}>{v['alt']}"
            if key in highlights:
                examples['Highlighted loci'].append((alt, displayed, None))
            if alt:
                examples["Normal-call overlap" if row["nomination"].get("origin_conflict") else
                         category if changed else "ALT without changed protein"].append((alt, displayed, None))
            if changed_nc:
                candidate = max(changed_nc, key=lambda p: (
                    p["coding_support"]["allele_observation_span_support"]["fragments"],
                    p["coding_support"]["minimum_fragment_depth"]))
                reference = next(m["reference"] for m in nc["models"]
                                 if m["model_id"] == candidate["model_id"])
                candidate = dict(candidate, reference_evidence=reference.get("evidence", {}),
                                 reference_source=reference.get("source", ""))
                if not row["nomination"].get("origin_conflict"):
                    examples["ncORF"].append((candidate["coding_support"]["minimum_fragment_depth"], displayed, candidate))
            # Bound report memory, without discarding any analysis export.
            for key in examples:
                if key != 'Highlighted loci' and len(examples[key]) > 20:
                    examples[key] = sorted(examples[key], key=lambda r: -r[0])[:10]
    return counts, categories, depths, examples


def _landscape(sample, analyses, directory, index=1):
    from matplotlib.figure import Figure

    figure = Figure(figsize=(7.2, max(3.2, 1.6 + len(analyses) * .42)), layout="constrained")
    ax = figure.subplots()
    labels = []
    for i, (result, counts, _, _, _) in enumerate(analyses):
        labels.append(result["source_id"] + " / " + result["pass_id"])
        for offset, key, label, color in [(-.23, "nominated", "Nominated", "#a4b6bd"),
                                         (0, "alt", "RNA ALT", TEAL),
                                         (.23, "changed_protein", "Changed protein context", MUTANT)]:
            value = counts[key]
            ax.barh(i + offset, value, height=.21, color=color, label=label if i == 0 else None)
            ax.text(max(value, 1), i + offset, f" {value:,}", fontsize=7, va="center", color=INK)
    ax.set_yticks(range(len(labels)), labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xscale("symlog", linthresh=1)
    ax.set_xlabel("Variant counts within each library (log scale above 1)", fontsize=8)
    ax.set_title(sample + " | RNA evidence funnel", loc="left", fontsize=12, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    prefix = directory / ("landscape_" + str(index))
    for extension in ("svg", "pdf", "png"):
        figure.savefig(str(prefix) + "." + extension, dpi=180, facecolor="white")
    return Path(str(prefix) + ".png"), figure.get_size_inches()



def _category_landscape(result, categories, depths, directory):
    from matplotlib.figure import Figure
    import numpy as np

    fig = Figure(figsize=(7.2, 3.6), layout="constrained")
    left, right = fig.subplots(1, 2)
    labels = list(categories)
    positions = range(len(labels))
    left.barh(positions, [categories[k] for k in labels], color="#a4b6bd", label="Evaluated")
    left.barh(positions, [sum(n >= 20 for n in depths[k]) for k in labels], color=TEAL, label="ALT >=20 fragments")
    left.set_yticks(positions, labels, fontsize=7)
    left.set_xscale("symlog", linthresh=1)
    left.set_xlabel("Calls (log scale above 1)", fontsize=8)
    left.legend(frameon=False, fontsize=7, loc="lower right")
    palette = [INK, MUTANT, TEAL, GERMLINE, COSOMATIC]
    for color, label in zip(palette, labels):
        values = np.array(sorted(n for n in depths[label] if n > 0))
        if len(values):
            right.step(values, np.arange(len(values), 0, -1), where="post", label=label,
                       color=color, linewidth=1.4)
    right.set_xscale("log")
    right.set_yscale("log")
    right.set_xlabel("Minimum RNA ALT fragments", fontsize=8)
    right.set_ylabel("Calls at or above support", fontsize=8)
    right.legend(frameon=False, fontsize=7)
    for ax in (left, right):
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=7)
    fig.suptitle(result["source_id"] + " | variant landscape", fontsize=11, color=INK, x=.02, ha="left")
    prefix = directory / (result["source_id"] + "-" + result["pass_id"] + "-categories")
    for ext in ("png", "svg", "pdf"):
        fig.savefig(str(prefix) + "." + ext, dpi=180, facecolor="white")
    return Path(str(prefix) + ".png")

def aligned_columns(reference, protein, interval, *, allele_length_delta=0, frameshift=False,
                    stop=False, reference_stop=False):
    """Display codon-anchored differences without re-aligning a frameshift tail.

    In-frame focal gaps come from the declared allele length and translated
    mutation interval. This preserves the nucleotide anchor in repeats.
    """
    begin, end = interval
    delta = allele_length_delta // 3 if not frameshift and allele_length_delta % 3 == 0 else 0
    reference_length = len(reference)
    reference = reference[:len(protein) + max(delta, 0)]
    columns = []

    def add(r, q, deletion=False):
        columns.append(dict(reference=reference[r] if r is not None and r < len(reference) else "-",
                            protein=protein[q] if q is not None and q < len(protein) else "-",
                            position=q, reference_position=r, focal_deletion=deletion))

    for i in range(min(begin, len(protein))):
        add(i, i)
    reference_end = max(begin, end + delta)
    for i in range(max(reference_end - begin, end - begin)):
        r = begin + i if begin + i < reference_end else None
        q = begin + i if begin + i < end else None
        add(r, q, deletion=q is None and r is not None)
    for q in range(max(begin, end), len(protein)):
        add(q + delta, q)
    if stop:
        r = len(protein) + delta
        columns.append(dict(reference="*" if reference_stop and r == reference_length else
                            reference[r] if r < len(reference) else "-", protein="*",
                            position=len(protein), reference_position=r, focal_deletion=False))
    return columns


class SequencePanel(Flowable):
    """Aligned residues, explicit focal gaps and exact original-RNA depth."""

    def __init__(self, columns, interval, depth, *, origins=None):
        super().__init__()
        self.columns, self.interval, self.depth = columns, interval, depth
        self.origins = origins or {}
        self.width, self.height = 510, 97

    def draw(self):
        canvas = self.canv
        canvas.setFillColor(colors.HexColor(PALE))
        canvas.roundRect(0, 0, self.width, self.height - 3, 5, fill=1, stroke=0)
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(colors.HexColor(UNKNOWN))
        for label, y in (("REF", 68), ("RNA", 49), ("Depth", 16)):
            canvas.drawString(8, y, label)
        positions = [c["position"] for c in self.columns if c["position"] is not None]
        if positions:
            canvas.drawString(45, 81, f"RNA protein-window residues {min(positions) + 1}-{max(positions) + 1}")
        canvas.setFont("Courier", 8.4)
        left, step = 45, 5.15
        begin, end = self.interval
        codons = [min(self.depth[i:i + 3]) for i in range(0, len(self.depth), 3) if self.depth[i:i + 3]]
        maximum = max(codons, default=0)
        for i, column in enumerate(self.columns):
            position = column["position"]
            canvas.setFillColor(colors.HexColor(INK))
            canvas.drawString(left + step * i, 68, column["reference"])
            same = column["protein"] == column["reference"]
            origin = self.origins.get(str(position), "unclassified")
            focal = (position is not None and begin <= position < end) or column["focal_deletion"]
            if column["protein"] == "*" and position == begin:
                focal = True
            color = (INK if same else MUTANT if focal else
                     GERMLINE if origin == "germline_variant" else
                     COSOMATIC if origin == "co_somatic_variant" else UNKNOWN)
            canvas.setFillColor(colors.HexColor(color))
            if column["focal_deletion"]:
                canvas.setStrokeColor(colors.HexColor(MUTANT))
                canvas.setLineWidth(1.6)
                x = left + step * i + step / 2
                canvas.line(x, 45, x, 59)
            else:
                canvas.drawString(left + step * i, 49, column["protein"])
            if position is not None and position < len(codons):
                count = codons[position]
                canvas.setFillColor(colors.HexColor(TEAL if count else UNKNOWN))
                height = 16 * math.log1p(count) / math.log1p(maximum) if maximum else 0
                canvas.rect(left + step * i, 14, step - .5, height, fill=1, stroke=0)
        canvas.setFont("Helvetica", 6.5)
        canvas.setFillColor(colors.HexColor(UNKNOWN))
        canvas.drawString(left, 4, f"Exact fragments per codon; min {min(codons, default=0)}, max {maximum}; log height")


class FusionSequencePanel(Flowable):
    """Two parental comparisons keep ordinary acceptor residues neutral."""

    def __init__(self, protein, parents, depth, junction, first, acceptor_in_frame=None):
        super().__init__()
        self.protein, self.parents, self.depth = protein, parents, depth
        self.junction, self.first = junction, first
        self.acceptor_in_frame = acceptor_in_frame
        self.width, self.height = 510, 110

    def draw(self):
        c = self.canv
        c.setFillColor(colors.HexColor(PALE))
        c.roundRect(0, 0, 510, 105, 5, fill=1, stroke=0)
        c.setFillColor(colors.HexColor(UNKNOWN))
        c.setFont('Helvetica', 7)
        c.drawString(45, 94, f'Fusion protein residues {self.first + 1}-{self.first + len(self.protein)}')
        for label, y in [('REF D', 79), ('REF A', 62), ('RNA', 45), ('Depth', 15)]:
            c.drawString(7, y, label)
        c.setFont('Courier', 8.4)
        maximum = max(self.depth, default=0)
        for i, amino in enumerate(self.protein):
            q = self.first + i
            donor, acceptor = (p[q] if q < len(p) else ' ' for p in self.parents)
            for aa, y in [(donor, 79), (acceptor, 62)]:
                c.setFillColor(colors.HexColor(INK))
                c.drawString(45 + 5.15 * i, y, aa)
            changed = amino not in (donor, acceptor)
            begin, end = self.junction
            focal = (q * 3 < end and (q + 1) * 3 > begin or
                     self.acceptor_in_frame is False and (q + 1) * 3 > end)
            c.setFillColor(colors.HexColor(MUTANT if changed and focal else UNKNOWN if changed else INK))
            c.drawString(45 + 5.15 * i, 45, amino)
            if q < len(self.depth):
                count = self.depth[q]
                c.setFillColor(colors.HexColor(TEAL))
                height = 16 * math.log1p(count) / math.log1p(maximum) if maximum else 0
                c.rect(45 + 5.15 * i, 12, 4.6, height, fill=1, stroke=0)
        for edge in sorted(set(self.junction)):
            boundary = edge / 3 - self.first
            if 0 <= boundary <= len(self.protein):
                c.setStrokeColor(colors.HexColor(MUTANT))
                c.setDash(2, 2)
                x = 45 + 5.15 * boundary
                c.line(x, 39, x, 89)
                c.setDash()
        c.setFont('Helvetica', 6.5)
        c.setFillColor(colors.HexColor(UNKNOWN))
        c.drawString(45, 3, 'Exact observed fragments per codon; no reference-filled sequence')


def _fusion_comparisons(result, supplied, depth, whole_fragments):
    """Project parental codons onto observed RNA without filling insertions."""
    from .fusion import FusionReference
    from .genetic_code import translate_cdna

    path = result['paths'][0]
    junction_start, junction_end = path['junctions'][0]['query_interval']
    refs = {(r['annotation'], r['transcript_id']): FusionReference(**r) for r in supplied['references']}
    id_counts = Counter(ref.transcript_id for ref in refs.values())

    def reference_name(ref):
        suffix = '--' + content_identifier('annotation', ref.annotation)[-12:] if id_counts[ref.transcript_id] > 1 else ''
        return 'reference_' + ref.transcript_id + suffix

    panels, fasta = [], []
    for number, translation in enumerate(path['translations'], 1):
        protein = translation['amino_acids']
        donor_frame = translation['frame_evidence'][0]
        donor = refs[(donor_frame['annotation'], donor_frame['transcript_id'])]
        acceptor_frames = translation['acceptor_frames']
        acceptors = [ref for ref in refs.values() if acceptor_frames and
                     ref.transcript_id == acceptor_frames[0]['transcript_id'] and ref.cds_start is not None]
        # Legacy acceptor-frame records do not identify the annotation. Do not
        # arbitrarily pick among references with the same transcript ID.
        acceptor = acceptors[0] if len(acceptors) == 1 else None
        acceptor_status = 'available' if acceptor else 'ambiguous' if acceptors else 'unavailable'
        donor_protein, _ = translate_cdna(donor.sequence[donor.cds_start:donor.cds_end])
        acceptor_protein = ''
        if acceptor is not None:
            acceptor_protein, _ = translate_cdna(acceptor.sequence[acceptor.cds_start:acceptor.cds_end])
            acceptor_start = acceptor.spliced_offset(path['junctions'][0]['right'][1]) - acceptor.cds_start
        q0 = translation['translation_start']
        donor_start = donor_frame['donor_reference_offset'] - donor.cds_start
        parent_donor, parent_acceptor = [], []
        for i in range(len(protein)):
            q = q0 + 3 * i
            d = (donor_start + q) // 3
            parent_donor.append(donor_protein[d] if q < junction_start and 0 <= d < len(donor_protein) else ' ')
            # A split codon must compare to the first contributing acceptor
            # base, not an extrapolated upstream codon absent from this path.
            a = (acceptor_start + max(q, junction_end) - junction_end) // 3 if acceptor is not None else -1
            parent_acceptor.append(acceptor_protein[a] if q + 3 > junction_end and 0 <= a < len(acceptor_protein) else ' ')
        codons = [min(depth[q:q + 3]) for q in range(q0, translation['translation_end'], 3)]
        panels.append(dict(translation=translation, parents=(parent_donor, parent_acceptor), depth=codons,
                           junction=[junction_start - q0, junction_end - q0], whole_fragments=whole_fragments,
                           donor=donor.transcript_id, donor_annotation=donor.annotation,
                           acceptor=acceptor.transcript_id if acceptor else None, acceptor_status=acceptor_status,
                           acceptor_annotation=acceptor.annotation if acceptor else None,
                           acceptor_in_frame=acceptor_frames[0]['in_frame'] if acceptor else None))
        fasta.extend([(f'fusion_{number}', protein), (reference_name(donor), donor_protein)])
        if acceptor is not None:
            fasta.append((reference_name(acceptor), acceptor_protein))
    return panels, dict(fasta)


def _fusion_display(result, supplied, directory):
    from matplotlib.figure import Figure

    path = result['paths'][0]
    sequence = path['sequence']
    coverage = [set() for _ in sequence]
    whole = set()
    for read in result['evidence']['observations']:
        start, end = read['cdna_start'], read['cdna_start'] + len(read['sequence'])
        if sequence[start:end] != read['sequence']:
            raise ValueError('Fusion observation does not reproduce its reported placement')
        fragment = (read['sample_id'], read['library_id'], read['fragment_id'])
        for q in range(start, end):
            coverage[q].add(fragment)
        if start == 0 and end == len(sequence):
            whole.add(fragment)
    depth = [len(v) for v in coverage]
    junctions = sorted(set(path['junctions'][0]['query_interval']))
    fig = Figure(figsize=(7.2, 3.2), layout='constrained')
    ax = fig.subplots()
    ax.fill_between(range(len(depth)), depth, color=TEAL, alpha=.8, step='mid')
    for i, junction in enumerate(junctions):
        ax.axvline(junction, color=MUTANT, linewidth=1.6, label='Fusion junction' if i == 0 else None)
    ax.set_yscale('symlog', linthresh=1)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Observed fusion cDNA position (nt)', fontsize=8)
    ax.set_ylabel('Exact original fragments', fontsize=8)
    ax.set_title(result['event_id'].replace('--', ' / ') + ' | observed RNA coverage', loc='left', fontsize=11)
    ax.legend(frameon=False, fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    prefix = directory / 'coverage'
    for ext in ('png', 'svg', 'pdf'):
        fig.savefig(str(prefix) + '.' + ext, dpi=180, facecolor='white')
    panels, fasta = _fusion_comparisons(result, supplied, depth, len(whole))
    with (directory / 'proteins.fasta').open('w') as handle:
        for name, sequence in fasta.items():
            handle.write('>' + name + '\n' + '\n'.join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + '\n')
    return Path(str(prefix) + '.png'), panels


def _sv_path_figure(candidate, assessment, directory, identifier):
    """Show exact path matching without equating tiling with joint phase."""
    from matplotlib.figure import Figure
    import numpy as np

    length = len(candidate['sequence'])
    intervals = defaultdict(list)
    for match in assessment.get('exact_matches', []):
        intervals[tuple(match['fragment_id'])].append(match['path_interval'])
    changes = np.zeros(length + 1, dtype=int)
    for spans in intervals.values():
        merged = []
        for start, end in sorted(spans):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(end, merged[-1][1])
            else:
                merged.append([start, end])
        for start, end in merged:
            changes[max(0, start)] += 1
            changes[min(length, end)] -= 1
    depth = np.cumsum(changes[:-1])
    figure = Figure(figsize=(7.2, 2.7), layout='constrained')
    blocks, coverage = figure.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [1, 2]})
    for index, block in enumerate(candidate['blocks']):
        a, b = block['query_start'], block['query_end']
        blocks.broken_barh([(a, b - a)], (.3, .3), facecolors=TEAL)
        if index in (0, len(candidate['blocks']) - 1):
            label = f"{block['contig']}:{block['reference_start'] + 1}-{block['reference_end']} ({block['strand']})"
            blocks.text(a if index == 0 else b, .78, label, fontsize=6.5,
                        ha='left' if index == 0 else 'right', color=INK)
    blocks.set_ylim(0, 1)
    blocks.axis('off')
    coverage.fill_between(np.arange(length), depth, color=TEAL, alpha=.8)
    for junction in candidate['junctions']:
        if junction['relation'] in ('breakpoint_junction', 'event_compatible_junction'):
            a, b = junction['query_interval']
            for ax in (blocks, coverage):
                ax.axvline(a, color=MUTANT, linewidth=1)
                if b > a:
                    ax.axvspan(a, b, color=MUTANT, alpha=.2)
    coverage.set_yscale('symlog', linthresh=1)
    coverage.set_ylim(bottom=0)
    coverage.set_ylabel('Exact matching fragments', fontsize=8)
    coverage.set_xlabel('Position in the selected RNA path (nt)', fontsize=8)
    coverage.tick_params(labelsize=7)
    coverage.spines[['top', 'right']].set_visible(False)
    prefix = directory / identifier
    for extension in ('png', 'svg', 'pdf'):
        figure.savefig(str(prefix) + '.' + extension, dpi=180, facecolor='white')
    return Path(str(prefix) + '.png')


def _additional_sections(manifest, summary, sample_id, figures, p, table, styles):
    from .antigen_report import fusion_summary, viral_summary, isoform_summary

    story = []
    for entry in summary.get('additional', []):
        task = entry['task']
        if task.get('sample_id') != sample_id:
            continue
        heading = Paragraph(escape(task.get("title", task["id"])), styles["Heading2"])
        if entry['status'] != 'recorded':
            story.append(KeepTogether([heading, p('This analysis has not completed. '
                'No absence conclusion is available. ' + escape(entry.get('reason', '')))]))
            continue
        path = Path(entry['path'])
        if file_sha256(path) != entry['sha256']:
            raise ValueError('Additional report input digest mismatch')
        data = read_json(path)
        family = entry['family']
        if family == 'fusions':
            counts = fusion_summary(data)
            if counts['translated_candidates']:
                story.append(PageBreak())
            block = [heading]
            if not data['paths']:
                block.append(p('No retained RNA path. ' + escape('; '.join(data.get('reasons', [])))))
                story.append(KeepTogether(block))
                continue
            j = data['paths'][0]['junctions'][0]
            coords = ' | '.join(f"chr{point[0].removeprefix('chr')}:{point[1] + 1}-{point[1] + 1}"
                                for point in (j['left'], j['right']))
            block.extend([p(f"{counts['direct_junction_fragments']:,} direct-junction fragments; "
                            f"{counts['translated_candidates']} coding candidates. Original RNA library: " +
                            escape(task['source_id'])),
                          p(escape(data['reference_name'] + ' IGV anchors: ' + coords), 'Small')])
            descriptor = manifest['files'][task['input']]
            input_path = Path(summary['manifest']).parent / descriptor['path']
            if file_sha256(input_path) != descriptor['sha256']:
                raise ValueError('Fusion reference input digest mismatch')
            task_figures = figures / 'fusions' / task['id']
            task_figures.mkdir(parents=True, exist_ok=True)
            plot, panels = _fusion_display(data, read_json(input_path), task_figures)
            block.append(Image(str(plot), width=510, height=226))
            if not panels:
                block.append(p('Coding frame unresolved; no protein sequence is inferred. ' +
                               escape('; '.join(data.get('reasons', []))), 'Small'))
            story.append(KeepTogether(block))
            for panel in panels:
                translation = panel['translation']
                protein = translation['amino_acids']
                story.append(p(f"<b>{len(protein)} aa</b>; " +
                               ('annotated start observed' if translation['complete_5prime'] else 'partial 5-prime sequence') +
                               ('; observed stop.' if translation['ends_with_stop_codon'] else '; stop not reached.') +
                               f" Minimum coding depth: {min(panel['depth'], default=0)} fragments per codon. "
                               f"{panel['whole_fragments']} original fragments individually span the entire assembled cDNA. "
                               "A fully tiled coding candidate is not a single observed full-length molecule."))
                acceptor_label = (panel['acceptor'] + ' (' + panel['acceptor_annotation'] + ')' if panel['acceptor'] else
                                  'ambiguous (transcript ID occurs in multiple annotations)' if panel['acceptor_status'] == 'ambiguous' else
                                  'unavailable (no supported coding acceptor model)')
                if panel['acceptor_in_frame'] is False:
                    acceptor_label += ' (shifted reading frame)'
                story.append(p('REF D: ' + escape(panel['donor'] + ' (' + panel['donor_annotation'] + ')') + '; REF A: ' +
                               escape(acceptor_label) +
                               '. Each available parental protein is shown only on its contributing side. Dashed boundaries mark '
                               'the novel adjacency or inserted interval; matching parental residues remain neutral. Available parental proteins '
                               'and the complete candidate are in the companion FASTA.', 'Small'))
                if panel['acceptor_in_frame'] is False:
                    story.append(p('Reference rows follow source codons; a split reference codon can span two '
                                   'columns when the reading frames differ.', 'Small'))
                for start in range(0, len(protein), 80):
                    story.append(FusionSequencePanel(protein[start:start + 80], panel['parents'],
                                                     panel['depth'], panel['junction'], start,
                                                     acceptor_in_frame=panel['acceptor_in_frame']))
                    story.append(Spacer(1, 7))
        elif family == 'sv_audits':
            story.extend([PageBreak(), heading])
            runs = data['runs']
            annotations = {}
            if task.get('locus_annotations'):
                descriptor = manifest['files'][task['locus_annotations']]
                annotation_path = Path(summary['manifest']).parent / descriptor['path']
                if file_sha256(annotation_path) != descriptor['sha256']:
                    raise ValueError('SV locus annotation digest mismatch')
                annotations = read_json(annotation_path)
            if data.get('status') == 'in_progress':
                story.append(p('<b>Interim SV snapshot.</b> Only completed event/source/orientation views are shown.'))
            by_source = defaultdict(list)
            for row in runs:
                by_source[row['source_id']].append(row)
            rows = []
            for source, values in by_source.items():
                ids = {r['event_id'] for r in values}
                linked = {r['event_id'] for r in values if r['paths_with_direct_junction']}
                coding = {r['event_id'] for r in values if r['paths_with_event_linked_protein']}
                unavailable = {r['event_id'] for r in values if r['status'] in ('unavailable', 'unassessed')}
                rows.append((source, len(ids), len(linked), len(coding), len(unavailable)))
            story.extend([p('Source VCF records and orientation views remain separate. Different records can describe '
                            'the same biological event. Direct junction, exon linkage and coding frame are independent evidence states.'),
                          table(('Source', 'Records', 'Direct junction', 'Event-linked protein', 'Unassessed'),
                                rows, (170, 70, 90, 95, 85))])
            strongest = {}
            for row in (r for r in runs if r['paths_with_direct_junction']):
                snapshot = read_json(path.parent / row['artifact'])
                reconstruction = snapshot['reconstruction']
                supported = sorted(reconstruction['paths'], key=lambda x: -max(
                    (j['direct_support'].get('fragments', 0) for j in x['junctions']
                     if j['relation'] in ('breakpoint_junction', 'event_compatible_junction')), default=0))
                if not supported:
                    continue
                candidate = supported[0]
                assessed = next(p for p in snapshot['assessment']['paths'] if p['path_id'] == candidate['path_id'])
                count = max((j['direct_support'].get('fragments', 0) for j in candidate['junctions']
                             if j['relation'] in ('breakpoint_junction', 'event_compatible_junction')), default=0)
                score = (any(p['event_linked'] for p in assessed['proteins']), count)
                key = (row['event_id'], row['source_id'])
                if key not in strongest or score > strongest[key][0]:
                    strongest[key] = (score, row, candidate, assessed, reconstruction.get('reference_name', ''))
            selected = sorted(strongest.values(), key=lambda x: (-x[0][0], -x[0][1], x[1]['event_id']))[:5]
            for display_index, (_, row, candidate, assessed, reference_name) in enumerate(selected):
                junctions = [j for j in candidate['junctions'] if j['relation'] in
                             ('breakpoint_junction', 'event_compatible_junction')]
                coord = []
                for j in junctions:
                    for point in (j.get('left'), j.get('right')):
                        if point:
                            coord.append(f"chr{point[0].removeprefix('chr')}:{point[1] + 1}-{point[1] + 1}")
                transcripts = sorted({t for model in candidate.get('reference_readings', [])
                                      for t in model.get('transcript_ids', [])})
                count = max((j['direct_support'].get('fragments', 0) for j in junctions), default=0)
                locus = annotations.get(row['event_id'], {})
                title = (' / '.join(locus.get('genes', [])) or 'Unassigned gene') + ' | ' + row['event_id']
                block = [p('<b>' + escape(title) + '</b>'),
                    p('Nominated GRCh38 event: ' + escape(' | '.join(locus.get('nominal_igv', []))), 'Small'),
                    p(escape(row['source_id'] + ' / ' + row['view_id']) + f'; strongest event-linked junction: {count} fragments. '
                      + 'Frame state: ' + escape(candidate['frame_status']), 'Small'),
                    p(escape(reference_name + ' IGV: ' + ' | '.join(dict.fromkeys(coord))), 'Small'),
                    p('Reference reading models: ' + escape(', '.join(transcripts) or 'unresolved') +
                      '. All retained paths and alternatives remain in the SV audit export.', 'Small')]
                plot = _sv_path_figure(candidate, assessed, figures, task['id'] + '-path-' + str(display_index + 1))
                block.extend([Image(str(plot), width=510, height=191),
                              p('Orange marks the event-associated junction. Coverage counts exact sequence matches '
                                'from distinct fragments at each position; tiling does not establish long-range phase. '
                                f"{len(assessed['full_sequence_support']['complete_fragments'])} fragments span the entire selected path; "
                                f"{len(assessed['exact_exon_anchors'])} exact reference-exon anchors. "
                                'Orientation views are consolidated for display only.', 'Small')])
                story.append(KeepTogether(block))
        elif family in ('viral', 'isoforms'):
            counts = viral_summary(data) if family == 'viral' else isoform_summary(data)
            story.append(KeepTogether([heading, p(escape(json.dumps(counts, sort_keys=True)))]))
        else:
            story.append(KeepTogether([heading,
                p('Prior assessment retained with its original source scope; this is not a new whole-library screen.'),
                p(escape(task.get('summary', json.dumps(data, sort_keys=True)[:1200])))]))
    return story

def write_cohort_reports(manifest, summary, output):
    """Render per-sample PDFs and standalone figures without rerunning analysis."""
    output = Path(output)
    directory = output / "pdf"
    directory.mkdir(exist_ok=True)
    styles = getSampleStyleSheet()
    for name, size, leading in [("Body", 9, 13), ("Small", 7.5, 10), ("Cell", 7.2, 9.5)]:
        styles.add(ParagraphStyle(name="Cohort" + name, fontName="Helvetica", fontSize=size,
                                  leading=leading, textColor=colors.HexColor(INK), spaceAfter=7))
    styles["Title"].textColor = colors.HexColor(INK)
    styles["Heading2"].textColor = colors.HexColor(TEAL)
    styles["Heading2"].keepWithNext = True

    def p(text, style="Body"):
        return Paragraph(str(text), styles["Cohort" + style])

    def table(headers, rows, widths):
        values = [[p("<b>" + escape(str(x)) + "</b>", "Cell") for x in headers]]
        values.extend([[p(escape(str(x)), "Cell") for x in row] for row in rows])
        result = Table(values, colWidths=widths, repeatRows=1, hAlign="LEFT")
        result.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(PALE)),
            ("LINEBELOW", (0, 0), (-1, 0), .7, colors.HexColor(TEAL)),
            ("LINEBELOW", (0, 1), (-1, -1), .3, colors.HexColor("#d3dddd")),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("TOPPADDING", (0, 0), (-1, -1), 5), ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ]))
        return result

    def footer(canvas, doc):
        canvas.setStrokeColor(colors.HexColor("#ccd9db"))
        canvas.line(42, 36, 570, 36)
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(colors.HexColor(UNKNOWN))
        canvas.drawString(42, 24, "RNA-supported sequence hypotheses; translation and presentation are not established.")
        canvas.drawRightString(570, 24, str(doc.page))

    paths = []
    for sample in manifest["samples"]:
        sample_id = sample["sample_id"]
        source_ids = set(sample.get("source_ids", [s["source_id"] for s in summary["sources"]
                                                   if s["sample_id"] == sample_id]))
        sources = [s for s in summary["sources"] if s["source_id"] in source_ids]
        results = [r for r in summary["small_variants"] if r["source_id"] in source_ids]
        analyses = [(r, *_collect(r, sample.get('highlight_variants', []))) for r in results
                    if r.get('status') in ('completed', 'completed_with_errors', 'in_progress')]
        figures = output / "figures" / sample_id
        figures.mkdir(parents=True, exist_ok=True)
        write_json(figures / "summary.json", [dict(source_id=r["source_id"], pass_id=r["pass_id"],
                   counts=dict(c), categories=dict(categories)) for r, c, categories, _, _ in analyses])
        report = directory / (sample_id + "_RNA_evidence.pdf")
        story = [Paragraph(escape(sample.get("title", sample_id)) + "<br/>RNA evidence report", styles["Title"]),
                 p(escape(manifest.get("description", "Source-separated Isovar cohort analysis."))),
                 p("<b>Scope.</b> Each library and analysis pass retains its own denominator. Alternative "
                   "transcripts, processing products and sensitivity passes are not independent replication."),
                 Paragraph("Source coverage", styles["Heading2"])]
        pending = [r for r in manifest.get("pending_inputs", []) if r.get("sample_id") == sample_id]
        if pending or any(r.get("status") == "in_progress" for r in results):
            story.append(p("<b>Interim report.</b> Processing and source acquisition are still in progress. "
                           "The figures count only evaluated observations; pending data cannot establish absence."))
        for item in pending:
            story.append(p(escape(item.get("source", "Additional sources") + ": " + item["state"]), "Small"))
        source_rows = []
        for source in sources:
            matching = [r for r in results if r["source_id"] == source["source_id"]]
            status = "; ".join(r["pass_id"] + ": " + r["status"] for r in matching)
            if not status:
                status = '; '.join(entry['family'] + ': artifact available'
                    for entry in summary.get('additional', [])
                    if entry['task'].get('source_id') == source['source_id'] and entry['status'] == 'recorded')
            source_rows.append((source.get("label", source["source_id"]), source["sample_id"],
                                status or source.get("unavailable_reason", "Not yet run"), source["scope"]))
        story.append(table(("RNA source", "Sample / timepoint", "Analysis state", "Input scope"),
                           source_rows, (130, 100, 130, 150)))
        if analyses:
            for index in range(0, len(analyses), 8):
                path, size = _landscape(sample_id, analyses[index:index + 8], figures, index=index // 8 + 1)
                story.append(Image(str(path), width=510, height=510 * size[1] / size[0]))
            story.extend([p("Counts remain per library and pass. A protein context is a reconstructed "
                            "local sequence; it is not automatically a complete mutant protein."),
                          table(("Source / pass", "Calls", "ALT", "ALT >=20", "Changed protein", "Changed ncORF", "Errors"),
                                [(r["source_id"] + " / " + r["pass_id"], c["nominated"], c["alt"], c["alt20"],
                                  c["changed_protein"], c["changed_nc"], c["errors"])
                                 for r, c, _, _, _ in analyses], (150, 55, 50, 60, 55, 85, 55))])
            story.extend([p("<b>At least 20 supporting fragments.</b> ALT support counts the focal allele. "
                            "Exact spans additionally support the displayed coding path; single-read spans exclude merged mates."),
                          table(("Source / pass", "ALT", "Changed protein spans", "Single-read spans", "Changed ncORF spans"),
                                [(r["source_id"] + " / " + r["pass_id"], c["alt20"],
                                  c["changed_protein_span20"], c["changed_protein_single20"], c["changed_nc_span20"])
                                 for r, c, _, _, _ in analyses], (170, 55, 95, 95, 95))])
        for analysis_index, (result, counts, categories, depths, examples) in enumerate(analyses):
            if analysis_index and counts['processed'] > 1:
                story.append(PageBreak())
            story.extend([Paragraph(escape(result["source_id"] + " / " + result["pass_id"]),
                                                styles["Heading2"]),
                          p(f"{counts['processed']:,}/{counts['nominated']:,} nominations evaluated; {counts['alt']:,} RNA ALT-supported; "
                            f"{counts['nc_overlap']:,} overlap ncORF coding blocks. "
                            f"{counts['nc_unresolved_models']:,} overlapping model bindings remain unresolved. "
                            f"{counts['unresolved_chromosomes']:,} nominations have unresolved chromosomes: "
                            "their RNA support is unmeasured and excluded from depth distributions.")])
            if counts['processed'] > 1:
                story.extend([Image(str(_category_landscape(result, categories, depths, figures)), width=510, height=255),
                          p("These are RNA-informed groups: frameshift labels use annotation predictions or reconstructed "
                            "coding contexts; in-frame labels require a reconstructed context. Other non-SNV includes indels "
                            "without resolved coding consequences. Consequence prediction is run on ALT-supported calls. "
                            "Unprocessed calls have no inferred support state.", "Small"),
                              table(("Variant category", "Evaluated calls"), sorted(categories.items()), (330, 180))])
            displayed_variants = set()
            for category in ("SNV", "Frameshift", "In-frame indel", "Other non-SNV", "ncORF", "ALT without changed protein", "Normal-call overlap", "Highlighted loci"):
                selected = sorted(examples.get(category, []), key=lambda row: -row[0])[:3 if category == "SNV" else 5]
                if category == 'Highlighted loci':
                    selected = [entry for entry in selected if entry[1]['igv'] not in displayed_variants]
                if not selected:
                    continue
                for exemplar_index, (_, row, nc_candidate) in enumerate(selected):
                    displayed_variants.add(row['igv'])
                    event, alt = _counts(row)
                    title = "/".join(row["gene_names"]) or "Unassigned gene"
                    title += " | " + row["igv"]
                    intro = [p("<b>" + escape(title) + "</b>"),
                             p(escape(event["variant"]["description"]) + f" | {alt:,} ALT fragments", "Small")]
                    if exemplar_index == 0:
                        intro.insert(0, Paragraph(escape(category) + " exemplars", styles["Heading2"]))
                    prediction = event.get('reference_prediction')
                    if prediction:
                        intro.append(p('Annotation prediction: ' + escape(prediction['effect_class'] + ' ' +
                            (prediction.get('description') or '')), 'Small'))
                    if row["nomination"].get("origin_conflict"):
                        intro.append(p("<b>Normal-call overlap:</b> " + escape(row["nomination"]["origin_conflict"]) +
                                       ". Tumor specificity is unresolved.", "Small"))
                    if row['nomination'].get('normal_evidence_status'):
                        intro.append(p(escape(row['nomination']['normal_evidence_status']), 'Small'))
                    proteins = event["protein_hypotheses"]
                    comparison = {}
                    if nc_candidate is not None:
                        candidate = nc_candidate
                        protein, interval = candidate["protein_sequence"], candidate["mutation_interval"]
                        reference, stop = candidate["reference_protein_window"], candidate["ends_with_stop_codon"]
                        coverage = candidate["coding_support"]
                        contributors = candidate["rna_support"]["fragments"]
                        intro.append(p(escape(candidate["orf_id"] + " | " + candidate["sequence_extent"]), "Small"))
                        evidence = candidate["reference_evidence"]
                        studies = evidence.get("ribosome_studies", [])
                        summaries = evidence.get("summaries", [])
                        evidence_sources = sorted({s.get("source", "") for s in summaries} - {""})
                        intro.append(p("Reference ORF: " + escape(candidate["reference_source"] +
                            ("; Ribo-seq: " + ", ".join(studies) if studies else "")) +
                            ("; catalogue evidence: " + escape(", ".join(evidence_sources)) if evidence_sources else "") +
                            f"; {len(evidence.get('peptides', []))} imported peptide evidence records. "
                            "These support the reference ORF in other experiments; mutant translation is untested.", "Small"))
                    elif proteins:
                        candidate = next((x for x in proteins if any(
                            t["reference_comparison"]["focal_amino_acid_change"] for t in x["translations"])), proteins[0])
                        translation = _best_translation(candidate)
                        protein, interval = candidate["amino_acids"], candidate["mutation_interval"]
                        coverage = translation["coding_support"]
                        comparison = translation["reference_comparison"]
                        reference = comparison["reference_protein_window"] if comparison else ""
                        stop = candidate["ends_with_stop_codon"]
                        contributors = candidate["rna_support"]["fragments"]
                        intro.append(p(escape(", ".join(candidate["transcript_ids"])) + " | " +
                                       escape(candidate["n_terminus"] + "; " + candidate["c_terminus"]), "Small"))
                    else:
                        intro.append(p(("RNA ALT is observed" if alt else "No RNA ALT support was observed under these settings") +
                                       "; no protein sequence was reconstructed. "
                                       "The full result retains ncORF model outcomes and RNA evidence.", "Small"))
                        story.append(KeepTogether(intro))
                        continue
                    spans = coverage["allele_observation_span_support"]["fragments"]
                    single = coverage["single_read_span_support"]["fragments"]
                    strength = ("High coverage of the displayed interval" if spans >= 20 else
                                "Limited full-window support" if spans >= 2 else
                                "Single-fragment full-window support" if spans == 1 else
                                "Tiling only: no exact full-window observation")
                    intro.append(p("<b>" + strength + "</b>", "Small"))
                    intro.append(p(f"{len(protein)} aa; {contributors:,} context contributors; "
                                   f"{spans:,} exact spans for the displayed nucleotide path, "
                                   f"including {single:,} single-original-read spans.", "Small"))
                    if nc_candidate is not None:
                        v = row["nomination"]
                        delta = len(v["ref"]) - len(v["alt"])
                    else:
                        delta = comparison["reference_allele_length"] - comparison["observed_allele_length"]
                    columns = aligned_columns(reference, protein, interval, allele_length_delta=delta,
                                              frameshift=candidate["frameshift"], stop=stop,
                                              reference_stop=comparison.get("reference_stop_in_context", False))
                    sequence_panels = [SequencePanel(columns[start:start + 80], interval,
                                                     coverage["fragment_depth"],
                                                     origins=comparison.get("amino_acid_origins"))
                                       for start in range(0, len(columns), 80)]
                    if comparison and not comparison["focal_amino_acid_change"]:
                        intro.append(p("The focal allele does not change this displayed amino-acid context.", "Small"))
                    story.append(KeepTogether(intro + sequence_panels[:1]))
                    for panel in sequence_panels[1:]:
                        story.append(panel)
                        story.append(Spacer(1, 7))
        story.extend(_additional_sections(manifest, summary, sample_id, figures, p, table, styles))
        story.extend([Spacer(1, 12), Paragraph("Evidence interpretation and reproducibility", styles["Heading2"]),
                      p(f'<font color="{MUTANT}">&#9632;</font> Focal difference / deletion boundary &nbsp; '
                        f'<font color="{GERMLINE}">&#9632;</font> Supplied germline &nbsp; '
                        f'<font color="{COSOMATIC}">&#9632;</font> Co-somatic &nbsp; '
                        f'<font color="{UNKNOWN}">&#9632;</font> Unclassified difference'),
                      p("Reference-matching residues remain neutral. Context contributors need not span the "
                        "displayed coding window. Exact allele-observation spans can use overlapping merged mates; "
                        "single-original-read spans are reported separately. Coverage is measured on the selected "
                        "nucleotide path; it is not summed over alternatives."),
                      p("Reference ORF translation evidence comes from other experiments and does not establish "
                        "translation of this sample's mutant sequence. No unobserved RNA sequence is filled from reference."),
                      p("Every call, alternative, error, original support identity, input hash and parameter is retained "
                        "in the companion source exports and receipts. Reproduce with the same cohort manifest; "
                        "use --render-only to regenerate figures and PDFs from verified results.")])
        for index, source in enumerate(sources):
            block = [p("<b>" + escape(source.get("label", source["source_id"])) + "</b>", "Small"),
                     p(escape(json.dumps(source["provenance"], sort_keys=True)), "Small")]
            if index == 0:
                block.insert(0, Paragraph("Source provenance", styles["Heading2"]))
            story.append(KeepTogether(block))
        normal_scopes = {result['matched_normal_scope']['receipt_sha256']: result['matched_normal_scope']
                         for result in results if result.get('matched_normal_scope')}
        for digest, scope in normal_scopes.items():
            story.extend([Paragraph('Matched-normal evidence scope', styles['Heading2']),
                          p(escape(scope['interpretation']), 'Small'),
                          p('Assessment receipt: ' + escape(scope['receipt_input']) + ' [' + digest[:16] + ']. '
                            'Unavailable territory and absent variant-only calls do not establish reference. '
                            'Caller products from the same library are not independent specimens.', 'Small')])
            story.append(table(('Allele state', 'Nominations'),
                sorted(scope['allele_states'].items()), (410, 100)))
            selected = [row for row in scope['source_inventory'] if row['selected']]
            if selected:
                story.append(table(('Source / specimen / library', 'Caller / evidence', 'Confidence policy'),
                    [(row['source_id'] + '\n' + row['specimen_id'] + '\n' + row['library_id'],
                      row['caller'] + '\n' + row['evidence_kind'],
                      row['policy']['genotype_quality_field'] + ' >= ' + str(row['policy']['minimum_genotype_quality']) +
                      '; DP or block MIN_DP >= ' + str(row['policy']['minimum_depth']) +
                      ('; unfiltered accepted' if row['policy']['allow_unfiltered'] else '; PASS required'))
                     for row in selected], (230, 130, 150)))
            excluded = Counter(row['selection_reason'] for row in scope['source_inventory'] if not row['selected'])
            if excluded:
                story.append(table(('Excluded or unavailable product: reason', 'Files'), sorted(excluded.items()), (450, 60)))
            story.append(p(f"{scope['background_alleles']:,} genotype-indexed background alleles; "
                           f"query flank {scope['background_flank']:,} bases. The pinned normal receipt contains "
                           'every source ID, exact VCF sample, processing relationship, acquisition interval, '
                           'file/index digest, rejected record and background-allele provenance.', 'Small'))
        input_keys = set()
        for result in results:
            input_keys.update(result.get('identity', {}).get('inputs', {}))
        for entry in summary.get('additional', []):
            if entry['task'].get('sample_id') == sample_id and entry['status'] == 'recorded':
                input_keys.update(entry['identity']['inputs'])
                if entry['task'].get('locus_annotations'):
                    input_keys.add(entry['task']['locus_annotations'])
        if input_keys:
            story.append(Paragraph('Reference and analysis inputs', styles['Heading2']))
            story.append(p('SHA-256 prefixes identify the pinned files. Full digests, original source records '
                           'and paths are retained in the cohort manifest and per-source receipts.', 'Small'))
            story.append(table(('Input', 'Description', 'SHA-256 prefix'),
                [(key, manifest['files'][key].get('description', Path(manifest['files'][key]['path']).name),
                  manifest['files'][key]['sha256'][:16]) for key in sorted(input_keys)], (150, 260, 100)))
        if sample.get('source_inventory'):
            descriptor = manifest['files'][sample['source_inventory']]
            inventory_path = Path(summary['manifest']).parent / descriptor['path']
            if file_sha256(inventory_path) != descriptor['sha256']:
                raise ValueError('Source inventory digest mismatch')
            inventory = read_json(inventory_path)
            story.extend([Paragraph('Complete source inventory', styles['Heading2']),
                          p(escape(inventory['scope']), 'Small')])
            for kind in dict.fromkeys(row['kind'] for row in inventory['sources']):
                selected = [row for row in inventory['sources'] if row['kind'] == kind]
                story.append(Paragraph(escape(kind) + f' ({len(selected)} files)', styles['Heading2']))
                story.append(table(('Source / source identifier', 'Timepoint', 'Assessment state'),
                    [(row['label'] + ' [' + row['source_id'][:12] + ']', row.get('timepoint', ''), row['state'])
                     for row in selected], (250, 55, 205)))
            story.append(p('Original URLs, full source IDs, available checksums and snapshot identities are in the '
                           'pinned source-inventory JSON beside the cohort manifest.', 'Small'))
        doc = SimpleDocTemplate(str(report), pagesize=letter, leftMargin=42, rightMargin=42,
                                topMargin=42, bottomMargin=48, title=sample_id + " RNA evidence")
        doc.build(story, onFirstPage=footer, onLaterPages=footer)
        paths.append(report)
    write_json(output / "report_files.json", [dict(path=str(p), sha256=file_sha256(p)) for p in paths])
    return paths
