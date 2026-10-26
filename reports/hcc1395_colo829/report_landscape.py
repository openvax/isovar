"""Report-only landscape figures and source-qualified exemplar chapters."""

from collections import Counter, defaultdict
from math import log10

from reportlab.lib import colors
from reportlab.platypus import Flowable, KeepTogether, PageBreak, Spacer

from isovar.evidence_pdf import ProteinComparison
from isovar.evidence_report import variant_type
from reports.hcc1395_colo829.build_evidence_pdfs import (
    REPORT, NAVY, TEAL, GRAY, PALE, ORANGE, data_table, para, read_tsv, section,
)


class NeutralReferenceComparison(ProteinComparison):
    """Color only differing mutant residues; keep reference residues neutral."""

    def __init__(self, reference, mutant, mutation_start, mutation_end, focal=True):
        Flowable.__init__(self)
        self.reference, self.mutant = reference, mutant
        self.mutation_start, self.mutation_end = mutation_start, mutation_end
        self.focal = focal
        self.width, self.height = 512, 47

    def draw(self):
        c = self.canv
        c.setFillColor(PALE)
        c.roundRect(0, 0, self.width, self.height, 5, fill=1, stroke=0)
        for label, sequence, y in [('REF', self.reference, 28), ('RNA', self.mutant, 10)]:
            c.setFillColor(NAVY)
            c.setFont('Helvetica-Bold', 8)
            c.drawString(10, y, label)
            c.setFont('Courier-Bold', 9)
            for i, residue in enumerate(sequence):
                different = label == 'RNA' and (i >= len(self.reference) or residue != self.reference[i])
                focal = self.focal and self.mutation_start <= i < self.mutation_end
                c.setFillColor(ORANGE if different and focal else GRAY if different else NAVY)
                c.drawString(53 + i * 5.4, y, residue)


class SVRoutes(Flowable):
    """Parallel evidence channels, deliberately not a nested SV funnel."""

    def __init__(self, data):
        super().__init__()
        self.data, self.width, self.height = data, 512, 172

    def draw(self):
        c, d = self.canv, self.data
        c.setFont('Helvetica', 9)
        c.setFillColor(NAVY)
        covered = sum(int(r['end1_alignment_records']) > 0 or int(r['end2_alignment_records']) > 0 for r in d['sv'])
        c.drawString(0, 158, f"{len(d['sv']):,} nominated SV records; {covered:,} have RNA near an end")
        channels = [('Split alignments', 'Observed placed segments'),
                    ('Cross-end mates', 'Template linkage / context'),
                    ('Soft clips / unmapped', 'Unplaced sequence for assembly')]
        for i, (title, description) in enumerate(channels):
            x = i * 174
            c.setFillColor(PALE)
            c.roundRect(x, 96, 164, 48, 5, fill=1, stroke=0)
            c.setFillColor(NAVY)
            c.setFont('Helvetica-Bold', 9)
            c.drawString(x + 9, 124, title)
            c.setFont('Helvetica', 7.6)
            c.drawString(x + 9, 110, description)
            c.setStrokeColor(GRAY)
            c.line(x + 82, 94, 256, 76)
        c.setFillColor(TEAL)
        c.roundRect(55, 42, 402, 34, 5, fill=1, stroke=0)
        c.setFillColor(colors.white)
        c.setFont('Helvetica-Bold', 10)
        c.drawCentredString(256, 62, 'Reconstruct competing RNA paths; check original witnesses')
        c.setFont('Helvetica', 8)
        c.drawCentredString(256, 49, 'Separate junction, local branch, exon, frame and complete-sequence evidence')
        c.setFillColor(GRAY)
        c.setFont('Helvetica', 8)
        c.drawString(0, 24, 'The legacy split/mate screen is one route, not an eligibility gate.')
        c.drawString(0, 10, 'Soft-clip counts were not enumerated in the old atlas; this channel is not zero.')


class Landscape(Flowable):
    """Inventory composition and ALT/window-support landscape on log axes."""

    def __init__(self, data):
        super().__init__()
        self.data, self.width, self.height = data, 512, 385

    def draw(self):
        c, d = self.canv, self.data
        classes = Counter(variant_type(r) for r in d['calls'])
        total = len(d['calls'])
        c.setFillColor(NAVY)
        c.setFont('Helvetica-Bold', 10)
        c.drawString(0, 371, 'DNA inventory composition and original-RNA outcomes')
        c.setFont('Helvetica', 8)
        c.drawString(280, 351, 'DNA rows')
        c.drawString(360, 351, 'ALT calls')
        c.drawString(435, 351, 'Local context')
        for i, (kind, count) in enumerate(classes.items()):
            y = 331 - 23 * i
            subset = [r for r in d['calls'] if variant_type(r) == kind]
            alt = sum(r['rna_status'] == 'alternate_observed' for r in subset)
            local = sum(r['call_id'] in d['summary']['top'] for r in subset)
            c.setFillColor(NAVY)
            c.setFont('Helvetica', 9)
            c.drawString(0, y, kind)
            c.setFillColor(PALE)
            c.rect(78, y - 2, 180, 11, fill=1, stroke=0)
            c.setFillColor(TEAL)
            c.rect(78, y - 2, 180 * count / total, 11, fill=1, stroke=0)
            c.setFillColor(NAVY)
            for x, value in [(317, count), (398, alt), (480, local)]:
                c.drawRightString(x, y, f'{value:,}')
        left, bottom, width, height = 50, 35, 440, 185
        calls = [r for r in d['calls'] if int(r['alt_fragments'])]
        top = d['summary']['top']
        maximum = max([int(r['alt_fragments']) for r in calls] + [1])
        bound = max(1, int(log10(maximum)) + 1)
        project = lambda n: log10(1 + n) / log10(1 + 10 ** bound)
        c.setFillColor(NAVY)
        c.setFont('Helvetica-Bold', 9)
        c.drawString(0, 239, 'RNA allele support versus support for the top local protein window')
        for power in range(bound + 1):
            tick = 0 if power == 0 else 10 ** power
            px, py = left + width * project(tick), bottom + height * project(tick)
            c.setStrokeColor(PALE)
            c.line(px, bottom, px, bottom + height)
            c.line(left, py, left + width, py)
            c.setFillColor(GRAY)
            c.setFont('Helvetica', 7)
            c.drawCentredString(px, bottom - 12, str(tick))
            c.drawRightString(left - 5, py - 2, str(tick))
        points = []
        for row in calls:
            n = int(row['alt_fragments'])
            context = top.get(row['call_id'])
            k = int(context['supporting_fragments']) if context else 0
            frameshift = 'FrameShift' in d['annotations'][row['call_id']]['top_effect_type']
            color = ORANGE if frameshift else TEAL if row['class'] != 'SNV' else GRAY
            x, y = left + width * project(n), bottom + height * project(k)
            c.setFillColor(color)
            c.circle(x, y, 2.3 if frameshift else 1.6, fill=1, stroke=0)
            points.append((n, k, row, x, y))
        for index, (_, _, row, x, y) in enumerate(sorted(points, reverse=True, key=lambda x: x[:2])[:3]):
            gene = d['loci'][row['call_id']]['protein_genes'] or d['loci'][row['call_id']]['genes'] or row['call_id']
            c.setFillColor(NAVY)
            c.setFont('Helvetica', 7)
            c.drawRightString(min(490, x + 12), min(227, y + 10 + 9 * index), gene.split(';')[0])
        c.setFillColor(GRAY)
        c.setFont('Helvetica', 8)
        c.drawCentredString(left + width / 2, 6, 'ALT fragments (log1p); vertical axis: local-window fragments (log1p)')


def landscape_story(sample, d):
    story = section('Variant landscape: depth, coding recovery and unresolved categories')
    story.append(para('Every nominated small variant is counted above. The scatter shows ALT-observed '
                      'calls; points on the baseline have no local protein context under the recorded '
                      'settings. A context may be reference-identical. Separate cultures and assays '
                      'are never pooled.'))
    story.append(Landscape(d))
    story.append(para('<font color="#526373">&#9632;</font> SNV &nbsp; '
                      '<font color="#ba5a29">&#9632;</font> Frameshift &nbsp; '
                      '<font color="#0f706d">&#9632;</font> Other indel', 'FootSmall'))
    story.append(data_table(['Category', 'What this report can establish'], [
        ('Small variants', 'Original-RNA allele observations and bounded local coding windows'),
        ('SVs / fusions', 'Original joins and competing RNA paths; annotated coding frame only where resolved'),
        ('ncORFs', 'Reference coding-block overlaps and RNA ALT; mutant ncORF proteins remain unresolved'),
        ('Viral / other antigen sources', 'Not surveyed in these host-focused inputs; no absence claim'),
    ], [145, 367]))
    story.append(PageBreak())
    return story

def exemplar_story(sample, d, k=3):
    top = d['summary']['top']
    contexts = {r['call_id']: r for r in d['example_references']}
    supported = sorted((r for r in d['calls'] if int(r['alt_fragments'])),
                       key=lambda r: (-int(r['alt_fragments']), r['call_id']))
    groups = [
        ('Highest-support SNVs with local coding context',
         [r for r in supported if r['class'] == 'SNV' and r['call_id'] in top]),
        ('Top RNA-supported frameshift nominations',
         [r for r in supported if 'FrameShift' in d['annotations'][r['call_id']]['top_effect_type']]),
        ('Top other insertions, deletions and replacements',
         [r for r in supported if r['class'] != 'SNV' and
          'FrameShift' not in d['annotations'][r['call_id']]['top_effect_type']]),
    ]
    story = []
    for title, candidates in groups:
        story += section(title)
        story.append(para(f'Top {min(k, len(candidates))} of {len(candidates)} ALT-supported calls in this '
                          'category, ranked by source-library ALT fragments. Window support is reported '
                          'separately. A reference-identical sequence is a coding context, not a neoORF. '
                          'DNA effect labels and the displayed reference may use different transcript isoforms.', 'FootSmall'))
        if not candidates:
            story.append(para('No ALT-supported nomination in this category. This is not evidence of biological absence.'))
        for call in candidates[:k]:
            ident = call['call_id']
            locus = d['loci'][ident]
            genes = locus['protein_genes'] or locus['genes'] or 'No overlapping gene'
            annotation = d['annotations'][ident]
            block = [para(f"<b>{genes.replace(';', ', ')} | {ident}</b>"),
                     para(f"IGV (GRCh38): {locus['igv_loci']} | {annotation['top_effect_type']} "
                          f"{annotation['top_effect_description']} | {int(call['alt_fragments']):,} ALT fragments", 'FootSmall')]
            if ident in top:
                product = top[ident]
                sequence = product['protein_sequence']
                block.append(para(f"{int(product['supporting_fragments']):,} reported window fragments; "
                                  f"{len(sequence)} aa. Transcript(s): {product['compatible_transcript_ids']}", 'FootSmall'))
                if ident in contexts:
                    reference = contexts[ident]
                    block.append(NeutralReferenceComparison(reference['reference_context'], sequence,
                                 int(product['mutation_start_idx']), int(product['mutation_end_idx']),
                                 focal=annotation['top_effect_type'] != 'Silent' and
                                 reference['placement_method'] == 'unique_shared_protein_prefix'))
                    same = reference['reference_context'] == sequence
                    caveat = ('Reference-identical local context. ' if same else
                              'RNA differences have unclassified provenance; DNA-codon-anchored comparison only. '
                              if reference['placement_method'] != 'unique_shared_protein_prefix' else
                              'Differing focal RNA residues highlighted. ')
                    block.append(para(caveat
                                      + f"Reference: {reference['transcript_id']}, residue {reference['protein_start_1based']}. "
                                      'No unobserved sequence is filled from reference.', 'FootSmall'))
                else:
                    block.append(para(sequence, 'TinyMono'))
                    block.append(para('No unambiguous pinned reference-window placement; comparison remains unresolved.', 'FootSmall'))
            else:
                block.append(para('RNA ALT observed, but no local protein context was reconstructed. No protein sequence is inferred.', 'FootSmall'))
            block.append(Spacer(1, 10))
            story.append(KeepTogether(block))
        if sample != 'HCC1395' or title != 'Top RNA-supported frameshift nominations':
            story.append(PageBreak())
    return story


def sv_exemplars(sample, d, k=3):
    grouped = defaultdict(list)
    for row in read_tsv(REPORT / 'sv_transcript_evidence.tsv'):
        if row['sample'] == sample:
            grouped[row['call_id']].append(row)
    ranked = sorted(grouped.items(), key=lambda item: -max(int(r['direct_junction_fragments']) for r in item[1]))[:k]
    story = section('SV-derived RNA transcripts: junctions, branches and unresolved frames')
    story.append(para('These are sequence-supported RNA paths. A breakpoint-associated transcript is '
                      'not automatically a coding fusion or a DNA-causal product. Split placement, '
                      'soft clips, original query sequence and exon anchors are complementary evidence.'))
    for ident, rows in ranked:
        locus = d['loci'][ident]
        genes = locus['genes'] or 'No gene at nominated ends'
        best = max(rows, key=lambda r: (int(r['local_exact_fragments']), int(r['direct_junction_fragments'])))
        block = [para(f'<b>{genes.replace(";", ", ")} | {ident}</b>'),
                 para('IGV (GRCh38): ' + locus['igv_loci'], 'FootSmall'),
                 para(f"{best['direct_junction_fragments']} original junction fragments; "
                      f"{len(rows)} pinned example RNA paths; representative local-window support "
                      f"{best['local_exact_fragments']} fragments. These measures are not pooled.", 'FootSmall'),
                 para('RNA junction context (nucleotides; no protein frame assigned):', 'FootSmall'),
                 para(best['junction_left'] + '<font color="#ba5a29">' + best['junction_inserted'] + '</font>'
                      + best['junction_right'], 'TinyMono')]
        if 'SV3_DEL' in ident:
            block.append(para('The newer shared SV3 audit assesses all retained paths and original witnesses: '
                              '3 placed junction fragments and 5 exact junction-window fragments in the expanded '
                              'same-library source. Exon-linked sequence lies 280 nt downstream; no annotated '
                              'frame is resolved. Split/mate-only screening missed this evidence.', 'FootSmall'))
        else:
            block.append(para('The pinned exon-path audit retains ordinary-transcript and ambiguous-frame '
                              'interpretations. No event-linked protein is established by these paths.', 'FootSmall'))
        block.append(Spacer(1, 12))
        story.append(KeepTogether(block))
    story += section('Coding-fusion exemplar')
    if sample == 'HCC1395':
        story.append(para('<b>EIF3K-CYP39A1:</b> separate supplied-fusion input resolves a 531-aa coding '
                          'candidate with 149 Illumina junction fragments and independent ONT support. '
                          'The appended deep dive shows full parental reference alignments and exact '
                          'per-codon/assembly coverage. It is outside the benchmark-SV funnel.'))
    else:
        story.append(para('No event-linked annotated-frame protein was resolved in the new COLO829 '
                          'inventory replay. SV3 exploratory ORFs are not promoted to coding proteins.'))
    story.append(PageBreak())
    return story


def other_categories(sample, d, k=3):
    candidates = sorted((r for r in d['calls'] if int(r['alt_fragments']) and
                         d['annotations'][r['call_id']]['orfdata_coding_block_overlap'] == 'True'),
                        key=lambda r: -int(r['alt_fragments']))
    story = section('ncORFs and other antigen-source categories')
    story.append(para('Reference ncORF overlap nominates an alternative coding model; it does not '
                      'establish that the ORF is translated, mutated in-frame, or covered over its '
                      'complete coding sequence. No mutant ncORF protein was reconstructed here.'))
    rows = []
    for call in candidates[:k]:
        annotation, locus = d['annotations'][call['call_id']], d['loci'][call['call_id']]
        rows.append((f"<b>{locus['genes'].replace(';', ', ') or 'No overlapping gene'}</b><br/>"
                     f"{call['call_id']}<br/>{locus['igv_loci']}",
                     f"{int(call['alt_fragments']):,}", annotation['orfdata_biotypes'],
                     'Unresolved mutant ncORF'))
    story.append(data_table(['Top overlap / IGV (GRCh38)', 'ALT fragments', 'Reference ORF class', 'Protein status'],
                            rows, [220, 65, 105, 122]))
    story.append(para(f'{len(candidates)} ncORF-overlap calls have ALT RNA; the table shows up to {k}. '
                      'Full ORF model identifiers remain in reference_coding_orfdata_annotations.tsv.gz. '
                      'A canonical protein context at the same locus does not validate the ncORF.', 'FootSmall'))
    story += section('Assessment coverage - unknown is not absent')
    story.append(data_table(['Category', 'State', 'What is missing'], [
        ('Viral RNA / viral ORFs', 'Not assessed', 'No viral-aligned RNA/reference panel in this report'),
        ('General cryptic initiation / retained introns', 'Not comprehensively assessed', 'Transcript discovery and initiation evidence beyond selected paths'),
        ('Targeted cancer/testis antigen isoforms', 'Not assessed', 'Source-qualified isoform coverage audit'),
        ('Germline / other somatic phasing', 'Not classified beyond available paths', 'Matched genotype and explicit linkage are required for colored attribution'),
        ('Translated protein / HLA presentation', 'Not measured by these RNA data', 'Independent proteomic / immunopeptidomic evidence'),
    ], [165, 135, 212]))
    story.append(PageBreak())
    return story
