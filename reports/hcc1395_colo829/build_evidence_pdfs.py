"""Build evidence PDFs for HCC1395 and COLO829 from pinned Isovar tables."""

import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path

from reportlab import rl_config
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import (
    HRFlowable, PageBreak, Paragraph, SimpleDocTemplate,
    Spacer, Table, TableStyle,
)

from isovar import __version__
from isovar.evidence_pdf import (
    EvidenceFunnel, altered_aa_figures, read_support_figures,
    support_by_class_figures,
)
from isovar.evidence_report import read_tsv as read_evidence_tsv, summarize


rl_config.invariant = 1
ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / 'reports/hcc1395_colo829'
ANNOTATIONS = REPORT / 'reference_coding_orfdata_annotations.tsv.gz'
SELF_SCORES = REPORT / 'reference_self_9mer_scores.tsv.gz'
SV_EXONS = REPORT / 'sv_exon_annotations.tsv.gz'
SV_JUNCTION_AUDIT = REPORT / 'sv_junction_audit.tsv'
EXAMPLE_REFERENCES = REPORT / 'example_reference_contexts.tsv'
LOCUS_ANNOTATIONS = REPORT / 'report_locus_annotations.tsv'
SELF_FASTA_SHA256 = 'd873c199d5dbf20869db3f2cbd6ac726259215d91c2782f9f1ddcad4bb117cf9'
OUT = ROOT / 'output/pdf'
NAVY = colors.HexColor('#18334b')
TEAL = colors.HexColor('#0f706d')
GRAY = colors.HexColor('#526373')
LIGHT = colors.HexColor('#eaf1f4')
PALE = colors.HexColor('#f5f8fa')
ORANGE = colors.HexColor('#ba5a29')
PAGE_W, PAGE_H = letter

styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name='ReportTitle', parent=styles['Title'],
                          fontName='Helvetica-Bold', fontSize=20,
                          leading=24, textColor=NAVY, spaceAfter=8))
styles.add(ParagraphStyle(name='ReportDeck', parent=styles['Normal'],
                          fontName='Helvetica', fontSize=9.5, leading=14,
                          textColor=GRAY, spaceAfter=12))
styles.add(ParagraphStyle(name='Section', parent=styles['Heading2'],
                          fontName='Helvetica-Bold', fontSize=12.5,
                          leading=16, textColor=TEAL, spaceBefore=13,
                          spaceAfter=6))
styles.add(ParagraphStyle(name='BodySmall', parent=styles['Normal'],
                          fontName='Helvetica', fontSize=8.6, leading=12.3,
                          textColor=NAVY, spaceAfter=6))
styles.add(ParagraphStyle(name='FootSmall', parent=styles['Normal'],
                          fontName='Helvetica', fontSize=7.4, leading=10.3,
                          textColor=GRAY, spaceAfter=5))
styles.add(ParagraphStyle(name='Cell', parent=styles['Normal'],
                          fontName='Helvetica', fontSize=8.2, leading=10.6,
                          textColor=NAVY))
styles.add(ParagraphStyle(name='CellHead', parent=styles['Normal'],
                          fontName='Helvetica-Bold', fontSize=8.2, leading=10.6,
                          textColor=colors.white))
styles.add(ParagraphStyle(name='TinyMono', parent=styles['Normal'],
                          fontName='Courier', fontSize=6.8, leading=9,
                          textColor=NAVY))


def para(text, style='BodySmall'):
    return Paragraph(text, styles[style])


def data_table(headers, rows, widths, *, numeric_columns=()):
    data = [[para(h, 'CellHead') for h in headers]]
    for row in rows:
        data.append([para(str(value), 'Cell') for value in row])
    table = Table(data, colWidths=widths, repeatRows=1, hAlign='LEFT')
    commands = [
        ('BACKGROUND', (0, 0), (-1, 0), NAVY),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, PALE]),
        ('LINEBELOW', (0, 0), (-1, 0), 0.8, TEAL),
        ('GRID', (0, 1), (-1, -1), 0.25, colors.HexColor('#d9e4e9')),
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ('LEFTPADDING', (0, 0), (-1, -1), 7),
        ('RIGHTPADDING', (0, 0), (-1, -1), 7),
        ('TOPPADDING', (0, 0), (-1, -1), 5),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
    ]
    for col in numeric_columns:
        commands.append(('ALIGN', (col, 1), (col, -1), 'RIGHT'))
    table.setStyle(TableStyle(commands))
    return table


def read_gzip(path):
    with gzip.open(path, 'rt', newline='') as source:
        return list(csv.DictReader(source, delimiter='\t'))


def read_tsv(path):
    with path.open(newline='') as source:
        return list(csv.DictReader(source, delimiter='\t'))


def collect(sample, annotations_path=ANNOTATIONS):
    prefix = 'hcc1395_direct_rna' if sample == 'HCC1395' else 'colo829'
    calls = read_gzip(REPORT / f'{prefix}_atlas_calls.tsv.gz')
    proteins = read_gzip(REPORT / f'{prefix}_atlas_proteins.tsv.gz')
    sv_file = ('hcc1395_direct_rna_sv_screen.tsv' if sample == 'HCC1395'
               else 'colo829_sv_rna_screen.tsv')
    sv = read_tsv(REPORT / sv_file)
    sv_audit = [row for row in read_tsv(SV_JUNCTION_AUDIT) if row['sample'] == sample]
    sv_ids = {row['call_id'] for row in sv}
    if not sv_audit or len({row['call_id'] for row in sv_audit}) != len(sv_audit) or any(
            row['call_id'] not in sv_ids for row in sv_audit):
        raise ValueError(f'Invalid selected SV junction audit for {sample}')
    sv_annotations = {(row['dna_source'], row['call_id']): row
                      for row in read_gzip(SV_EXONS) if row['sample'] == sample}
    if len(sv_annotations) != len(sv):
        raise ValueError(f'Incomplete SV exon annotations for {sample}')
    for row in sv:
        key = row['dna_source'], row['call_id']
        annotation = sv_annotations[key]
        if any(row[field] != annotation[field]
               for field in ('class', 'chrom1', 'pos1', 'chrom2', 'pos2')):
            raise ValueError(f'SV annotation mismatch: {key}')
        row['exon_overlap'] = annotation['exon_overlap']
        row['exon_rule'] = annotation['exon_rule']
    annotations = {row['call_id']: row for row in read_gzip(annotations_path)
                   if row['sample'] == sample}
    assert len(calls) == len(annotations)
    assert {c['call_id'] for c in calls} == set(annotations)
    assert all(not row['annotation_error'] for row in annotations.values())
    summary = summarize(calls, annotations, proteins, sv)
    counts = Counter()
    rna = Counter(row['rna_status'] for row in calls)
    classes = Counter(row['class'] for row in calls)
    filters = Counter(row['dna_filter'] for row in calls)
    sv_classes = Counter(row['class'] for row in sv)
    sv_status = Counter(row['rna_status'] for row in sv)
    for call in calls:
        annotation = annotations[call['call_id']]
        coding = annotation['reference_coding_effect'] == 'True'
        ncorf = annotation['orfdata_coding_block_overlap'] == 'True'
        if coding:
            counts['coding'] += 1
        if ncorf:
            counts['ncorf'] += 1
        if coding or ncorf:
            counts['coding_union'] += 1
        if coding and ncorf:
            counts['coding_both'] += 1
        if call['rna_status'] == 'alternate_observed':
            if coding:
                counts['coding_alt'] += 1
            if ncorf:
                counts['ncorf_alt'] += 1
        if int(call['protein_count']):
            counts['protein_calls'] += 1
            assert coding
        counts['alt_reads'] += int(call['alt_reads'])
        counts['alt_fragments'] += int(call['alt_fragments'])
    counts['protein_contexts'] = len(proteins)
    counts['sv_split_records'] = sv_status['split_alignment_candidate_unresolved']
    counts['sv_split_witnesses'] = sum(int(row['observed_split_alignment_candidates']) for row in sv)
    counts['sv_exonic'] = sum(row['exon_overlap'] == 'True' for row in sv)
    counts['sv_exon_unknown'] = sum(row['exon_overlap'] == 'Unknown' for row in sv)
    exonic_sv = [row for row in sv if row['exon_overlap'] == 'True']
    counts['sv_exonic_covered'] = sum(
        int(row['end1_alignment_records']) > 0 or int(row['end2_alignment_records']) > 0
        for row in exonic_sv)
    counts['sv_exonic_candidate'] = sum(
        int(row['observed_split_alignment_candidates']) > 0 or
        int(row.get('cross_end_mate_pair_candidates') or 0) > 0
        for row in exonic_sv)
    exonic_ids = {row['call_id'] for row in exonic_sv}
    counts['sv_exonic_junction'] = sum(
        row['call_id'] in exonic_ids and int(row['junction_fragments_all_mapq']) > 0
        for row in sv_audit)
    counts['sv_exonic_product'] = sum(
        row['protein_status'] == 'called' for row in exonic_sv)
    for key, value in summary['counts'].items():
        if key in counts and counts[key] != value:
            raise ValueError(f'Inconsistent {key} count for {sample}')
        counts[key] = value
    self_scores = {}
    if SELF_SCORES.exists():
        self_scores = {row['call_id']: row for row in read_evidence_tsv(SELF_SCORES)
                       if row['sample'] == sample}
        if set(self_scores) != set(summary['top']) or any(
                row['reference_fasta_sha256'] != SELF_FASTA_SHA256
                for row in self_scores.values()):
            raise ValueError(f'Incomplete or incompatible reference-self scores for {sample}')
    examples = [row for row in read_tsv(EXAMPLE_REFERENCES)
                if row['sample'] == sample]
    if not examples or any(row['call_id'] not in summary['top'] or
                           row['ensembl_release'] != '87' for row in examples):
        raise ValueError(f'Incomplete example reference contexts for {sample}')
    return {
        'calls': calls, 'proteins': proteins, 'sv': sv, 'sv_audit': sv_audit,
        'counts': counts,
        'rna': rna, 'classes': classes, 'filters': filters,
        'sv_classes': sv_classes, 'sv_status': sv_status,
        'summary': summary, 'annotations': annotations, 'self_scores': self_scores,
        'example_references': examples,
        'loci': {r['call_id']: r for r in read_tsv(LOCUS_ANNOTATIONS) if r['sample'] == sample},
    }


def page_footer(canvas, doc, sample):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor('#c6d5dc'))
    canvas.line(50, PAGE_H - 39, PAGE_W - 50, PAGE_H - 39)
    canvas.setFont('Helvetica-Bold', 7.5)
    canvas.setFillColor(TEAL)
    canvas.drawString(50, PAGE_H - 31, f'ISOVAR {__version__}  /  {sample} EVIDENCE AUDIT')
    canvas.setFont('Helvetica', 7.5)
    canvas.setFillColor(GRAY)
    canvas.drawRightString(PAGE_W - 50, 32, f'Updated 2026-10-26  |  page {doc.page}')
    canvas.restoreState()


def section(title, body=None):
    result = [para(title, 'Section')]
    if body:
        result.append(para(body))
    return result


def main_story(sample, d):
    from .report_landscape import SVRoutes
    c = d['counts']
    small = len(d['calls'])
    sv_total = len(d['sv'])
    story = [
        Spacer(1, 4),
        para(f'{sample}: somatic mutation and RNA evidence', 'ReportTitle'),
        para('Whole source-call audit | matched-normal DNA baseline | original tumor RNA | '
             'pinned reference ORFs. Counts below describe the specified public call sets and '
             'RNA run, not every biological mutation or a validated full-length isoform.', 'ReportDeck'),
        HRFlowable(width='100%', thickness=1.5, color=TEAL),
    ]
    story += section('Small variants: DNA to RNA and local protein')
    story.append(EvidenceFunnel(c))
    story += section('Structural variants: complementary routes to RNA reconstruction')
    story.append(SVRoutes(d))
    story.append(para(
        'Breakend-window coverage is locus coverage. Split reads, cross-end mates and soft clips '
        'are complementary discovery evidence; none alone resolves a coding frame. Exon overlap '
        'is an annotation, not an eligibility filter. The older exon-restricted split/mate screen '
        'and selected junction replay are documented later without treating their negatives as '
        'absence. The newer COLO829 workflow uses all retained paths and exon-anchored clips.', 'FootSmall'))
    story.append(PageBreak())
    story += section('At a glance')
    story.append(data_table(
        ['Evidence layer', 'Count', 'Interpretation'],
        [
            ('Ingested small variants', f'{small:,}',
             f"{d['classes']['SNV']:,} SNVs; {d['classes']['indel']:,} indels"),
            ('Ingested SV source records', f'{sv_total:,}',
             'DNA breakpoint records/events, counted in their source sets'),
            ('SVs involving an annotated exon', f"{c['sv_exonic']:,}",
             'Span for DEL/DUP/INV; nominated breakpoint or insertion flanks otherwise'),
            ('Any annotated exon involved', f"{c['exonic']:,}",
             'At least one changed base or insertion flank overlaps an Ensembl 87 exon'),
            ('Ensembl 87 coding effect', f"{c['coding']:,}",
             'Any annotated protein-coding transcript, reference model'),
            ('orfdata ncORF coding-block overlap', f"{c['ncorf']:,}",
             'Candidate reference ORF overlap; distinct from mutant translation'),
            ('Original RNA ALT-observed calls', f"{d['rna']['alternate_observed']:,}",
             f"{c['alt_reads']:,} ALT reads; {c['alt_fragments']:,} ALT fragments"),
            ('Atlas local protein contexts', f"{c['protein_calls']:,}",
             f"{c['protein_contexts']:,} ranked contexts; new ncORF pass reported separately"),
        ], [179, 68, 265], numeric_columns=(1,)))
    return story


def coding_story(sample, d):
    c = d['counts']
    small = len(d['calls'])
    story = []
    story += section('Coding annotation and its denominator')
    story.append(para(
        f"Of {small:,} DNA small-variant records, {c['coding']:,} change the coding sequence "
        f"of at least one Ensembl 87 transcript. The pinned orfdata 0.34.1 GRCh38 "
        f"catalog (7,264 reference ncORF models; release 2026-10-11.1) has exact coding-block "
        f"overlaps at {c['ncorf']:,} calls. The groups overlap at "
        f"{c['coding_both']:,} calls, so the "
        f"union is {c['coding_union']:,}. These are annotation opportunities. The catalog "
        'does not assess a cell-line-specific mutant ORF, RNA expression, or peptide evidence.'))
    story.append(data_table(
        ['Reference annotation', 'All DNA calls', 'With original RNA ALT'],
        [
            ('Ensembl coding effect', f"{c['coding']:,}", f"{c['coding_alt']:,}"),
            ('orfdata coding-block overlap', f"{c['ncorf']:,}", f"{c['ncorf_alt']:,}"),
            ('Union of either category', f"{c['coding_union']:,}", 'see per-call table'),
        ], [234, 131, 147], numeric_columns=(1, 2)))
    story.append(para(
        'SNVs query their changed reference base. Deletions query removed reference bases; '
        'insertions use the VCF anchor as a candidate overlap. The per-call table preserves '
        'that rule. A coding-block hit does not prove the inserted sequence is in frame.',
        'FootSmall'))
    story.append(para(
        f"Any annotated exon is involved in {c['exonic']:,} calls ({100 * c['exonic'] / small:.1f}% "
        'of the small-variant source set). This includes UTR and noncoding-transcript exons. '
        'Deletions query every removed reference base, including exon crossings with intronic '
        'endpoints; insertions query both neighboring reference bases as candidate exon involvement.',
        'FootSmall'))
    return story


def read_support_story(d):
    story = section('RNA alternate-read support by variant class')
    for figure in read_support_figures(d['summary']['reads_by_class']):
        story.append(figure)
    story.append(para(
        'The positive-support panel excludes zero-ALT records. Bars stack SNV, '
        'insertion, deletion and other-replacement record counts. Matching '
        'colored columns print each class count directly. These are variant '
        'counts, not read counts.', 'FootSmall'))
    story.append(PageBreak())
    return story


def exonic_read_support_story(d):
    story = section('RNA ALT reads: exon-involved small variants')
    for figure in read_support_figures(d['summary']['exonic_reads_by_class'],
                                       scope='exon-involved variant'):
        story.append(figure)
    story.append(para(
        'The input is the exon-involved subset of the small-variant source calls. '
        'The first panel retains zero-ALT calls; the second conditions on positive '
        'ALT reads. Insertion flanks are candidate exon involvement. SVs have a '
        'separate plot because a split alignment is not an ALT-spanning read.',
        'FootSmall'))
    story.append(PageBreak())
    return story


def sv_exon_story(d):
    sv = d['sv']
    included = [row for row in sv if row['exon_overlap'] == 'True']
    values = {}
    for row in included:
        values.setdefault(row['class'], []).append(
            int(row['observed_split_alignment_candidates']) +
            int(row.get('cross_end_mate_pair_candidates') or 0))
    split_records = sum(int(row['observed_split_alignment_candidates']) > 0
                        for row in included)
    mate_records = sum(int(row.get('cross_end_mate_pair_candidates') or 0) > 0
                       for row in included)
    story = section('Exon-involved nominated SVs: RNA breakpoint screen')
    story.append(para(
        f"{len(included):,} of {len(sv):,} nominated SV source records involve an Ensembl 87 "
        'exon. Intrachromosomal DEL/DUP/INV records query their full nominated '
        'span; INS queries insertion flanks; TRA/BND query nominated breakpoints. '
        f"{d['counts']['sv_exon_unknown']:,} records have insufficient breakpoint "
        'coordinates to rule out exon involvement and are excluded from this plot.'))
    story.append(data_table(
        ['RNA screen', 'Exonic SV records', 'Candidate reads or pairs'],
        [
            ('Same-read split alignments', f'{split_records:,}',
             f"{sum(int(row['observed_split_alignment_candidates']) for row in included):,}"),
            ('Cross-end primary mates', f'{mate_records:,}',
             f"{sum(int(row.get('cross_end_mate_pair_candidates') or 0) for row in included):,}"),
        ], [230, 130, 152], numeric_columns=(1, 2)))
    for figure in support_by_class_figures(
            values, ('DEL', 'INS', 'DUP', 'INV', 'TRA', 'BND'),
            title='RNA link candidates per exon-involved SV'):
        story.append(figure)
    story.append(para(
        'Bars count SV source records by reported class, binned on split-read '
        'plus cross-end mate candidate counts. Same-read supplementary alignments '
        'and opposite primary mates are different kinds of screening evidence; '
        'neither validates a junction sequence or translated product. '
        'The screen uses 500-base windows and mapping quality at least 20. '
        'The COLO829 benchmark call sets remain separate in the underlying table.',
        'FootSmall'))
    story.append(PageBreak())
    return story


def sv_junction_story(d):
    story = section('Original-read replay: junction versus RNA path')
    story.append(para(
        'This selected-event replay uses oriented DNA adjacencies and actual RNA '
        'alignment records, including soft-clipped bases and low-MAPQ pieces. '
        'Direct junction fragments count reads crossing the same RNA join, '
        'even when their wider transcript sequences differ. MAPQ 20 is a '
        'sensitivity subset, not a gate on whether the junction sequence was '
        'observed. Isovar extends these joins through overlapping original '
        'read sequences; one read need not span the assembled context. Context '
        'length is a bounded assembly result, not an isoform-support count.'))
    rows = []
    for row in d['sv_audit']:
        label = row['call_id'].replace('COLO829_SOMATIC_', '')
        parts = label.split('_')
        if len(parts) == 5 and parts[0] == parts[2] and parts[1].isdigit() and parts[3].isdigit():
            label = f'{parts[0]}:{int(parts[1]):,}-{int(parts[3]):,}'
        locus = d['loci'][row['call_id']]
        label += '<br/><b>' + (locus['genes'].replace(';', ', ') or 'No overlapping gene') + '</b>'
        label += '<br/>' + locus['igv_loci'].replace(' ', '<br/>')
        rows.append((label, row['class'], row['junction_fragments_all_mapq'],
                     row['junction_fragments_mapq20_subset'],
                     row['max_assembled_context_bases']))
    story.append(data_table(
        ['Nominated DNA record', 'Class', 'RNA join, all MAPQ',
         'Join, MAPQ 20 subset', 'Longest context (bp)'],
        rows, [203, 49, 83, 88, 89], numeric_columns=(2, 3, 4)))
    story.append(Spacer(1, 8))
    story.append(para('Gene names overlap the source breakends (Ensembl 87); they do not establish a '
                      'coding fusion. IGV loci are one-based inclusive, GRCh38, with no digit separators.', 'FootSmall'))
    if d['sv_audit'][0]['sample'] == 'HCC1395':
        story.append(para(
            'The chr9 duplication has seven direct junction fragments; three '
            'event paths survive assembly and the search reaches its ten-path '
            'cap. The chr13 deletion has four direct fragments, including one '
            'lost at MAPQ 20; chr21 has two. Overlapping reads extend each '
            'candidate context; their local branch support still needs '
            'evaluation. The three chr1 '
            'screen rows reuse one complex read and were not oriented here.'))
    else:
        story.append(para(
            'The chr1 SV3 deletion has six opposite-mate candidate pairs and '
            'three original fragments with an identical split junction sequence '
            'containing 19 unplaced bases. Its linked alignment pieces span '
            'MAPQ 1-52, including low-MAPQ supplementary pieces; the MAPQ 20 '
            'subset loses all three joins. Assembly extends this sequence to '
            '1,213 bases through overlapping reads and reaches an extension '
            'limit. SV42 and SV49 have one mate candidate each '
            'but no directly observed join. SV3 is outside the exon-involved '
            'funnel above.'))
    story.append(para(
        'The HCC1395 replay BAM is a source-derived 309-record regional '
        'fixture; COLO829 uses the 129,534-record selected-site BAM. Both '
        'were extracted after aligning the complete original RNA run. '
        'No transcript models were supplied in this geometry audit, so '
        'frame and SV-derived protein status remain unresolved. Other '
        'splicing and DNA-to-RNA causation have not been excluded.', 'FootSmall'))
    story.append(PageBreak())
    return story


def evidence_story(sample, d):
    c, rna, sv_status = d['counts'], d['rna'], d['sv_status']
    run = 'SRR23456709, ONT direct RNA' if sample == 'HCC1395' else 'SRR8615617, Illumina paired RNA'
    story = section('Original RNA and local protein evidence',
                    f'RNA assay: {run}. Each nominated small variant has a status from this run. '
                    'An ALT read is a sequence observation; it is not a full-length mutant isoform.')
    status_order = (
        'alternate_observed', 'reference_only', 'reference_and_other_only',
        'other_allele_only', 'alignment_without_allele_call',
        'no_alignment_overlap',
    )
    story.append(data_table(['RNA allele status', 'Small-variant records'],
                            [(key.replace('_', ' '), f'{rna[key]:,}') for key in status_order],
                            [340, 172], numeric_columns=(1,)))
    story.append(para(
        f"The {rna['alternate_observed']:,} ALT-observed calls contribute "
        f"{c['alt_reads']:,} allele-read and {c['alt_fragments']:,} fragment "
        "observations. Isovar resolves at least one "
        f"read-compatible <b>local</b> protein context for {c['protein_calls']:,} "
        f"calls ({c['protein_contexts']:,} ranked rows). The remaining "
        f"{rna['alternate_observed'] - c['protein_calls']:,} ALT-observed calls have "
        'no local translation under the stated method. No-overlap and uncallable statuses '
        'cannot establish absence of RNA expression.'))
    if sample == 'HCC1395':
        story.append(para(
            'OSTC F9L and PRDX5 F157L have original-read local contexts. An independent '
            'proteomic study observed the mutant residues, not these entire contexts. '
            'The bounded ONT cDNA fixture is excluded from whole-run counts.', 'FootSmall'))
    else:
        story.append(para(
            'At CDKN2A chr9:21971154:CCG&gt;C, all 6,984 ALT allele reads were evaluated. '
            'The top local frameshift context is <b>MMMGSARVAELLLLHGGAQLRRPR</b> '
            'with 5,338 supporting reads. Its 1,416 ranked local contexts reflect '
            'read-spanning alternatives, not full-length isoforms. TAGET polished '
            'representatives are an independent culture and assay; their cluster '
            'counts do not measure mutant allele support.', 'FootSmall'))
    story += section('Structural-variant RNA screen')
    status_names = {
        'no_breakend_window_alignment': 'No breakend-window alignment',
        'one_end_covered_no_split_candidate': 'One end covered; no split candidate',
        'both_ends_covered_no_split_candidate': 'Both ends covered; no split candidate',
        'split_alignment_candidate_unresolved': 'Split candidate; path unresolved',
        'cross_end_pair_context_unresolved': 'Cross-end mate context; path unresolved',
        'unresolved_partner_coordinate': 'Partner coordinate unresolved',
    }
    story.append(data_table(
        ['Breakend screen outcome', 'Source records'],
        [(status_names.get(key, key), f'{sv_status[key]:,}')
         for key in status_names if sv_status[key]],
        [340, 172], numeric_columns=(1,)))
    if sample == 'HCC1395':
        story.append(para(
            'Six Talsania records have same-read split-alignment candidates (17 '
            'record-level observations). Published POU2F1-NME7 has one candidate '
            'read. Junction bases, orientation, splice choice and ORF remain '
            'unresolved; no fusion protein is called.', 'FootSmall'))
    else:
        story.append(para('The two SV benchmarks remain separate; no junction '
                          'or SV-derived protein is called.', 'FootSmall'))
    story.append(PageBreak())
    return story


def altered_aa_story(d):
    story = section('Isovar-marked amino-acid span')
    story.append(para(
        'For each RNA ALT-observed call, count amino-acid positions in the '
        'mutation interval of its top-ranked local protein context. This is '
        'not an independently verified count of changed residues: synonymous '
        'DNA calls can produce a reference-identical amino acid window. '
        'Calls with no local context are shown '
        'explicitly. The second panel omits the one-residue category to expose '
        'the longer tail; its denominator is smaller. A zero-length interval can '
        'represent a deletion junction. These windows do not determine a complete '
        'frameshift tail or a full-length isoform.'))
    for figure in altered_aa_figures(d['summary']['aa_counts']):
        story.append(figure)
    story.append(PageBreak())
    return story


def provenance_story(sample, d):
    story = section('Sources, processing and limits')
    story.append(para('All RNA libraries, DNA nominations, polished inputs and references across the '
                      'report set are listed in <link href="https://github.com/openvax/isovar/blob/master/'
                      'reports/hcc1395_colo829/SEQUENCING_SOURCES.md">SEQUENCING_SOURCES.md</link>. '
                      'Focused CDKN2A and EIF3K-CYP39A1 reports show highlighted reference alignments '
                      'and original-RNA support separately from these atlas summaries.'))
    if sample == 'HCC1395':
        sources = [
            ('Matched-normal small variants', 'SEQC2 v1.2.1 SNV + indel, GRCh38.d1.vd1; HCC1395BL normal'),
            ('SV benchmark', 'Talsania et al. 2022 Table S4, HCC1395 vs HCC1395BL'),
            ('Original RNA', 'PRJNA934933 / SRR23456709; 1,049,452 ONT direct RNA reads'),
            ('Alignment', 'minimap2 2.26 -ax splice -uf -k14; NCBI GRCh38 no-alt'),
            ('ncORF witness realignment', '36 original ALT reads; GRCh38 full analysis set + decoy/HLA; same library'),
        ]
        hashes = [
            ('SEQC2 SNV VCF gzip', '0c49f7e1519fc7860b773f79bdcf07453fe4def8303707f34542471d92494c9d'),
            ('SEQC2 indel VCF gzip', 'db7dc8be2ec04bab830ff1411146c7d538c217ff684e874b85d98a3dbbda91cf'),
            ('Talsania Table S4 xlsx', '11e8d6393b83f0e0b7fdb0b13a57376663fc9662116d78c89a0d5fa530df096b'),
            ('SRR23456709 SRA archive', 'd69db8e48bbb20665c3e03a1fa35742be818a0127994e015cb3a560afae39ac8'),
        ]
        caveat = (
            'ONT raw-read errors can create other-base calls and competing local contexts. '
            'The RNA and DNA libraries come from separate studies. A 500-base breakend '
            'window and split-alignment candidate are screening evidence only. '
            'Selected junctions were resolved with original reads, but complete '
            'isoform assignment and SV protein reconstruction remain unresolved.'
        )
    else:
        sources = [
            ('Matched-normal small variants', 'SMaHT COLO829/COLO829BL, GRCh38 no-alt'),
            ('SV benchmarks', 'Paulin Paz 2024 v2 (54 records) and Valle-Inclan 2022 lifted (134 records)'),
            ('Original RNA', 'SRR8615617; 115,269,780 paired Illumina spots / 230,539,560 reads'),
            ('Alignment', 'minimap2 2.26 -ax sr --frag=yes; NCBI GRCh38 no-alt'),
        ]
        hashes = [
            ('SMaHT SNV VCF gzip', 'e57f77c8bc280bbad2255856c4062348f92f25da38726f367322c52736578f31'),
            ('SMaHT indel VCF gzip', 'e8d0fb67925a2365b485cd65feab8dfdc12bbe54c0bc0bdec44c1846b518b685'),
            ('SV v2 GRCh38 VCF', '8ad075bf43aa71d89c591396569e813959ada44a30f171982c2d299411290bac'),
            ('Earlier SV lifted VCF', 'a8f38454367cf8837a7dba1d29c9bc338a8e62cadfee6d2489931f8105f0dcdc'),
            ('SRR8615617 SRA archive', '9d89bb4f5a74fba7a8ba9495fd44046c35f4abb60d10c7bd9261f07569520b4f'),
        ]
        caveat = (
            'The SMaHT DNA, CCLE RNA and TAGET polished isoforms are separate cultures. '
            'The short-read genomic alignment is not splice-aware; junction-crossing '
            'reads may be lost. Three sequence-supported SV3 junction fragments '
            'have low-MAPQ linking pieces; placement, splice alternatives, '
            'full-length isoform linkage and SV-derived proteins remain unresolved.'
        )
    story.append(data_table(['Input', 'Pinned scope'], sources, [160, 352]))
    story += section('Reference model and interpretation')
    story.append(para(
        'Small-variant coding effects use Varcode with Ensembl release 87 on GRCh38. '
        'The human orfdata 0.34.1 catalog uses Ensembl 101 GRCh38 reference ORFs '
        '(7,264 models, dataset 2026-10-11.1). The new ORFdata 0.39.1 replay validates '
        'nine complete catalogue ORF paths/proteins against Ensembl 87, then reconstructs '
        'their frames on mutation-supporting RNA through the shared library API. '
        'The ncORF chapter separates default and single-read results. Genomic transcript overlap is not '
        'RNA-compatible isoform support. Original-read allele counts and local '
        'protein contexts are separate observations.'))
    story += section('SHA-256 input pins')
    for name, digest in hashes:
        story.append(para(f'<b>{name}</b>', 'FootSmall'))
        story.append(para(digest, 'TinyMono'))
    story.append(para('<b>GRCh38 no-alt FASTA (uncompressed)</b>', 'FootSmall'))
    story.append(para('9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7',
                      'TinyMono'))
    if sample == 'HCC1395':
        story.append(para('<b>ncORF realignment GRCh38 full analysis + decoy/HLA FASTA</b>', 'FootSmall'))
        story.append(para('3b103f4742abfd54938fb0333e19ad067635c8eb86f1dbf0ce44b165c4292b50', 'TinyMono'))
    story.append(para('<b>orfdata bundled GRCh38 catalog</b>', 'FootSmall'))
    story.append(para('06bc0bf3f81c25a7356f71a38581ff5704c05dcb75752454bc1ca75075ab8c09',
                      'TinyMono'))
    story.append(para('<b>Ensembl 87 GRCh38 protein FASTA gzip (self screen)</b>', 'FootSmall'))
    story.append(para(SELF_FASTA_SHA256, 'TinyMono'))
    story += section('What remains unresolved')
    story.append(para(caveat))
    story.append(para(
        'Machine-readable per-call annotation, allele, protein, and SV tables '
        'and complete reconstruction commands are published with the Isovar '
        'HCC1395/COLO829 reports. These PDFs summarize selected public source '
        'call sets, not a complete catalogue of every cell-line mutation.',
        'FootSmall'))
    return story


def build_one(sample, output_dir=OUT, annotations_path=ANNOTATIONS):
    from .report_landscape import landscape_story, exemplar_story, sv_exemplars, other_categories
    from .report_ddx17 import ddx17_story
    from pypdf import PdfWriter
    data = collect(sample, annotations_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f'{sample}_somatic_RNA_evidence.pdf'
    doc = SimpleDocTemplate(
        str(path), pagesize=letter, leftMargin=50, rightMargin=50,
        topMargin=52, bottomMargin=51, title=f'{sample} somatic mutation and RNA evidence',
        author='Isovar evidence audit',
    )
    intro = main_story(sample, data)
    first_break = next(i for i, item in enumerate(intro) if isinstance(item, PageBreak)) + 1
    intro = intro[:first_break] + landscape_story(sample, data) + intro[first_break:]
    story = (intro + coding_story(sample, data)
             + [PageBreak()] + evidence_story(sample, data) + read_support_story(data)
             + exonic_read_support_story(data) + sv_exon_story(data)
             + sv_junction_story(data)
             + altered_aa_story(data) + exemplar_story(sample, data)
             + ddx17_story(sample)
             + sv_exemplars(sample, data) + other_categories(sample, data)
             + provenance_story(sample, data))
    doc.build(story, onFirstPage=lambda c, d: page_footer(c, d, sample),
              onLaterPages=lambda c, d: page_footer(c, d, sample))
    focused = REPORT / 'focused_proteins' / (
        'HCC1395_EIF3K_CYP39A1_fusion.pdf' if sample == 'HCC1395' else 'COLO829_CDKN2A_frameshift.pdf')
    combined = path.with_suffix('.combined.pdf')
    with PdfWriter() as writer:
        writer.append(path, outline_item='Variant landscape and category exemplars')
        writer.append(focused, outline_item='Focused protein alignment, coverage and sources')
        with combined.open('wb') as stream:
            writer.write(stream)
    combined.replace(path)
    print(path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--annotations', type=Path, default=ANNOTATIONS)
    parser.add_argument('--output-dir', type=Path, default=OUT)
    args = parser.parse_args()
    for sample in ('HCC1395', 'COLO829'):
        build_one(sample, args.output_dir, args.annotations)
