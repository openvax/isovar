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
    EvidenceFunnel, ProteinComparison, altered_aa_figures, read_support_figures,
    support_by_class_figures,
)
from isovar.evidence_report import read_tsv as read_evidence_tsv, summarize, variant_type


rl_config.invariant = 1
ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / 'reports/hcc1395_colo829'
ANNOTATIONS = REPORT / 'reference_coding_orfdata_annotations.tsv.gz'
SELF_SCORES = REPORT / 'reference_self_9mer_scores.tsv.gz'
SV_EXONS = REPORT / 'sv_exon_annotations.tsv.gz'
SV_JUNCTION_AUDIT = REPORT / 'sv_junction_audit.tsv'
EXAMPLE_REFERENCES = REPORT / 'example_reference_contexts.tsv'
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
        row['call_id'] in exonic_ids and int(row['mapq20_junction_fragments']) > 0
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
    canvas.drawRightString(PAGE_W - 50, 32, f'2026-10-22  |  page {doc.page}')
    canvas.restoreState()


def section(title, body=None):
    result = [para(title, 'Section')]
    if body:
        result.append(para(body))
    return result


def main_story(sample, d):
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
    story += section('Structural variants: DNA to RNA breakpoint evidence')
    story.append(EvidenceFunnel(stages=(
        ('SV records in DNA source', c['sv_total']),
        ('Any annotated exon involved', c['sv_exonic']),
        ('RNA alignment near either breakend', c['sv_exonic_covered']),
        ('Split-read or cross-end mate candidate', c['sv_exonic_candidate']),
        ('Junction with MAPQ 20 reads (selected)', c['sv_exonic_junction']),
        ('SV-derived protein called', c['sv_exonic_product']),
    )))
    story.append(para(
        'Each funnel uses its own DNA source-record denominator. Breakend-window '
        'coverage is locus coverage, not evidence for the SV. A split-read or '
        'cross-end mate candidate links nominated windows but does not resolve '
        'junction bases, orientation, splice path or a translated product. '
        'The junction stage is a targeted original-read replay of selected events; '
        'it is a lower bound for the full source set.', 'FootSmall'))
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
            ('Calls with local protein contexts', f"{c['protein_calls']:,}",
             f"{c['protein_contexts']:,} ranked contexts; full-length isoforms unresolved"),
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
        'alignment records. The MAPQ 20 column requires every linked segment '
        'to pass that threshold. Direct junction fragments count reads crossing '
        'the same join, even when their wider transcript sequences differ. '
        'Full-path fragments require one read to match the entire retained '
        'candidate path and its junction. A zero in this table is a result '
        'for the listed alignment subset, not absence of expression.'))
    rows = []
    for row in d['sv_audit']:
        label = row['call_id'].replace('COLO829_SOMATIC_', '')
        parts = label.split('_')
        if len(parts) == 5 and parts[0] == parts[2] and parts[1].isdigit() and parts[3].isdigit():
            label = f'{parts[0]}:{int(parts[1]):,}-{int(parts[3]):,}'
        rows.append((label, row['class'], row['mapq1_junction_fragments'],
                     row['mapq20_junction_fragments'],
                     row['mapq20_max_full_path_fragments']))
    story.append(data_table(
        ['Nominated DNA record', 'Class', 'Junction MAPQ 1',
         'Junction MAPQ 20', 'Full path MAPQ 20'],
        rows, [203, 49, 83, 88, 89], numeric_columns=(2, 3, 4)))
    if d['sv_audit'][0]['sample'] == 'HCC1395':
        story.append(para(
            'The chr9 duplication has seven MAPQ 20 junction fragments but '
            'only one exact full-path witness for each of two retained paths. '
            'Path search reached its ten-path cap, so additional isoforms may '
            'exist. The chr13 deletion has three MAPQ 20 junction fragments '
            'and one exact whole-path witness; the chr21 duplication has two '
            'junction fragments and one per retained path. The three chr1 '
            'screen rows reuse one complex read and were not oriented in this '
            'targeted replay.'))
    else:
        story.append(para(
            'The chr1 SV3 deletion has six high-MAPQ opposite-mate candidate '
            'pairs and three apparent split junction fragments at MAPQ 1. '
            'Their linking segments have MAPQ 1-6; no junction survives MAPQ '
            '20. SV42 and SV49 each have one mate candidate and no observed '
            'junction at either threshold. The separate older benchmark '
            'records reuse those mate names and must not be added as new '
            'fragments. SV3 does not meet the annotated-exon rule used by '
            'the preceding funnel.'))
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


def example_story(d):
    c, top = d['counts'], d['summary']['top']
    calls = d['calls']
    story = section('Mutation classes and resolved outcomes')
    story.append(data_table(
        ['Category', 'DNA source', 'RNA ALT', 'Local context / product'],
        [
            ('DNA-predicted frameshift', f"{c['frameshift']:,}",
             f"{c['frameshift_alt']:,}", 'See per-call local contexts'),
            ('Reference ncORF coding-block overlap', f"{c['ncorf']:,}",
             f"{c['ncorf_alt']:,}", 'Mutated ncORF not reconstructed'),
            ('SV source record', f"{c['sv_total']:,}",
             'Breakend screen only', f"{c['sv_product']:,} resolved products"),
        ], [179, 80, 88, 165]))
    story.append(para(
        'Frameshifts are Varcode top-effect predictions from DNA, not a count of '
        'RNA-confirmed full-length frameshift proteins. SV split/mate candidates '
        'do not count as productive until an RNA junction and protein path resolve.',
        'FootSmall'))

    translated = [call for call in calls if call['call_id'] in top]
    translated.sort(key=lambda row: (-int(row['alt_reads']), row['call_id']))

    def examples(title, selected, *, score=False):
        result = section(title)
        if not selected:
            result.append(para('No qualifying local product was observed.', 'FootSmall'))
            return result
        rows = []
        for call in selected[:3]:
            ident = call['call_id']
            protein = top[ident]
            sequence = protein['protein_sequence']
            preview = sequence[:30] + ('...' if len(sequence) > 30 else '')
            if score:
                record = d['self_scores'][ident]
                measure = f"{record['novel_9mers']}/{record['total_9mers']} novel"
            else:
                measure = f"{int(call['alt_reads']):,} ALT reads"
            rows.append((ident, variant_type(call), measure, preview))
        result.append(data_table(['Source call', 'Type', 'Support / self screen',
                                  'Top local sequence'], rows,
                                 [165, 45, 115, 187]))
        return result

    story += examples('Highest original-RNA ALT support with local contexts', translated)
    if d['self_scores']:
        scored = [call for call in translated if call['call_id'] in d['self_scores']
                  and int(d['self_scores'][call['call_id']]['total_9mers'])]
        scored.sort(key=lambda row: (
            -int(d['self_scores'][row['call_id']]['novel_9mers']) /
            int(d['self_scores'][row['call_id']]['total_9mers']),
            -int(d['self_scores'][row['call_id']]['novel_9mers']),
            -int(row['alt_reads']), row['call_id']))
        story += examples('Least reference-self-like local contexts', scored, score=True)
        story.append(para(
            'The score is the fraction of distinct mutation-overlapping 9-mers absent '
            'from the pinned Ensembl 87 reference protein FASTA. It is an exact-match '
            'screen of these local windows, not patient-specific self, HLA binding, '
            'presentation or immunogenicity.', 'FootSmall'))
    else:
        story += section('Least reference-self-like local contexts')
        story.append(para('Reference-proteome comparison was not evaluated.', 'FootSmall'))
    non_snv = [call for call in translated if call['class'] != 'SNV']
    story += examples('Highest supported non-SNV local products', non_snv)
    story.append(para(
        'The tables rank source records, not independent molecules. A displayed '
        'local context may represent only part of a translated protein and can '
        'be reference-identical for a synonymous allele.', 'FootSmall'))
    story.append(PageBreak())
    return story


def comparison_story(d):
    """Show pinned reference and read-supported mutant local protein windows."""
    calls = {row['call_id']: row for row in d['calls']}
    top = d['summary']['top']
    story = section('Reference and mutant protein windows')
    story.append(para(
        'The <font color="#ba5a29">orange reference</font> and '
        '<font color="#0f706d">teal mutant</font> residues differ at the '
        'marked boundary; dark residues are shared flanks. Each reference '
        'window comes from the named Ensembl 87 transcript. The mutant is '
        'Isovar\'s top-ranked local context supported by tumor RNA reads. '
        'For frameshifts, downstream columns are anchored at the change, '
        'not a homology alignment or a complete protein. A dash means the '
        'reference protein ends before the mutant window.'))
    for row in d['example_references']:
        call = calls[row['call_id']]
        product = top[row['call_id']]
        mutation_start = int(product['mutation_start_idx'])
        mutation_end = int(product['mutation_end_idx'])
        point_change = (variant_type(call) == 'SNV' and
                        mutation_end == mutation_start + 1)
        story.append(para(
            f"<b>{row['call_id']}</b> | {variant_type(call)} | "
            f"{int(call['alt_reads']):,} ALT reads | "
            f"{row['transcript_id']} from amino acid {row['protein_start_1based']}",
            'FootSmall'))
        story.append(ProteinComparison(
            row['reference_context'], product['protein_sequence'],
            mutation_start, mutation_end, point_change=point_change))
        story.append(Spacer(1, 7))
    story.append(PageBreak())
    return story


def provenance_story(sample, d):
    story = section('Sources, processing and limits')
    if sample == 'HCC1395':
        sources = [
            ('Matched-normal small variants', 'SEQC2 v1.2.1 SNV + indel, GRCh38.d1.vd1; HCC1395BL normal'),
            ('SV benchmark', 'Talsania et al. 2022 Table S4, HCC1395 vs HCC1395BL'),
            ('Original RNA', 'PRJNA934933 / SRR23456709; 1,049,452 ONT direct RNA reads'),
            ('Alignment', 'minimap2 2.26 -ax splice -uf -k14; NCBI GRCh38 no-alt'),
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
            'reads may be lost. A weak MAPQ 1 join at SV3 fails the MAPQ 20 gate; '
            'full-length isoform linkage and SV-derived proteins remain unresolved.'
        )
    story.append(data_table(['Input', 'Pinned scope'], sources, [160, 352]))
    story += section('Reference model and interpretation')
    story.append(para(
        'Small-variant coding effects use Varcode with Ensembl release 87 on GRCh38. '
        'The human orfdata 0.34.1 catalog uses Ensembl 101 GRCh38 reference ORFs '
        '(7,264 models, dataset 2026-10-11.1). An orfdata coding-block overlap '
        'is a reference annotation candidate; this analysis did not translate '
        'mutated ncORFs from the RNA reads. Genomic transcript overlap is not '
        'RNA-compatible isoform support. Original-read allele counts and local '
        'protein contexts are separate observations.'))
    story += section('SHA-256 input pins')
    for name, digest in hashes:
        story.append(para(f'<b>{name}</b>', 'FootSmall'))
        story.append(para(digest, 'TinyMono'))
    story.append(para('<b>GRCh38 no-alt FASTA (uncompressed)</b>', 'FootSmall'))
    story.append(para('9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7',
                      'TinyMono'))
    story.append(para('<b>orfdata bundled GRCh38 catalog</b>', 'FootSmall'))
    story.append(para('06bc0bf3f81c25a7356f71a38581ff5704c05dcb75752454bc1ca75075ab8c09',
                      'TinyMono'))
    story.append(para('<b>Ensembl 87 GRCh38 protein FASTA gzip (self screen)</b>', 'FootSmall'))
    story.append(para(SELF_FASTA_SHA256, 'TinyMono'))
    story += section('Other antigen sources in this run')
    story.append(data_table(
        ['Category', 'Assessment'],
        [
            ('Viral RNA and viral coding paths', 'Not assessed: no viral-aligned RNA panel'),
            ('Supplied fusion transcripts', 'Not assessed with direct junction reconstruction'),
            ('Targeted antigen-gene isoforms', 'Not assessed with original-read isoform coverage'),
            ('Mutated ncORF proteins', 'Not reconstructed; reference coding-block overlaps only'),
        ], [210, 302]))
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
    data = collect(sample, annotations_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f'{sample}_somatic_RNA_evidence.pdf'
    doc = SimpleDocTemplate(
        str(path), pagesize=letter, leftMargin=50, rightMargin=50,
        topMargin=52, bottomMargin=51, title=f'{sample} somatic mutation and RNA evidence',
        author='Isovar evidence audit',
    )
    story = (main_story(sample, data) + coding_story(sample, data)
             + [PageBreak()] + evidence_story(sample, data) + read_support_story(data)
             + exonic_read_support_story(data) + sv_exon_story(data)
             + sv_junction_story(data)
             + altered_aa_story(data) + example_story(data)
             + comparison_story(data)
             + provenance_story(sample, data))
    doc.build(story, onFirstPage=lambda c, d: page_footer(c, d, sample),
              onLaterPages=lambda c, d: page_footer(c, d, sample))
    print(path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--annotations', type=Path, default=ANNOTATIONS)
    parser.add_argument('--output-dir', type=Path, default=OUT)
    args = parser.parse_args()
    for sample in ('HCC1395', 'COLO829'):
        build_one(sample, args.output_dir, args.annotations)
