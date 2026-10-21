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


rl_config.invariant = 1
ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / 'reports/hcc1395_colo829'
ANNOTATIONS = REPORT / 'reference_coding_orfdata_annotations.tsv.gz'
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
    annotations = {row['call_id']: row for row in read_gzip(annotations_path)
                   if row['sample'] == sample}
    assert len(calls) == len(annotations)
    assert {c['call_id'] for c in calls} == set(annotations)
    assert all(not row['annotation_error'] for row in annotations.values())
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
    return {
        'calls': calls, 'proteins': proteins, 'sv': sv, 'counts': counts,
        'rna': rna, 'classes': classes, 'filters': filters,
        'sv_classes': sv_classes, 'sv_status': sv_status,
    }


def page_footer(canvas, doc, sample):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor('#c6d5dc'))
    canvas.line(50, PAGE_H - 39, PAGE_W - 50, PAGE_H - 39)
    canvas.setFont('Helvetica-Bold', 7.5)
    canvas.setFillColor(TEAL)
    canvas.drawString(50, PAGE_H - 31, f'ISOVAR 1.59.7  /  {sample} EVIDENCE AUDIT')
    canvas.setFont('Helvetica', 7.5)
    canvas.setFillColor(GRAY)
    canvas.drawRightString(PAGE_W - 50, 32, f'2026-10-21  |  page {doc.page}')
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
    story += section('At a glance')
    story.append(data_table(
        ['Evidence layer', 'Count', 'Interpretation'],
        [
            ('Ingested small variants', f'{small:,}',
             f"{d['classes']['SNV']:,} SNVs; {d['classes']['indel']:,} indels"),
            ('Ingested SV source records', f'{sv_total:,}',
             'DNA breakpoint records/events, counted in their source sets'),
            ('Ensembl 87 coding effect', f"{c['coding']:,}",
             'Any annotated protein-coding transcript, reference model'),
            ('orfdata ncORF coding-block overlap', f"{c['ncorf']:,}",
             'Candidate reference ORF overlap; distinct from mutant translation'),
            ('Original RNA ALT-observed calls', f"{d['rna']['alternate_observed']:,}",
             f"{c['alt_reads']:,} ALT reads; {c['alt_fragments']:,} ALT fragments"),
            ('Calls with local mutant protein contexts', f"{c['protein_calls']:,}",
             f"{c['protein_contexts']:,} ranked contexts; full-length isoforms unresolved"),
        ], [179, 68, 265], numeric_columns=(1,)))
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
        f"read-compatible <b>local</b> mutant protein context for {c['protein_calls']:,} "
        f"calls ({c['protein_contexts']:,} ranked rows). The remaining "
        f"{rna['alternate_observed'] - c['protein_calls']:,} ALT-observed calls have "
        'no local translation under the stated method. No-overlap and uncallable statuses '
        'cannot establish absence of RNA expression.'))
    if sample == 'HCC1395':
        story += section('Representative mutant residues with original-read support')
        story.append(data_table(
            ['SEQC2 call', 'Top local amino-acid context', 'Mutant interval', 'Supporting reads'],
            [
                ('OSTC F9L / COSM32736', 'METLYRVPLLVLECPNLKLKKPPWL', '[8,9)', '5'),
                ('PRDX5 F157L / COSM32759', 'VRLLADPTGALGKETDLLLD', '[10,11)', '6'),
            ], [122, 208, 82, 100], numeric_columns=(3,)))
        story.append(para(
            'The SEQC2 tumor/normal DNA VAFs are 0.490/0.001 (OSTC) and 0.442/0.002 '
            '(PRDX5). An independent proteomic study observed the mutant residues, not '
            'these entire local sequence contexts. A separate first-100,000-read ONT '
            'cDNA fixture has a PRDX5 local context with six supporting reads; its '
            'bounded counts are not included in the whole-run denominators.', 'FootSmall'))
    else:
        story += section('High-depth indel and polished sequence comparison')
        story.append(para(
            'At CDKN2A chr9:21971154:CCG&gt;C, all 6,984 ALT allele reads were evaluated. '
            'The top local frameshift context is <b>MMMGSARVAELLLLHGGAQLRRPR</b> '
            'with 5,338 supporting reads. Its 1,416 ranked local contexts reflect '
            'read-spanning alternatives; none is asserted to be a full-length isoform.'))
        story.append(para(
            'TAGET polished representatives are an independent culture and assay: 132 '
            'SMaHT SNVs have at least one tumor-positive representative carrying ALT. '
            'The 0/1 representative-sequence observations and full-length cluster '
            'counts are separate from original allele-spanning read support. The '
            'polished count table has 51,387 tumor-positive transcript IDs, but those '
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
        story.append(para(
            'The 54 Paulin Paz v2 records and 134 Valle-Inclan lifted records '
            'remain separate benchmarks. Eleven VCF records have opposite-mate '
            'context, with no observed same-read split candidate in the retained '
            'BAM. No junction or SV-derived protein is called.', 'FootSmall'))
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
            'Full-length isoform linkage and SV protein reconstruction were not assessed.'
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
            'reads may be lost. Full-length isoform linkage, SV junction paths and '
            'SV-derived proteins were not resolved.'
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
    story = main_story(sample, data) + evidence_story(sample, data) + provenance_story(sample, data)
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
