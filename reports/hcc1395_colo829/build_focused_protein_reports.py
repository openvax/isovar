"""Render pinned CDKN2A and EIF3K-CYP39A1 sequence/support reports offline."""

import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import shutil
import textwrap

from reportlab.lib.colors import HexColor
from reportlab.pdfgen.canvas import Canvas

from isovar.genetic_code import translate_cdna
from reports.hcc1395_colo829.build_hcc1395_eif3k_fusion import (
    ACCEPTOR_ID, DONOR_ID, FIXTURE, load_references, report,
)

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'reports/hcc1395_colo829'
FOCUSED = BASE / 'focused_proteins'
OUT = ROOT / 'output/pdf'
REPO = 'https://github.com/openvax/isovar/blob/master/'
NAVY, TEAL, ORANGE, GRAY = map(HexColor, ['#18334b', '#087d8a', '#b84a13', '#526373'])
GERMLINE, OTHER_SOMATIC = map(HexColor, ['#275bcc', '#8d3fa8'])


def rows(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


class Page:
    """Small fixed-layout report with explicit page boundaries."""

    def __init__(self, path, title):
        self.canvas = Canvas(str(path), pagesize=(612, 792), invariant=1)
        self.canvas.setTitle(title)
        self.title, self.number, self.y = title, 0, 0

    def new(self, subtitle):
        if self.number:
            if self.y < 48:
                raise ValueError('Report page overflow')
            self.footer()
            self.canvas.showPage()
        self.number += 1
        self.y = 746
        self.text(self.title, 19, bold=True)
        self.text(subtitle, 10, color=GRAY)
        self.y -= 8

    def footer(self):
        c = self.canvas
        c.setFillColor(GRAY)
        c.setFont('Helvetica', 8)
        c.drawString(44, 28, 'OpenVax Isovar | source-specific RNA evidence | 2026-10-26')
        c.drawRightString(568, 28, str(self.number))

    def text(self, value, size=10, color=NAVY, bold=False):
        self.canvas.setFillColor(color)
        self.canvas.setFont('Helvetica-Bold' if bold else 'Helvetica', size)
        self.canvas.drawString(44, self.y, value)
        self.y -= size + 7

    def paragraph(self, value, size=10, width=96):
        for line in textwrap.wrap(value, width):
            self.text(line, size)
        self.y -= 4

    def section(self, value):
        self.y -= 9
        self.text(value, 12, color=TEAL, bold=True)

    def sequence(self, label, sequence, highlights=(), colors=None, size=11, step=None):
        c = self.canvas
        c.setFillColor(GRAY)
        c.setFont('Helvetica', 9)
        c.drawString(44, self.y, label)
        c.setFont('Courier', size)
        for i, aa in enumerate(sequence):
            color = colors[i] if colors is not None else ORANGE if any(a <= i < b for a, b in highlights) else NAVY
            c.setFillColor(color)
            c.drawString(165 + i * size * .6, self.y, aa)
        self.y -= step if step is not None else size + 9

    def source(self, title, detail, url=None):
        self.text(title, 10, bold=True)
        if url:
            self.canvas.linkURL(url, (44, self.y + 4, 568, self.y + 19), relative=0)
        self.paragraph(detail, 9, 106)

    def save(self):
        if self.y < 48:
            raise ValueError('Report page overflow')
        self.footer()
        self.canvas.save()


def variant_legend(p):
    p.text('Variant provenance', 10, bold=True)
    c = p.canvas
    for x, label, color in [(44, 'Focal change', ORANGE), (190, 'Germline', GERMLINE),
                             (325, 'Phased somatic', OTHER_SOMATIC)]:
        c.setFillColor(color)
        c.rect(x, p.y - 1, 8, 8, fill=1, stroke=0)
        c.setFillColor(NAVY)
        c.setFont('Helvetica', 9)
        c.drawString(x + 14, p.y, label)
    p.y -= 22
    p.paragraph('No additional non-reference residues occur in the displayed RNA sequence relative '
                'to its indicated reference contributor(s). No germline or other somatic allele is '
                'assigned without genotype and phasing evidence. Missing sequence stays unobserved.', 9)


def cdkn2a():
    call_id = 'chr9:21971154:CCG>C'
    calls_path = BASE / 'colo829_atlas_calls.tsv.gz'
    proteins_path = BASE / 'colo829_atlas_proteins.tsv.gz'
    reference_path = BASE / 'example_reference_contexts.tsv'
    full_reference_path = FOCUSED / 'cdkn2a_reference.json'
    full_reference = json.loads(full_reference_path.read_text())
    call, = [r for r in rows(calls_path) if r['call_id'] == call_id]
    proteins = sorted((r for r in rows(proteins_path) if r['call_id'] == call_id), key=lambda r: int(r['rank']))
    reference, = [r for r in rows(reference_path) if r['call_id'] == call_id]
    top = proteins[0]
    start, end = int(top['mutation_start_idx']), int(top['mutation_end_idx'])
    wt, mutant = reference['reference_context'], top['protein_sequence']
    assert wt[:start] == mutant[:start] and len(wt) == len(mutant) == 24
    assert int(reference['protein_start_1based']) + start == 68
    offset = int(reference['protein_start_1based']) - 1
    full_wt = full_reference['sequence']
    assert full_wt[offset:offset + len(wt)] == wt
    assert full_reference['protein_sha256'] == reference['reference_protein_sha256']
    value = dict(sample='COLO829', call=call, reference=reference, top_protein=top,
                 retained_contexts=len(proteins), inputs={str(p.relative_to(ROOT)): digest(p)
                 for p in (calls_path, proteins_path, reference_path, full_reference_path)},
                 full_reference=full_reference,
                 residue_provenance=dict(focal_variant=[start, end], germline=[], phased_other_somatic=[],
                                         unclassified_nonreference=[], unobserved_outside_local_window=True))
    pdf = OUT / 'COLO829_CDKN2A_frameshift.pdf'
    p = Page(pdf, 'COLO829 | CDKN2A frameshift')
    p.new('High-support local mutant protein window - not a complete mutant protein')
    p.paragraph('GRCh38 chr9:21971154 CCG>C. A two-base deletion changes the coding frame in the '
                'Ensembl 87 CDKN2A reference shown below (ENST00000304494 / ENSP00000307101).')
    p.text('IGV (GRCh38): chr9:21971154-21971156', 10, bold=True)
    p.section('RNA support: allele and complete window are different measurements')
    p.text(f"{int(call['alt_fragments']):,} ALT fragments / {int(call['alt_reads']):,} ALT allele-read observations", 15, bold=True)
    p.text(f"{int(top['supporting_fragments']):,} fragments support this entire 24-aa window", 15, color=TEAL, bold=True)
    p.paragraph(f"The window has {int(top['supporting_reads']):,} supporting allele-read observations. "
                f"The locus also has {call['ref_fragments']} REF fragments and {call['other_fragments']} other-allele fragments. "
                'These are template counts in one Illumina library, not UMI-confirmed molecules.')
    p.section('Reference versus observed mutant window')
    p.text('Reference residues 52-75; the altered tail starts at reference residue 68.', 9)
    p.sequence('Reference 52-75', wt)
    p.sequence('RNA mutant', mutant, [(start, end)])
    p.sequence('Change', ' ' * start + '^' * (end - start))
    p.paragraph('The highlighted mutant tail begins after a 16-aa shared prefix; altered tail: '
                '8 aa (GAQLRRPR). Columns after the change are boundary-aligned, not a homology '
                'alignment. No stop or full-length mutant protein is resolved by this local window.', 9)
    p.section('Retained alternatives - never add their support')
    for row in proteins[:3]:
        p.text(f"Rank {row['rank']}: {len(row['protein_sequence'])} aa, "
               f"{int(row['supporting_fragments']):,} complete-window fragments", 10)
    p.paragraph(f'{len(proteins):,} overlapping local contexts remain in the atlas. Eight transcript IDs '
                'are compatible with the displayed window; this does not select a unique isoform. '
                'All ALT reads were evaluated by the recorded high-depth read-spanning method.', 9)
    p.section('Interpretation')
    p.paragraph('This is a strongly RNA-supported frameshift window. RNA support does not establish '
                'translated protein abundance, tumor specificity in another culture, or antigen presentation.', 9)
    p.new('Reference protein, observed mutant window, and assembly extent')
    p.paragraph('The full reference protein is shown for orientation. Only the 24-aa RNA window is '
                'resolved in this report; dots in the mutant row mean unobserved sequence, not gaps '
                'filled from reference. Coordinates are anchored to the reference before the frameshift.', 9)
    full_mutant = '.' * offset + mutant + '.' * (len(full_wt) - offset - len(mutant))
    for left in range(0, len(full_wt), 60):
        right = min(left + 60, len(full_wt))
        p.text(f'Reference coordinates {left + 1}-{right}', 8, color=GRAY)
        colors = [ORANGE if offset + start <= i < offset + end else
                  TEAL if offset <= i < offset + len(mutant) else GRAY for i in range(left, right)]
        p.sequence('Reference', full_wt[left:right], size=9, step=15)
        p.sequence('RNA mutant', full_mutant[left:right], colors=colors, size=9, step=22)
    p.section('Assembly / reference annotation track')
    c, y, width = p.canvas, p.y - 18, 510
    c.setFillColor(HexColor('#dce3e8'))
    c.rect(44, y, width, 13, fill=1, stroke=0)
    c.setFillColor(TEAL)
    c.rect(44 + width * offset / len(full_wt), y, width * len(mutant) / len(full_wt), 13, fill=1, stroke=0)
    c.setFillColor(ORANGE)
    c.rect(44 + width * (offset + start) / len(full_wt), y,
           width * (end - start) / len(full_wt), 13, fill=1, stroke=0)
    p.y = y - 22
    for x, label, color in [(44, 'Reference only', HexColor('#dce3e8')),
                             (212, 'RNA observed', TEAL), (380, 'Altered tail', ORANGE)]:
        c.setFillColor(color)
        c.rect(x, p.y - 1, 8, 8, fill=1, stroke=0)
        c.setFillColor(NAVY)
        c.setFont('Helvetica', 9)
        c.drawString(x + 14, p.y, label)
    p.y -= 23
    p.paragraph('RNA coverage is restricted to reference-anchored positions 52-75 in this selected '
                'window, supported by 5,248 reported fragments. The boundaries at 52 and 76 '
                'mark the observed window, not validated full-ORF ends. No reference sequence is '
                'stitched into the mutant sequence. Upstream initiation and downstream stop remain unresolved.', 9)
    variant_legend(p)
    p.new('All inputs used in this focused report')
    p.source('Original RNA: SRR8615617 / PRJNA523380 / SAMN10988208',
             'COLO829 CCLE Illumina paired cDNA RNA-seq: 115,269,780 pairs, 101 bases per end. '
             'The whole run was aligned to GRCh38 before selected-site extraction. The atlas BAM '
             'is a derivative of this same library, not independent evidence.',
             'https://www.ncbi.nlm.nih.gov/sra/SRR8615617')
    p.source('DNA nomination: SMaHT COLO829BLT50 reference set',
             'Truthset_Indel_COLO829.vcf.gz; PASS, COLO829T-origin call relative to COLO829BL. '
             'DNA and RNA cultures come from different studies.',
             'https://github.com/parklab/SMaHT_SNV_COLO829BLT50_HAPMAP/tree/main/Resource/COLO829BLT50_reference_set')
    p.source('Genome and protein reference',
             'NCBI GRCh38 no-alt analysis FASTA; Ensembl 87 transcript ENST00000304494. '
             'Pinned reference window starts at residue 52. Protein SHA-256: ' + reference['reference_protein_sha256'],
             REPO + 'reports/hcc1395_colo829/example_reference_contexts.tsv')
    p.source('Processing and exact original witnesses',
             'minimap2 2.26 -ax sr --frag=yes (genomic, not splice-aware); Isovar 1.59.6 atlas, '
             'default read handling and high-depth direct contexts. The call, protein and ALT-witness '
             'tables retain the source-run IDs. Re-rendering does not rerun the atlas.',
             REPO + 'reports/hcc1395_colo829/COLO829_ATLAS.md')
    p.source('Complete study/source index and pinned inputs',
             'SEQUENCING_SOURCES.md lists all sources across the three-sample reports. The adjacent '
             'COLO829_CDKN2A_frameshift.json pins the call, reference, sequence and table SHA-256 digests.',
             REPO + 'reports/hcc1395_colo829/SEQUENCING_SOURCES.md')
    p.section('Compatible transcript IDs')
    p.paragraph(top['compatible_transcript_ids'].replace(';', ', '), 9)
    p.save()
    return pdf, value, wt, mutant


def fusion():
    paths = {name: FIXTURE / name for name in ('original_unmapped.bam', 'models.json.gz',
             'ont_candidate_sequences.fasta.gz', 'coding_tiling.bam', 'ont_candidate_alignments.bam',
             'ont_exon_models.json', 'dna_tumor_loci.bam', 'dna_normal_loci.bam')}
    v = report(*paths.values())
    references = {r.transcript_id: r for r in load_references(paths['models.json.gz'])}
    donor, acceptor = (references[i] for i in (DONOR_ID, ACCEPTOR_ID))
    dp, _ = translate_cdna(donor.sequence[donor.cds_start:donor.cds_end])
    ap, _ = translate_cdna(acceptor.sequence[acceptor.cds_start:acceptor.cds_end])
    protein, = v['proteins']
    junction = v['fusion'].junction_start
    codon = junction // 3
    acceptor_offset = 517 - acceptor.cds_start
    acceptor_codon = acceptor_offset // 3
    assert junction == 499 and codon == 166 and acceptor_codon == 104
    assert protein[:codon] == dp[:codon]
    assert protein[codon:] == ap[acceptor_codon:]
    assert len(protein) == 531 and len(v['fusion'].sequence) == 1596
    coverage = [len({r.fragment_id for r in v['reads']
                    if r.cdna_start <= i and r.cdna_start + len(r.sequence) >= i + 3})
                for i in range(0, len(v['fusion'].sequence), 3)]
    assert min(coverage) > 0
    value = dict(sample='HCC1395', protein=protein, donor_protein=dp, acceptor_protein=ap,
                 junction_cdna_offset=junction, junction_protein_position_1based=codon + 1,
                 acceptor_reference_position_1based=acceptor_codon + 1,
                 exact_fragment_coverage_per_codon_including_stop=coverage,
                 ont=v['ont'], dna=v['dna'], inputs={str(p.relative_to(ROOT)): digest(p) for p in paths.values()})
    value['residue_provenance'] = dict(focal_junction=[codon, codon + 1], germline=[],
                                     phased_other_somatic=[], unclassified_nonreference=[], reference_stitched_bases=0)
    support = v['result']['paths'][0]['junctions'][0]['direct_support']['fragments']
    value['junction_fragments'] = support
    pdf = OUT / 'HCC1395_EIF3K_CYP39A1_fusion.pdf'
    p = Page(pdf, 'HCC1395 | EIF3K-CYP39A1 fusion')
    p.new('531-aa complete coding candidate - RNA reconstruction, not measured translation')
    p.text(f'{support} exact Illumina junction fragments', 16, bold=True)
    p.text('9 ONT coding-span reads; 4 also have the exact junction', 13, color=TEAL, bold=True)
    p.paragraph('The reference-guided 1,596-nt coding RNA path is tiled by exact original reads, '
                'including 19 mapped Illumina reads with overlaps of at least 20 nt. The Illumina '
                'input is fusion-enriched. Counts do not measure whole-library prevalence.')
    p.text('IGV (GRCh37): chr19:39123318-39123318 chr6:46607405-46607405', 8.5, bold=True)
    p.section('Protein architecture and reference alignment')
    p.paragraph('Fusion residues 1-166 match EIF3K residues 1-166. Residue 167 is the split junction '
                'codon: one EIF3K base plus two CYP39A1 bases. Its alanine matches CYP39A1 residue '
                '105. Fusion residues 168-531 match CYP39A1 residues 106-469, followed by a stop.', 10)
    p.text('Junction zoom: aligned by the pinned exon/CDS map', 10, bold=True)
    lo, hi = codon - 16, codon + 20
    p.sequence('EIF3K 151-186', dp[lo:hi], size=10)
    p.sequence('Fusion 151-186', protein[lo:hi], [(16, 17)], size=10)
    p.sequence('CYP39A1 89-124', ap[acceptor_codon - 16:acceptor_codon + 20], size=10)
    p.paragraph('The mutant highlight marks junction codon 167. The novel feature is the joined protein context; '
                'the CYP39A1-derived suffix is reference sequence, not a novel frameshift tail. '
                'Other reported fusion isoforms remain possible.', 9)
    p.section('Exact original RNA support across the entire coding path')
    c, y, height = p.canvas, p.y - 105, 90
    c.setStrokeColor(TEAL)
    for i, count in enumerate(coverage):
        x = 60 + i * 490 / (len(coverage) - 1)
        c.line(x, y, x, y + height * count / max(coverage))
    c.setStrokeColor(ORANGE)
    x = 60 + codon * 490 / (len(coverage) - 1)
    c.line(x, y, x, y + height + 7)
    c.setFillColor(GRAY)
    c.setFont('Helvetica', 8)
    c.drawString(60, y - 14, '1')
    c.drawString(x - 15, y - 14, '167: join')
    c.drawRightString(550, y - 14, '531 + stop')
    c.drawString(60, y + height + 9, f'Fragments per codon: min {min(coverage)}, max {max(coverage)} (selected exact-read subset)')
    c.setFillColor(TEAL)
    c.rect(60, y - 30, 490, 6, fill=1, stroke=0)
    p.y = y - 48
    p.text('Assembly track: all 1,596 coding bases RNA-witnessed; 0 reference-filled bases.', 8.5)
    p.paragraph('Each bar counts distinct original fragment IDs with an exact read covering the '
                'whole codon. It is not raw BAM depth. The RNA library and ONT library remain separate. '
                'Six of 21 ONT candidate reads contain the exact 20+20-base junction; nine span '
                'coding start to stop in local collinear alignments. Four satisfy both checks.', 9)
    p.new('Complete protein aligned to its two reference contributors')
    p.paragraph('D = retained EIF3K reference; F = fusion candidate; A = retained CYP39A1 reference. '
                'Dashes mean the contributor is absent at that fusion position, not a deletion from '
                'the other gene. The fusion row distinguishes the CYP39A1 suffix and split junction codon. Coordinates '
                'at left refer to the fusion protein. The terminal * is the observed stop.', 9)
    donor_row = dp[:codon] + '-' * (len(protein) - codon)
    acceptor_row = '-' * codon + ap[acceptor_codon:]
    alignment = []
    for left in range(0, len(protein), 60):
        right = min(left + 60, len(protein))
        p.text(f'Fusion amino acids {left + 1}-{right}', 8, color=GRAY)
        colors = [ORANGE if i == codon else TEAL if i > codon else NAVY for i in range(left, right)]
        for label, seq in [('D', donor_row), ('F', protein), ('A', acceptor_row)]:
            p.sequence(label, seq[left:right], colors=colors if label == 'F' else None, size=9, step=12)
        alignment.extend([f'D {donor_row[left:right]}', f'F {protein[left:right]}', f'A {acceptor_row[left:right]}', ''])
        p.y -= 1
    p.text('Terminal stop: *  |  Annotated EIF3K start and CYP39A1 stop retained.', 9)
    p.new('Full reference proteins and provenance of amino-acid changes')
    for name, seq, interval in [('EIF3K / ' + DONOR_ID, dp, (0, codon)),
                                ('CYP39A1 / ' + ACCEPTOR_ID, ap, (acceptor_codon, len(ap)))]:
        p.section(name)
        for left in range(0, len(seq), 60):
            right = min(left + 60, len(seq))
            colors = [TEAL if interval[0] <= i < interval[1] else GRAY for i in range(left, right)]
            p.sequence(f'{left + 1}-{right}', seq[left:right], colors=colors, size=9, step=17)
    p.paragraph('Highlighted reference residues contribute to the fusion. Dimmed residues are not retained. '
                'The fusion codon at residue 167 contains bases from both partners and encodes '
                'CYP39A1 alanine 105; it is not counted as a complete EIF3K-derived codon.', 9)
    variant_legend(p)
    p.paragraph('There is no reference-filled region in the determined ORF. The gene junction at '
                'coding offset 499 is an RNA-supported splice/fusion boundary, not an assembly-to-reference '
                'stitching point. Support falls outside the enriched junction region, as shown by '
                'the codon track; no junction count is assigned to every residue.', 9)
    p.new('All sequencing and reference inputs used in this focused report')
    p.source('Illumina tumor RNA: public INTEGRATE test-data.tar.gz',
             'Original unmapped.bam (185 fusion-enriched records) and accepted_hits.bam '
             '(19 mapped records retained for coding tiling). Archive SHA-256: '
             '9a37130e768e424935008f1d0a67050f6e11a242a263ed2f20ba071d94b767cb. '
             'The archive is the source identifier; no SRA run accession is asserted for this subset.',
             'https://sourceforge.net/p/integrate-fusion/wiki/example/')
    p.source('Independent ONT direct RNA: SRR23456709 / PRJNA934933',
             'BioSample SAMN33292477; 1,049,452 reads in the source run. The report checks '
             '21 original candidates and 42 local alignment records. This is a local candidate '
             'screen, not a genome-wide uniqueness test.', 'https://www.ncbi.nlm.nih.gov/sra/SRR23456709')
    p.source('Tumor and matched-normal DNA: INTEGRATE example archive',
             'All original alignments in two GRCh37 flanking intronic windows were retained. '
             '18 tumor cross-gene fragments versus 0 normal. This supports a rearrangement, '
             'without resolving its exact DNA breakpoint.', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC4691743/')
    p.source('Exon/cDNA and coding references',
             'Pinned GRCh37-lite archive models: EIF3K ENST00000248342 and CYP39A1 ENST00000275016. '
             'Ensembl 87 CDS offsets transfer only after exact whole-transcript sequence identity. '
             'ONT placements use GRCh38 and separately pinned Ensembl 87 exon coordinates; '
             'genomic coordinates from the two assemblies are never mixed.',
             REPO + 'tests/data/hcc1395_eif3k_cyp39a1/README.md')
    p.source('Derived files, checksums and broader report source index',
             'The adjacent HCC1395_EIF3K_CYP39A1_fusion.json pins all eight BAM/FASTA/model inputs '
             'and contains the complete protein, reference contributors and per-codon support. '
             'SEQUENCING_SOURCES.md lists every library, DNA source and reference used in the '
             'Sid/HCC1395/COLO829 report set.', REPO + 'reports/hcc1395_colo829/SEQUENCING_SOURCES.md')
    p.paragraph('The first page highlights strong junction support. It does not assign 149 fragments '
                'to every amino acid or claim that 149 reads span the entire protein. Complete coding '
                'sequence support is the separately displayed overlap tiling, corroborated by the '
                'source-separated ONT checks.', 9)
    p.save()
    return pdf, value, '\n'.join(alignment)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FOCUSED.mkdir(parents=True, exist_ok=True)
    pdf, value, wt, mutant = cdkn2a()
    (FOCUSED / 'COLO829_CDKN2A_frameshift.md').write_text(
        '# COLO829 CDKN2A frameshift\n\n'
        '[Open the focused PDF](COLO829_CDKN2A_frameshift.pdf).\n\n'
        'The 24-aa local window has **5,248 complete-window fragments** (5,338 allele-read observations); '
        'the allele has 6,967 ALT fragments (6,984 allele-read observations).\n\n'
        'Ensembl 87 ENST00000304494 reference residues 52-75; frameshift begins at residue 68.\n\n'
        '**IGV (GRCh38): `chr9:21971154-21971156`.**\n\n'
        f'```text\nReference  {wt}\nRNA mutant {mutant}\n           {" " * 16}{"^" * 8}\n```\n\n'
        'Boundary-aligned local window; no full mutant protein or stop is resolved. '
        'All 1,416 competing contexts are retained in the atlas.\n\n'
        'The PDF also shows the complete reference protein and an assembly track: only residues '
        '52-75 are represented by this RNA window. No unobserved flank is filled from reference. '
        'Orange denotes the focal change; blue is reserved for genotype-supported germline changes '
        'and purple for other somatic changes phased onto the same RNA. Neither additional class '
        'is observed in this window.\n\n'
        'Sources: SRR8615617 (Illumina RNA), SMaHT COLO829BLT50 indel truth set, NCBI GRCh38 no-alt, '
        'Ensembl 87. See the [complete source index](../SEQUENCING_SOURCES.md) and adjacent JSON input hashes.\n')
    (FOCUSED / (pdf.stem + '.json')).write_text(json.dumps(value, indent=2) + '\n')
    (FOCUSED / 'candidate_proteins.fasta').write_text(
        '>COLO829_CDKN2A local_window=24aa focal_interval_0based=16:24 complete_window_fragments=5248 full_protein_resolved=false\n'
        + mutant + '\n')
    shutil.copyfile(pdf, FOCUSED / pdf.name)
    pdf, value, alignment = fusion()
    (FOCUSED / 'HCC1395_EIF3K_CYP39A1_fusion.md').write_text(
        '# HCC1395 EIF3K-CYP39A1 fusion\n\n'
        '[Open the focused PDF](HCC1395_EIF3K_CYP39A1_fusion.pdf).\n\n'
        '**531-aa complete coding candidate**: 149 exact Illumina junction fragments; exact coding '
        'tiling includes 19 mapped reads. Separately, 9 ONT coding-span reads, 4 with an exact junction; '
        '6/21 ONT candidates have an exact junction overall.\n\n'
        '**IGV (GRCh37): `chr19:39123318-39123318 chr6:46607405-46607405`.**\n\n'
        'EIF3K residues 1-166 are retained. Fusion residue 167 is a split codon (1 EIF3K base + '
        '2 CYP39A1 bases) encoding the same alanine as CYP39A1 residue 105; the remaining sequence '
        'matches CYP39A1 residues 106-469. The joined context is novel, not the entire suffix.\n\n'
        'D = retained EIF3K; F = fusion; A = retained CYP39A1. Dashes mean a contributor is absent.\n\n'
        f'```text\n{alignment}\nTerminal stop: *\n```\n\n'
        'The PDF includes both complete reference proteins, a per-codon exact-fragment coverage '
        'track, and an assembly track. All 1,596 coding bases are RNA-witnessed; no bases are '
        'filled from reference. The junction is highlighted separately from retained parental '
        'sequence. No additional germline or phased somatic amino-acid change is observed.\n\n'
        'Sources: INTEGRATE public example RNA and tumor/normal DNA, independent SRR23456709 ONT RNA, '
        'archive GRCh37-lite transcripts with sequence-verified Ensembl 87 CDS, and separate GRCh38 '
        'ONT exon models. See the [complete source index](../SEQUENCING_SOURCES.md) and adjacent '
        'JSON for every input hash and per-codon exact-fragment support.\n')
    (FOCUSED / (pdf.stem + '.json')).write_text(json.dumps(value, indent=2) + '\n')
    with (FOCUSED / 'candidate_proteins.fasta').open('a') as stream:
        stream.write('>HCC1395_EIF3K_CYP39A1 complete_coding_candidate=531aa junction_codon_1based=167 RNA_translation_observed=false\n'
                     + value['protein'] + '\n')
    shutil.copyfile(pdf, FOCUSED / pdf.name)
    print('Wrote both focused PDFs, Markdown alignments and provenance JSON.')


if __name__ == '__main__':
    main()
