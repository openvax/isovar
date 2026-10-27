"""Focused, source-qualified DDX17 intronic RNA evidence for the COLO829 PDF."""

import json
from math import log10

from reportlab.platypus import Flowable, PageBreak

from .build_evidence_pdfs import REPORT, NAVY, TEAL, GRAY, PALE, ORANGE, para, section, data_table


class DDX17Coverage(Flowable):
    def __init__(self, result):
        super().__init__()
        self.result = result
        self.width, self.height = 512, 210

    def draw(self):
        c, d = self.canv, self.result
        left, right = 38487800, 38492300
        x = lambda position: 49 + 452 * (position - left) / (right - left)
        c.setFont('Helvetica', 8)
        c.setFillColor(NAVY)
        c.drawString(0, 200, 'Fragment depth: locally unique placements; maximum per 25-base bin')
        for level in (0, 10, 100, 1000):
            y = 110 + 70 * log10(1 + level) / log10(3001)
            c.setStrokeColor(PALE)
            c.line(49, y, 501, y)
            c.setFillColor(GRAY)
            c.drawRightString(43, y - 2, str(level))
        for values, color in ((d['coverage_max_per_bin'], TEAL),
                              (d['alt_coverage_max_per_bin'], ORANGE)):
            c.setStrokeColor(color)
            c.setLineWidth(1)
            path = c.beginPath()
            for index, value in enumerate(values):
                px = x(left + index * 25 + 12.5)
                py = 110 + 70 * log10(1 + value) / log10(3001)
                if index == 0:
                    path.moveTo(px, py)
                else:
                    path.lineTo(px, py)
            c.drawPath(path)
        for label, y, intervals in (
                ('Coding', 78, [(38487873, 38488115), (38492056, 38492115)]),
                ('NMD', 54, [(38488026, 38488115), (38489228, 38489458), (38492056, 38492115)])):
            c.setFillColor(NAVY)
            c.drawString(0, y - 2, label)
            c.setStrokeColor(GRAY)
            c.line(x(intervals[0][0]), y, x(intervals[-1][1]), y)
            for start, end in intervals:
                c.rect(x(start), y - 4, x(end) - x(start), 8, fill=1, stroke=0)
        c.setFillColor(GRAY)
        c.drawString(x(38488200), 86, '449 fragments: ordinary exon-to-exon join')
        c.setStrokeColor(TEAL)
        c.setLineWidth(1.4)
        c.line(x(38490436), 30, x(38492055), 30)
        for pos in (38490436, 38492055):
            c.line(x(pos), 26, x(pos), 34)
        c.setFillColor(TEAL)
        c.drawRightString(501, 15, '10 fragments: intronic acceptor to coding exon')
        c.setStrokeColor(ORANGE)
        c.setDash(2, 2)
        c.line(x(38489857), 40, x(38489857), 190)
        c.setDash()
        c.setFillColor(ORANGE)
        c.drawCentredString(x(38489857), 191, 'ALT deletion')
        c.setFillColor(GRAY)
        c.drawString(49, 0, 'chr22:38,487,800')
        c.drawRightString(501, 0, '38,492,300; GRCh38; DDX17 transcription runs left')


def ddx17_story(sample):
    if sample != 'COLO829':
        return []
    d = json.loads((REPORT / 'ddx17_audit/results.json').read_text())
    story = section('DDX17: intronic ALT RNA, alternative splicing and the missing protein')
    story.append(para('<b>chr22:38489856:GA&gt;G</b> | IGV (GRCh38): chr22:38489856-38489857 | '
                      '221 original ALT fragments. The old protein run found zero annotated coding '
                      'reference contexts and returned before RNA assembly. It did not test whether '
                      'cryptic exons or retained introns existed.'))
    story.append(para('The follow-up recovers both mates for 603 selected fragments and expands to '
                      '27,309 pairs from the same SRR8615617 library. Splice-aware STAR alignment to a '
                      '50-kb DDX17 reference finds continuous intronic coverage (minimum 115 fragments), '
                      'ordinary exon joins, NMD-model joins and alternative intronic acceptors.'))
    story.append(DDX17Coverage(d))
    story.append(para('<font color="#0f706d">&#9632;</font> All recovered fragments &nbsp; '
                      '<font color="#ba5a29">&#9632;</font> Original ALT-fragment subset. '
                      'Coding and NMD boxes are Ensembl 87 exons. Coverage alone does not phase a transcript.', 'FootSmall'))
    story.append(data_table(['Evidence / GRCh38 interval', 'Result and interpretation'], [
        ('Ordinary join: chr22:38488116-38492055',
         '449 locally unique fragments; this splice removes the mutation-containing intron.'),
        ('Alternative join: chr22:38490437-38492055',
         '10 locally unique fragments; canonical minus-strand splice motif; no direct ALT-fragment linkage.'),
        ('ALT-linked alternative: chr22:38486940-38489622',
         'One ambiguous fragment; 9-base anchor, two local placements, noncanonical motif, neither end an annotated exon.'),
        ('All 221 original ALT fragments',
         '214 locally unique, 5 multiply placed, 2 unaligned in the expanded local run; no locally unique ALT-fragment splice.'),
    ], [221, 291]))
    story.append(para('<b>Why no defensible mutant protein yet?</b> The deletion is intronic in both '
                      'complete coding models; a retained-intron annotation covers it but lacks a complete CDS. '
                      'Continuing either annotated coding frame into the retained intron predicts a stop at '
                      '<b>chr22:38492042-38492044</b>. The 10-fragment alternative path predicts a stop at '
                      '<b>chr22:38490360-38490362</b>. Both stops precede the deletion in the minus-strand '
                      'translation direction. These are reference-path predictions, not ALT-phased proteins; '
                      'cryptic initiation and other paths remain unresolved.'))
    story.append(para('Scope: SRA Toolkit 3.4.1; STAR 2.7.11b without supplied splice annotations; '
                      'GRCh38 local sequence verified against all four Ensembl 87 transcript sequences. '
                      'The regional SRA index uses GRCh37; recovered reads are realigned to GRCh38. '
                      'Local uniqueness is not genome-wide uniqueness. This bounded selection cannot establish '
                      'whole-library absence or mutation causality. Original pairs, every alternative placement, '
                      'source coordinates, hashes and commands: ddx17_audit/manifest.json and README.md. '
                      'Reusable intronic assessment remains tracked in '
                      '<link href="https://github.com/openvax/isovar/issues/583">#583</link>.', 'FootSmall'))
    story.append(PageBreak())
    return story
