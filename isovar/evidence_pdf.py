"""Vector figures shared by source audits and the VCF/BAM evidence command.

ReportLab is optional; callers import this module only when writing a PDF.
"""

from collections import Counter
from html import escape

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import (
    Flowable, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle,
)

from .evidence_report import AA_BINS, READ_BINS, aa_bin, read_bin, variant_type


NAVY = colors.HexColor("#18334b")
TEAL = colors.HexColor("#0f706d")
ORANGE = colors.HexColor("#ba5a29")
GRAY = colors.HexColor("#526373")
PALE = colors.HexColor("#eaf1f4")


class EvidenceFunnel(Flowable):
    """Nested source-record counts; widths are schematic, numbers authoritative."""

    def __init__(self, counts, width=512):
        super().__init__()
        self.counts = counts
        self.width = width
        small = counts["total"]
        self.stages = (
            (("All DNA source records", small + counts["sv_total"]),
             ("Small-variant records", small))
            if counts["sv_total"] else (("All small-variant records", small),)
        ) + (
            ("Any annotated exon involved", counts["exonic"]),
            ("Exonic with original RNA ALT", counts["exonic_alt"]),
            ("Exonic with local mutant protein", counts["exonic_product"]),
        )
        self.height = 11 + len(self.stages) * 48

    def draw(self):
        c = self.canv
        total = self.stages[0][1]
        for index, (label, value) in enumerate(self.stages):
            box_width = self.width - index * 60
            left = (self.width - box_width) / 2
            bottom = self.height - 42 - index * 48
            c.setFillColor(TEAL if index == len(self.stages) - 1
                           else NAVY if index == 0 else PALE)
            c.roundRect(left, bottom, box_width, 39, 7, fill=1, stroke=0)
            dark = index in (0, len(self.stages) - 1)
            c.setFillColor(colors.white if dark else NAVY)
            c.setFont("Helvetica-Bold", 9)
            c.drawString(left + 12, bottom + 22, label)
            c.setFont("Helvetica", 8.2)
            pct = 100 * value / total if total else 0
            c.drawString(left + 12, bottom + 9, f"{value:,}  /  {pct:.1f}% of input")
            if index < len(self.stages) - 1:
                c.setStrokeColor(GRAY)
                c.line(self.width / 2, bottom - 3, self.width / 2, bottom - 9)
        c.setFillColor(GRAY)
        c.setFont("Helvetica-Oblique", 7)
        c.drawString(0, 2, "Boxes show subset order; widths are schematic, not scaled to count.")


class BinnedBars(Flowable):
    """Compact vector histogram with count labels and optional two-class bars."""

    def __init__(self, title, bins, values, *, classes=None, width=512):
        super().__init__()
        self.title = title
        self.bins = tuple(bins)
        self.values = values
        self.classes = classes
        self.width = width
        self.height = 37 + len(self.bins) * 24 + (16 if classes else 0)

    def draw(self):
        c = self.canv
        c.setFillColor(NAVY)
        c.setFont("Helvetica-Bold", 9.5)
        c.drawString(0, self.height - 12, self.title)
        if self.classes:
            c.setFont("Helvetica", 7.5)
            for index, (name, color) in enumerate(zip(
                    self.classes, (TEAL, ORANGE, NAVY, GRAY))):
                c.setFillColor(color)
                c.drawString(4 + index * 105, self.height - 27, name)
        top = self.height - (42 if self.classes else 28)
        bar_left, bar_width = 105, self.width - 225
        maximum = max((sum(self.values.get(label, ())) if self.classes
                       else self.values.get(label, 0) for label in self.bins), default=0)
        for index, label in enumerate(self.bins):
            y = top - index * 24
            c.setFillColor(GRAY)
            c.setFont("Helvetica", 7.7)
            c.drawRightString(bar_left - 9, y + 2, label)
            c.setFillColor(PALE)
            c.roundRect(bar_left, y - 3, bar_width, 13, 3, fill=1, stroke=0)
            if self.classes:
                values = self.values.get(label, ())
                scale = bar_width / maximum if maximum else 0
                offset = 0
                for number, color in zip(values, (TEAL, ORANGE, NAVY, GRAY)):
                    c.setFillColor(color)
                    c.rect(bar_left + offset, y - 3, number * scale, 13, fill=1, stroke=0)
                    offset += number * scale
                shown = " / ".join(f"{number:,}" for number in values)
            else:
                value = self.values.get(label, 0)
                c.setFillColor(TEAL)
                c.rect(bar_left, y - 3, bar_width * value / maximum if maximum else 0,
                       13, fill=1, stroke=0)
                shown = f"{value:,}"
            c.setFillColor(NAVY)
            c.setFont("Helvetica", 7.4)
            c.drawRightString(self.width, y + 2, shown)


def read_support_figures(reads_by_class, *, width=512):
    """Distributions of ALT reads by type, with and without zero-support calls."""
    classes = tuple(name for name in ("SNV", "insertion", "deletion", "replacement")
                    if name in reads_by_class)
    counts = {name: Counter(read_bin(n) for n in reads_by_class.get(name, ()))
              for name in classes}
    figures = []
    for title, bins in (("ALT reads per variant: all calls", READ_BINS),
                        ("ALT reads per variant: positive support only", READ_BINS[1:])):
        values = {label: tuple(counts[name][label] for name in classes) for label in bins}
        figures.append(BinnedBars(title, bins, values, classes=classes, width=width))
    return figures


def altered_aa_figures(aa_counts, *, width=512):
    """Top local context per ALT-observed call, including and excluding one AA."""
    counts = Counter(aa_bin(n) for n in aa_counts)
    return [
        BinnedBars("Changed AA positions: all ALT-observed calls", AA_BINS,
                   counts, width=width),
        BinnedBars("Changed AA positions: excluding one-AA calls",
                   tuple(label for label in AA_BINS if label != "1"),
                   counts, width=width),
    ]


def write_generic_report(path, calls, annotations, summary, sources, self_scores,
                         *, catalog_assessed=False):
    """Write a portable, vector-figure PDF for one VCF and RNA alignment set."""
    from . import __version__

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="EvidenceTitle", parent=styles["Title"],
                              fontName="Helvetica-Bold", fontSize=19, leading=23,
                              textColor=NAVY, spaceAfter=10))
    styles.add(ParagraphStyle(name="EvidenceBody", parent=styles["Normal"],
                              fontSize=8.5, leading=12, textColor=NAVY,
                              spaceAfter=8))
    styles.add(ParagraphStyle(name="EvidenceSection", parent=styles["Heading2"],
                              fontSize=12, leading=15, textColor=TEAL,
                              spaceBefore=10, spaceAfter=7))
    styles.add(ParagraphStyle(name="EvidenceCell", parent=styles["Normal"],
                              fontSize=8, leading=10, textColor=NAVY))
    p = lambda value, style="EvidenceBody": Paragraph(escape(str(value)), styles[style])

    def table(headers, rows, widths):
        content = [[p(value, "EvidenceCell") for value in headers]]
        content.extend([[p(value, "EvidenceCell") for value in row] for row in rows])
        result = Table(content, colWidths=widths, repeatRows=1, hAlign="LEFT")
        result.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), PALE),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PALE]),
            ("LINEBELOW", (0, 0), (-1, 0), 1, TEAL),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 6),
            ("RIGHTPADDING", (0, 0), (-1, -1), 6),
            ("TOPPADDING", (0, 0), (-1, -1), 5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ]))
        return result

    counts = summary["counts"]
    title = f"{sources['sample_id']}: variant and RNA evidence"
    story = [Spacer(1, 8), p(title, "EvidenceTitle"),
             p(f"DNA variant set: {sources['dna_source']} | RNA alignment: {sources['rna_source']}"),
             p("Counts refer to loaded small-variant alleles and this RNA read set. "
               "Reference annotation, observed ALT reads and local Isovar protein "
               "contexts are distinct evidence layers."),
             p("Variant-to-RNA funnel", "EvidenceSection"), EvidenceFunnel(counts),
             p("Outcome counts", "EvidenceSection"),
             table(("Measure", "Count", "Interpretation"), (
                 ("Loaded small variants", f"{counts['total']:,}", "VCF alleles retained by the loader"),
                 ("Any annotated exon", f"{counts['exonic']:,}", "At least one changed base or insertion flank"),
                 ("Reference coding effect", f"{counts['coding']:,}", "Any coding-sequence-changing Varcode effect"),
                 ("RNA ALT-observed", f"{counts['alt']:,}", "Original allele-spanning read evidence"),
                 ("Local protein context", f"{counts['product']:,}", "At least one translated local window"),
             ), (190, 67, 255)), PageBreak(),
             p("RNA alternate-read support by variant type", "EvidenceSection")]
    story.extend(read_support_figures(summary["reads_by_class"]))
    story.extend((p("The first panel includes calls with zero ALT reads. The second "
                    "conditions on positive ALT support. Bars count variants, with "
                    "classes shown in legend order."), PageBreak(),
                  p("Changed amino acids in local protein contexts", "EvidenceSection")))
    story.extend(altered_aa_figures(summary["aa_counts"]))
    story.extend((p("One top-ranked local context is counted per ALT-observed call. "
                    "The second panel removes the one-residue category; calls with "
                    "no local product remain explicit. Zero denotes an empty "
                    "mutation interval such as a deletion junction. No full "
                    "frameshift tail or complete isoform is inferred."), PageBreak(),
                  p("Mutation classes and examples", "EvidenceSection")))
    story.append(table(("Category", "Count", "Scope"), (
        ("DNA-predicted frameshift", f"{counts['frameshift']:,}",
         "Top Varcode effect; RNA ALT at %s" % counts["frameshift_alt"]),
        ("ncORF coding-block overlap", f"{counts['ncorf']:,}" if catalog_assessed else
         "not assessed", "Reference candidate; mutant ncORF translation not called"),
        ("Productive SV", "not assessed", "This command accepts small-variant VCF alleles"),
    ), (190, 95, 227)))

    def example_table(heading, chosen, *, score=False):
        story.append(p(heading, "EvidenceSection"))
        if not chosen:
            story.append(p("No qualifying local product was observed or scored."))
            return
        rows = []
        for call in chosen[:3]:
            ident = call["call_id"]
            protein = summary["top"][ident]
            if score:
                novel, total = self_scores[ident]
                measure = f"{novel}/{total} novel 9-mers"
            else:
                measure = f"{int(call['alt_reads']):,} ALT reads"
            sequence = protein["protein_sequence"]
            rows.append((ident, variant_type(call), measure,
                         sequence[:29] + ("..." if len(sequence) > 29 else "")))
        story.append(table(("Call", "Type", "Support / self screen", "Local sequence"),
                           rows, (150, 45, 125, 192)))

    translated = [call for call in calls if call["call_id"] in summary["top"]]
    translated.sort(key=lambda call: (-int(call["alt_reads"]), call["call_id"]))
    example_table("Highest RNA ALT support with local products", translated)
    scored = [call for call in translated if call["call_id"] in self_scores
              and self_scores[call["call_id"]][1]]
    scored.sort(key=lambda call: (
        -self_scores[call["call_id"]][0] / self_scores[call["call_id"]][1],
        -self_scores[call["call_id"]][0],
        -int(call["alt_reads"]), call["call_id"]))
    if self_scores:
        example_table("Least reference-self-like local contexts", scored, score=True)
    else:
        story.append(p("Reference-self comparison", "EvidenceSection"))
        story.append(p("Not assessed. Supply --self-fasta to screen exact "
                       "mutation-overlapping 9-mers."))
    example_table("Highest supported non-SNV local products",
                  [call for call in translated if call["class"] != "SNV"])
    story.append(p("Self comparison, when requested, measures distinct "
                   "mutation-overlapping 9-mers absent from the supplied reference "
                   "protein FASTA. It is not patient-specific self, MHC binding, "
                   "presentation, or immunogenicity. Local products are not "
                   "full-length mutant isoforms."))
    if sources.get("self_fasta_sha256"):
        story.append(p("Reference protein FASTA: %s | SHA-256: %s" % (
            sources["self_fasta"], sources["self_fasta_sha256"])))

    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(GRAY)
        canvas.drawString(50, 30, f"ISOVAR {__version__} | {sources['sample_id']}")
        canvas.drawRightString(562, 30, f"page {doc.page}")
        canvas.restoreState()

    document = SimpleDocTemplate(str(path), pagesize=letter, leftMargin=50,
                                 rightMargin=50, topMargin=48, bottomMargin=48,
                                 title=title, author="Isovar evidence report")
    document.build(story, onFirstPage=footer, onLaterPages=footer)
