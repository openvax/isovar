"""Vector report for source-separated Isovar antigen evidence exports."""

from html import escape

from reportlab import rl_config
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import PageBreak, Paragraph, SimpleDocTemplate, Table, TableStyle

from .evidence_pdf import (
    GRAY, NAVY, PALE, TEAL, read_support_figures, viral_support_figures,
)


rl_config.invariant = 1
LABELS = {
    "small_variants": "Literal DNA variants",
    "sv": "Nominated DNA SVs",
    "fusion": "Supplied fusion RNAs",
    "viral": "Viral RNA panels",
    "isoforms": "Targeted antigen genes",
}
UNITS = {
    "small_variants": "VCF alleles",
    "sv": "DNA events",
    "fusion": "supplied transcripts",
    "viral": "selected contigs per export",
    "isoforms": "target genes per export",
}


def write_antigen_report(path, inventory):
    """Render one PDF; source counts stay in their native units."""
    from . import __version__

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="AntigenTitle", parent=styles["Title"],
                              fontName="Helvetica-Bold", fontSize=18, leading=22,
                              textColor=NAVY, spaceAfter=10))
    styles.add(ParagraphStyle(name="AntigenBody", parent=styles["Normal"],
                              fontSize=8.5, leading=12, textColor=NAVY,
                              spaceAfter=8))
    styles.add(ParagraphStyle(name="AntigenSection", parent=styles["Heading2"],
                              fontSize=12, leading=15, textColor=TEAL,
                              spaceBefore=10, spaceAfter=7))
    styles.add(ParagraphStyle(name="AntigenCell", parent=styles["Normal"],
                              fontSize=8, leading=10, textColor=NAVY))
    styles.add(ParagraphStyle(name="AntigenHeader", parent=styles["AntigenCell"],
                              fontName="Helvetica-Bold", textColor=colors.white))

    def p(value, style="AntigenBody"):
        return Paragraph(escape(str(value)), styles[style])

    def table(headers, rows, widths):
        values = [[p(x, "AntigenHeader") for x in headers]]
        values.extend([[p(x, "AntigenCell") for x in row] for row in rows])
        result = Table(values, colWidths=widths, repeatRows=1, hAlign="LEFT")
        result.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), NAVY),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PALE]),
            ("LINEBELOW", (0, 0), (-1, 0), 1, TEAL),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 6),
            ("RIGHTPADDING", (0, 0), (-1, -1), 6),
            ("TOPPADDING", (0, 0), (-1, -1), 5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ]))
        return result

    categories = inventory["categories"]
    overview = []
    for key, category in categories.items():
        count = len(category["records"])
        overview.append((LABELS[key], category["status"],
                         (f"{count} export(s)" if count else "—"), UNITS[key]))
    story = [p(f"{inventory['sample_id']}: antigen-source evidence", "AntigenTitle"),
             p("Each row has its own input unit. The report does not sum DNA "
               "alleles, SV events, fusion transcripts, viral contigs or target genes "
               "into one antigen count. A missing export is not assessed."),
             p("Source coverage", "AntigenSection"),
             table(("Category", "Status", "Inputs", "Denominator"),
                   overview, (157, 75, 85, 195)),
             p("Interpretation", "AntigenSection"),
             p("RNA sequence and coding hypotheses are candidate evidence, not "
               "measured protein expression or peptide presentation. Reference "
               "ncORF overlaps, exploratory SV ORFs and viral coding predictions "
               "remain separate from locally translated small-variant products.")]

    def new_page(title):
        story.extend((PageBreak(), p(title, "AntigenSection")))

    small = categories["small_variants"]
    if small["status"] == "assessed":
        item = small["records"][0]
        ncorf = item["ncorf_reference_candidates"]
        ncorf_text = "not assessed" if ncorf is None else f"{ncorf:,} candidates"
        new_page("Exon-involved small variants: original RNA ALT reads")
        story.append(p(
            f"{item['input_records']:,} loaded literal VCF alleles; "
            f"{item['exon_involved']:,} involve an annotated exon; "
            f"{item['rna_alt_observed']:,} have an original RNA ALT read; "
            f"{item['local_products']:,} have a local mutant protein context. "
            f"Reference ncORF coding-block overlap: {ncorf_text}. "
            "The plots include zero-support exon-involved calls "
            "and then condition on positive ALT support."))
        story.extend(read_support_figures(item["exonic_reads_by_class"],
                                          scope="exon-involved variant"))
        story.append(p("Insertion flanks count as candidate exon involvement. "
                       "Symbolic SVs use their own event and junction units."))

    sv = categories["sv"]
    fusion = categories["fusion"]
    if sv["status"] == "assessed" or fusion["status"] == "assessed":
        new_page("Structural variants and supplied fusions")
        if sv["records"]:
            story.append(p("Nominated DNA SV events", "AntigenSection"))
            rows = sorted(sv["records"], key=lambda row: (
                -row["max_direct_junction_fragments"], row["event_id"]))
            story.append(table(("Event", "RNA status", "Direct junction fragments",
                                "Linked translated paths"),
                               ((row["event_id"], row["status"],
                                 row["max_direct_junction_fragments"],
                                 row["event_linked_translated_paths"])
                                for row in rows[:12]), (145, 127, 116, 124)))
            story.append(p(
                f"{sum(row['exploratory_orf_candidates'] for row in rows):,} "
                "exploratory ORF candidates across supplied SV exports. These "
                "may be partial or lack full-interval read support. Event linkage, "
                "junction support and translation are separate evidence layers."))
        if fusion["records"]:
            story.append(p("Supplied fusion RNAs", "AntigenSection"))
            rows = sorted(fusion["records"], key=lambda row: (
                -row["direct_junction_fragments"], row["event_id"]))
            story.append(table(("Fusion", "Status", "Direct fragments",
                                "Translated candidates"),
                               ((row["event_id"], row["status"],
                                 row["direct_junction_fragments"],
                                 row["translated_candidates"])
                                for row in rows[:12]), (150, 135, 110, 117)))

    for number, row in enumerate(categories["viral"]["records"], 1):
        new_page(f"Viral RNA panel {number}: aligned fragments")
        story.append(p(
            f"{row['selected_contigs']:,} selected contigs; "
            f"{row['aligned_contigs']:,} with aligned fragments; "
            f"{row['assigned_contigs']:,} with sequence-assigned fragments; "
            f"{row['observed_paths']:,} observed paths; "
            f"{row['coding_predictions']:,} predicted coding paths "
            f"({row['complete_cds_predictions']:,} complete CDS). "
            "Bars count contigs, not fragments or viral genomes."))
        figures = viral_support_figures(row["contigs"])
        story.extend(figures[:2])
        new_page(f"Viral RNA panel {number}: sequence-assigned fragments")
        story.extend(figures[2:])
        story.append(p("The second panel excludes zero-support contigs. "
                       "Fragment support can be ambiguous between viral paths; "
                       "coding predictions do not establish translated proteins."))

    if categories["isoforms"]["records"]:
        new_page("Targeted antigen-gene isoforms")
        story.append(table(("Source", "Target genes", "RNA-supported genes",
                            "Fully covered protein groups"),
                           ((row["source"], row["target_genes"],
                             row["rna_supported_genes"],
                             row["fully_covered_protein_groups"])
                            for row in categories["isoforms"]["records"]),
                           (166, 88, 115, 143)))
        story.append(p("A selected target panel is not an exhaustive cancer-testis "
                       "antigen catalogue. RNA coverage of a coding sequence "
                       "supports a sequence hypothesis, not its translation."))

    new_page("Input provenance and scope")
    sources = [(LABELS[key], row["source_file"])
               for key, category in categories.items()
               for row in category["records"] if "source_file" in row]
    if small["status"] == "assessed":
        sources.insert(0, (LABELS["small_variants"], small["source"]))
    story.append(table(("Category", "Input export"), sources, (153, 359)))
    story.append(p("The JSON alongside this PDF retains every supplied export's "
                   "separate summary and input path. Unprovided categories are "
                   "explicitly not assessed."))

    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(GRAY)
        canvas.drawString(50, 30, f"ISOVAR {__version__} | {inventory['sample_id']}")
        canvas.drawRightString(562, 30, f"page {doc.page}")
        canvas.restoreState()

    document = SimpleDocTemplate(str(path), pagesize=letter, leftMargin=50,
                                 rightMargin=50, topMargin=48, bottomMargin=48,
                                 title=f"{inventory['sample_id']} antigen-source evidence",
                                 author="Isovar antigen evidence report")
    document.build(story, onFirstPage=footer, onLaterPages=footer)
