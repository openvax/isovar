"""Named, source-separated antigen evidence with producer RNA count units."""

from html import escape
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import Paragraph, SimpleDocTemplate, Table, TableStyle

from .evidence_pdf import GRAY, NAVY, PALE, TEAL, read_support_figures, viral_support_figures

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
    """Render every supplied export without pooling sources or truncating tables."""
    from . import __version__

    styles = getSampleStyleSheet()
    for name, parent, size, leading, color in (
            ("Title", "Title", 18, 22, NAVY), ("Body", "Normal", 8.5, 12, NAVY),
            ("Section", "Heading2", 12, 15, TEAL), ("Source", "Heading3", 9, 12, NAVY),
            ("Cell", "Normal", 7.5, 10, NAVY)):
        styles.add(ParagraphStyle(name="Antigen" + name, parent=styles[parent],
            fontSize=size, leading=leading, textColor=color,
            spaceAfter=7 if name != "Cell" else 0, keepWithNext=name in {"Section", "Source"}))

    def p(value, style="Body"):
        return Paragraph(escape(str(value)), styles["Antigen" + style])

    def shown(value):
        return "-" if value is None or value == "" else str(value)

    def words(value):
        return shown(value).replace("_", " ")

    def table(headers, rows, widths):
        values = [[p(value, "Cell") for value in headers]]
        values.extend([[p(shown(value), "Cell") for value in row] for row in rows])
        result = Table(values, colWidths=widths, repeatRows=1, hAlign="LEFT", splitInRow=1)
        result.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), PALE),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PALE]),
            ("LINEBELOW", (0, 0), (-1, 0), .8, TEAL),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 5), ("RIGHTPADDING", (0, 0), (-1, -1), 5),
            ("TOPPADDING", (0, 0), (-1, -1), 4), ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ]))
        return result

    categories = inventory["categories"]
    story = [p(f"{inventory['sample_id']}: antigen-source evidence", "Title"),
             p("Original reads are sequenced segments; fragments group mates. Each measurement names its "
               "subject: an ALT allele, a junction, RNA path contributors or a complete interval. "
               "Alternative paths can share reads and fragments; their counts are not additive. "
               "Unique means unique among retained paths under the producer's placement rules. "
               "A dash means unavailable or unmeasured; zero means assessed without support."),
             table(("Category", "Status", "Exports", "Input unit"),
                   [(LABELS[key], words(category["status"]), len(category["records"]), UNITS[key])
                    for key, category in categories.items()], (150, 95, 55, 212)),
             p("RNA and coding hypotheses do not establish protein expression or peptide presentation. "
               "Input samples and libraries retain their own counts within this report group.")]

    def source(item, label=None):
        identifier = label or item.get("source") or Path(item.get("source_file", "export")).name
        parts = [str(item.get("input_sample_id", inventory["sample_id"])), str(identifier)]
        filename = Path(item.get("source_file", "")).name
        if filename and filename != identifier:
            parts.append(filename)
        story.append(p(" | ".join(parts), "Source"))
        if item.get("discovery_truncated"):
            story.append(p("Discovery reached a configured limit; undiscovered alternative paths may remain."))

    def measurements(item):
        rows = item.get("measurements", [])
        if not rows:
            story.append(p("No path-level RNA measurements were supplied."))
            return
        scopes = sorted({words(row.get("scope")) for row in rows})
        story.append(p("Support scope: " + "; ".join(scopes) + "."))
        story.append(table(("Named RNA path", "Measurement", "Reads", "Fragments", "Status"),
            [(row["name"], row["subject"], row.get("reads"), row.get("fragments"), words(row.get("status")))
             for row in rows], (166, 137, 43, 54, 112)))

    small = categories["small_variants"]
    if small["status"] == "assessed":
        story.append(p("Small variants: named calls and RNA ALT support", "Section"))
        for item in small["records"]:
            source(item)
            ncorf = item["ncorf_reference_candidates"]
            ncorf_text = "not assessed" if ncorf is None else f"{ncorf:,} candidates"
            story.append(p(f"{item['input_records']:,} loaded literal alleles; {item['exon_involved']:,} exon-involved; "
                f"{item['rna_alt_observed']:,} with RNA ALT support; {item['rna_unmeasured']:,} unmeasured. "
                f"Reference ncORF coding-block overlap: {ncorf_text}."))
            rows = []
            for call in item.get("calls", []):
                name = " | ".join(str(value) for value in (call.get("gene"), call.get("top_effect"), call["call_id"]) if value)
                rows.append((name, shown(call.get("alt_reads")),
                    shown(call.get("alt_fragments")) + " / " + shown(call.get("total_fragments")),
                    shown(call.get("ref_reads")) + " / " + shown(call.get("other_reads")),
                    words(call.get("rna_status")) + "; " + shown(call.get("protein_count")) + " local proteins"))
            if rows:
                story.append(table(("Gene / effect / call", "ALT reads", "ALT / total fragments", "REF / other reads", "RNA / protein outcome"),
                                   rows, (170, 50, 80, 75, 137)))
            if len(rows) > 12 or not rows:
                story.extend(read_support_figures(item["exonic_reads_by_class"], scope="exon-involved variant"))

    if categories["sv"]["records"] or categories["fusion"]["records"]:
        story.append(p("Structural variants and supplied fusions", "Section"))
        for key in ("sv", "fusion"):
            for item in categories[key]["records"]:
                source(item, item["event_id"] + (" / " + item["source"] if item.get("source") else ""))
                story.append(p(LABELS[key] + ": " + words(item["status"]) + "."))
                if item.get("reasons"):
                    story.append(p("Reasons: " + "; ".join(item["reasons"])))
                if key == "sv":
                    story.append(p(f"{item['event_linked_translated_paths']} event-linked translated paths; "
                        f"{item['exploratory_orf_candidates']} exploratory ORF candidates. "
                        "Regional or annotated junctions remain context and do not establish support for the nominated SV."))
                measurements(item)

    if categories["viral"]["records"]:
        story.append(p("Viral contigs and RNA paths", "Section"))
        for item in categories["viral"]["records"]:
            source(item)
            story.append(p(f"{item['observed_paths']} observed paths; {item['coding_predictions']} coding predictions "
                           f"({item['complete_cds_predictions']} complete CDS)."))
            story.append(table(("Viral contig", "Aligned reads", "Aligned fragments", "Sequence-assigned fragments", "Status"),
                [(row["contig"], row.get("aligned_segments"), row["aligned_fragments"],
                  row["sequence_assigned_fragments"], words(row["status"])) for row in item["contigs"]],
                (110, 65, 75, 95, 167)))
            measurements(item)
            if len(item["contigs"]) > 12:
                story.extend(viral_support_figures(item["contigs"]))

    if categories["isoforms"]["records"]:
        story.append(p("Targeted antigen-gene isoforms", "Section"))
        for item in categories["isoforms"]["records"]:
            source(item)
            targets = item.get("targets", [])
            story.append(p("Targets: " + (", ".join(row.get("gene_name", row["gene_id"]) for row in targets) or "not named in export") + ". "
                f"{item['rna_supported_genes']} RNA-supported genes; "
                f"{item['fully_covered_protein_groups']} fully covered protein groups."))
            measurements(item)
            story.append(p("Whole retained RNA paths need not be full-length transcripts. Whole-CDS counts require "
                           "a complete coding model and original spanning witnesses; tiled coverage alone is insufficient."))

    story.append(p("Input provenance", "Section"))
    story.append(p("Input samples, sources and export filenames are named above. The accompanying JSON retains "
                   "their full source paths and detailed measurements. No library is pooled into another source's RNA denominator."))

    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(GRAY)
        canvas.drawString(50, 30, f"ISOVAR {__version__} | {inventory['sample_id']}")
        canvas.drawRightString(562, 30, f"page {doc.page}")
        canvas.restoreState()

    document = SimpleDocTemplate(str(path), pagesize=letter, leftMargin=50, rightMargin=50,
        topMargin=48, bottomMargin=48, title=f"{inventory['sample_id']}: antigen-source evidence",
        author="Isovar antigen report", invariant=1)
    document.build(story, onFirstPage=footer, onLaterPages=footer)
