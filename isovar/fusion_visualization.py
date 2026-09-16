"""Separate, publication-ready panels for an explicit fusion result."""
import json
from dataclasses import asdict
from pathlib import Path
import re
import textwrap

from .default_parameters import PLOT_DPI, PLOT_WIDTH
from .visualization import BLUE, GRAY, INK, ORANGE, _plot_imports, _side_note, _style_axis

GREEN = "#009E73"


def _canvas(title, subtitle, height=5.3):
    Figure, _, _ = _plot_imports()
    figure = Figure(figsize=(PLOT_WIDTH, height), facecolor="white")
    ax = figure.subplots()
    figure.suptitle(title, x=.035, y=.965, ha="left", fontsize=19, fontweight="bold")
    _style_axis(ax, subtitle)
    figure.subplots_adjust(left=.19, right=.78, bottom=.23, top=.76)
    return figure, ax


def fusion_figures(result, references=(), transcript_names=None):
    """Yield named panels; unresolved RNA never receives an invented protein track.

    Transcript panels are local genomic windows, not full-length isoform models.
    A straight thin line denotes the transcript locus; only annotated exons
    are thick. No inferred splice connector is drawn across a fusion.
    """
    _, _, rectangle = _plot_imports()
    references = tuple(references)
    names = transcript_names or {}
    j0, j1 = result["junction_interval"]
    sequence = result["cdna_sequence"]
    title = result["event_id"].replace("--", " / ") + "   " + result["provenance"]["sample_id"]
    figure, ax = _canvas(title, "Directly observed fusion RNA")
    observations = result["evidence"]["observations"]
    groups = {}
    for read in observations:
        key = (read["sample_id"], read["library_id"], read["fragment_id"])
        groups.setdefault(key, []).append(read)
    shown = list(groups.values())[:6]
    for y, rows in enumerate(shown[::-1], 1):
        for read in rows:
            start, end = read["cdna_start"], read["cdna_start"] + len(read["sequence"])
            for lo, hi, color in [(start, min(end, j0), BLUE), (max(start, j1), end, GREEN)]:
                if hi > lo:
                    ax.plot([lo, hi], [y, y], color=color, lw=5, solid_capstyle="butt")
            if j1 > j0:
                ax.plot([max(start, j0), min(end, j1)], [y, y], color=ORANGE, lw=5, solid_capstyle="butt")
        ax.text(-.025, y, "Fragment %d" % (len(shown)-y+1), transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=11)
    for boundary in sorted({j0, j1}):
        ax.axvline(boundary, color=ORANGE, lw=1.2, ls="--")
    ax.set(xlim=(-2, len(sequence)+2), ylim=(.3, max(2, len(shown)+.7)),
           xlabel="Supplied cDNA offset (nt, 5' to 3')")
    status = result["status"].replace("_", " ")
    note = ("%d spanning fragments\n%d shown\n\n%s\n\nBlue: donor RNA\nGreen: acceptor RNA" %
            (result["evidence"]["directly_spanning_fragments"], len(shown), status))
    if j1 > j0:
        note += "\nOrange: %d-nt insert" % (j1-j0)
    _side_note(ax, note)
    motif = sequence[max(0,j0-18):j0] + " | " + (sequence[j0:j1] + " | " if j1 > j0 else "") + sequence[j1:j1+18]
    figure.text(.19, .075, motif, fontfamily="monospace", fontsize=12, color=INK)
    figure.text(.19, .035, "Observed junction sequence; counts describe this supplied hypothesis, not total event abundance.",
                fontsize=10, color=GRAY)
    yield "junction", figure

    for side, color in (("donor", BLUE), ("acceptor", GREEN)):
        partner = result[side]
        position = partner["position"]
        models = [r for r in references if r.contig == partner["contig"] and r.reference_name == result["reference_name"]
                  and r.exons[0][0] < position+120 and r.exons[-1][1] > position-120]
        compatible = set(result["compatible_transcripts"][side])
        # Relevant same-strand models first; preserve all models in the JSON.
        models.sort(key=lambda r: (r.transcript_id not in compatible, r.strand != partner["strand"], r.transcript_id))
        shown_models = models[:6]
        figure, ax = _canvas(title, side.capitalize() + " reference context", 3.7 + .55*len(shown_models))
        lo, hi = position-120, position+120
        for y, ref in enumerate(shown_models[::-1], 1):
            a, b = max(lo, ref.exons[0][0]), min(hi, ref.exons[-1][1])
            if b > a:
                ax.plot([a-position,b-position], [y,y], color=GRAY, lw=.9)
            for start, end in ref.exons:
                start, end = max(start,lo), min(end,hi)
                if end > start:
                    ax.add_patch(rectangle((start-position,y-.13),end-start,.26,
                        facecolor=color if ref.transcript_id in compatible else GRAY, edgecolor="none"))
            label = ref.transcript_id + ("\n(" + names[ref.transcript_id] + ")" if names.get(ref.transcript_id) else "")
            ax.text(-.035,y,label, transform=ax.get_yaxis_transform(),ha="right",va="center",fontsize=10)
            ax.text(122,y,ref.strand,va="center",color=GRAY)
        ax.axvline(0,color=ORANGE,lw=1.3,ls="--")
        ax.set(xlim=(-120,120), ylim=(.3,max(2,len(shown_models)+.7)),
               xlabel="Genomic offset from %s:%s (interbase; forward strand)" % (partner["contig"],position))
        _side_note(ax, "%d of %d models\n%s event strand\n\nThick: annotated exon\nThin: transcript locus\n\n%s" %
                   (len(shown_models),len(models),partner["strand"],status))
        reasons = "; ".join(result["reasons"])
        if reasons:
            # Raw identifiers and all reasons remain in evidence.json.
            summary = ("Junction lies before the annotated donor CDS; no coding frame assigned."
                       if any("junction_before_donor_CDS" in r for r in result["reasons"]) else
                       "No validated coding frame for this observed RNA window; no protein is inferred.")
            figure.text(.19,.065,summary,fontsize=11,color=INK)
        yield side + "-context", figure

    for number, protein in enumerate(result["translations"],1):
        figure, ax = _canvas(title, "Protein hypothesis %d" % number)
        left = max(0,protein["junction_in_translated_cds"][0]//3-24)
        right = min(len(protein["amino_acids"]),left+49)
        amino_acids = protein["amino_acids"][left:right]
        for i, aa in enumerate(amino_acids,left):
            ax.text(i+.5,1,aa,ha="center",va="center",fontsize=14,fontfamily="monospace")
        for boundary in sorted(set(protein["junction_in_translated_cds"])):
            if left <= boundary/3 <= right:
                ax.axvline(boundary/3,color=ORANGE,ls="--",lw=1.3)
        ax.set(xlim=(left,max(left+1,right)),ylim=(0,2),xlabel="Protein offset (amino acids; local donor-junction window)")
        _side_note(ax, "%s\n%s\n%d junction peptides\n%s" %
                   (status, "Observed CDS start" if protein["complete_5prime"] else "Conditional partial CDS",
                    len(protein["junction_peptides"]), "Ends at stop" if protein["ends_with_stop_codon"] else "Sequence ends first"))
        figure.text(.19,.065,"No proteome-novelty or protein-expression claim. Alternative hypotheses are not ranked.",fontsize=10,color=GRAY)
        yield "protein-%d" % number, figure


def save_fusion_figures(result, output_dir, references=(), transcript_names=None, dpi=PLOT_DPI):
    """Save individual PNG/SVG panels and one vector PDF, without overwriting."""
    if isinstance(dpi,bool) or not isinstance(dpi,int) or dpi < 72:
        raise ValueError("dpi must be an integer >= 72")
    references = tuple(references)
    _, rc_context, _ = _plot_imports()
    from matplotlib.backends.backend_pdf import PdfPages

    slug = re.sub(r"[^A-Za-z0-9_.-]", "_", result["event_id"] + "-" + result["provenance"]["sample_id"])
    directory = Path(output_dir) / slug
    directory.mkdir(parents=True,exist_ok=False)
    with rc_context({"font.family":"DejaVu Sans","svg.fonttype":"none","pdf.fonttype":42}), \
            PdfPages(directory / "all-figures.pdf") as pdf:
        for name,figure in fusion_figures(result,references,transcript_names):
            figure.savefig(directory / (name+".png"),dpi=dpi,facecolor="white")
            figure.savefig(directory / (name+".svg"),facecolor="white",metadata={"Date":None})
            pdf.savefig(figure,facecolor="white")
            figure.clear()
    (directory / "evidence.json").write_text(json.dumps(result,indent=2)+"\n")
    (directory / "reference-models.json").write_text(json.dumps(dict(
        references=[asdict(r) for r in references], transcript_names=transcript_names or {}),indent=2)+"\n")
    (directory / "caption.md").write_text(textwrap.dedent("""\
        # Fusion evidence

        Blue/green denote donor/acceptor RNA; orange marks junction boundaries
        and any inserted sequence. RNA rows represent physical fragment groups,
        not supplementary records. Only six groups/models are displayed; all
        evidence and compatible transcript identifiers remain in evidence.json.
        Reference panels show local genomic context, not a full transcript or
        inferred fusion splice model. Thin lines are transcript loci; thick
        boxes are annotated exons. The event strand is separate from each
        reference model's strand. Unresolved RNA receives no protein track.
        Conditional partial CDSs and alternative translations must not be
        silently treated as a uniquely expressed coding fusion.
        """))
    return directory
