"""Paginated RNA protein alternatives across explicitly labelled source products.

Products are displayed separately, never pooled or treated as molecular
replicates. A frame is reference-supported relative to the nominated edit;
this is not unrestricted six-frame ORF prediction.
"""

import json
from pathlib import Path

from .default_parameters import PLOT_DPI, PLOT_PROTEIN_ROWS_PER_PAGE, PLOT_WIDTH
from .visualization import (
    BLUE, GRAY, INK, draw_protein_rows, _plot_imports,
    _protein_disagreements, _style_axis,
)


def protein_groups(proteins):
    """Compact exact covered subsequences, never synthesize a longer protein.

    Frame/transcript assignments must agree. A shorter ambiguous result can
    belong to several displayed alternatives; its counts are never added to
    any representative. Every original rank remains explicitly recoverable.
    """
    def covers(a, b):
        def contexts(p):
            return {(f["strand"], f["variant_codon_phase"], tid)
                    for f in p["frames"] for tid in f["transcript_ids"]}
        if contexts(a) != contexts(b) or a["frameshift"] != b["frameshift"]:
            return False
        left = a["mutation_start"] - b["mutation_start"]
        sa = a["amino_acids"] + ("*" if a["ends_with_stop_codon"] else "")
        sb = b["amino_acids"] + ("*" if b["ends_with_stop_codon"] else "")
        return left >= 0 and sa[left:left + len(sb)] == sb

    representatives = [i for i, p in enumerate(proteins) if not any(
        j != i and covers(other, p) and (not covers(p, other) or j < i)
        for j, other in enumerate(proteins))]
    return [(i + 1, [j + 1 for j, p in enumerate(proteins) if covers(proteins[i], p)])
            for i in representatives]


def comparison_rows(products):
    """Return every recovered protein plus unavailable RNA/Varcode outcomes.

    Parameters
    ----------
    products : sequence of dict
        Each contains unique ``source``, a display ``label`` and
        ``visualization`` from collect_visualization_data with an uncapped
        creator. Products must describe the same variant and annotation.
    """
    if not products or len({p["source"] for p in products}) != len(products):
        raise ValueError("Expected nonempty, uniquely identified source products")
    first = products[0]["visualization"]
    rows, predictions, transcripts = [], {}, {}
    for product in products:
        data = product["visualization"]
        if data["variant"] != first["variant"] or data["varcode_version"] != first["varcode_version"]:
            raise ValueError("Cannot compare different variants, assemblies or Varcode versions")
        for model in data["transcripts"]:
            reference = {k: model[k] for k in ("id", "strand", "exons")}
            if model["id"] in transcripts and transcripts[model["id"]] != reference:
                raise ValueError("Transcript models differ between source products")
            transcripts[model["id"]] = reference
        for mode in sorted(data["modes"], key=lambda m: m["assembly"]):
            if "proteins" not in mode or mode["settings"]["max_protein_sequences_per_variant"]:
                raise ValueError("All-alternative comparison requires uncapped protein collection")
            groups = protein_groups(mode["proteins"])
            for rank, covered in groups or [(None, [])]:
                protein = mode["proteins"][rank - 1] if rank else None
                frames = sorted({(f["strand"], f["variant_codon_phase"]) for f in protein["frames"]}) if protein else []
                status = ("reconstructed" if protein else "no_callable_reads" if not
                          sum(c["templates"] for c in data["counts"].values()) else
                          "no_alternate_support" if not data["counts"]["alt"]["templates"] else "no_translated_protein")
                rows.append(dict(source=product["source"], label=product["label"], assembly=mode["assembly"],
                                 rank=rank if protein else None, protein=protein, frames=frames, status=status,
                                 covered_ranks=covered,
                                 transcript_ids=protein["transcript_ids"] if protein else []))
                quality = data.get("provenance", {}).get("input_quality", {})
                rows[-1]["input_quality"] = quality
        for prediction in data["reference_predictions"]:
            key = prediction["transcript_id"]
            if key in predictions and predictions[key] != prediction:
                raise ValueError("Reference predictions differ between source products")
            predictions[key] = prediction
    # Group identical predictions for display, but keep every transcript and
    # missing-prediction reason. Stops/frame semantics are part of identity.
    grouped = {}
    for prediction in predictions.values():
        p = prediction["protein"]
        identity = {k: p[k] for k in ("amino_acids", "mutation_start", "mutation_end", "frameshift", "ends_with_stop_codon")} if p else None
        key = json.dumps((identity, prediction.get("unavailable_reason")), sort_keys=True)
        grouped.setdefault(key, []).append(prediction)
    for rank, group in enumerate(grouped.values(), 1):
        rows.append(dict(source="Varcode", label="Varcode", assembly=None, rank=rank,
                         protein=group[0]["protein"], frames=[],
                         status="predicted" if group[0]["protein"] else "no_concrete_prediction",
                         transcript_ids=sorted(p["transcript_id"] for p in group), predictions=group))
    if not predictions:
        rows.append(dict(source="Varcode", label="Varcode", assembly=None, rank=None,
                         protein=None, frames=[], transcript_ids=[], predictions=[], status="no_reference_prediction_supplied"))
    return rows


def protein_comparison_figures(products, rows_per_page=PLOT_PROTEIN_ROWS_PER_PAGE):
    """Yield white-background pages without silently dropping alternatives."""
    if isinstance(rows_per_page, bool) or not isinstance(rows_per_page, int) or rows_per_page < 2:
        raise ValueError("rows_per_page must be an integer >= 2")
    rows = comparison_rows(products)
    Figure, rc_context, rectangle = _plot_imports()
    differences = _protein_disagreements(row["protein"] for row in rows)
    title = " / ".join(products[0]["visualization"]["genes"]) + " | Protein alternatives"
    references = [r for r in rows if r["source"] == "Varcode"]
    rna = [r for r in rows if r["source"] != "Varcode"]
    reference_size = min(len(references), max(1, rows_per_page // 2))
    rna_size = rows_per_page - reference_size
    pages = [(rna[start:start + rna_size] + references[ref:ref + reference_size])
             for start in range(0, len(rna), rna_size)
             for ref in range(0, len(references), reference_size)]
    for number, page in enumerate(pages, 1):
        with rc_context({"svg.fonttype": "none", "font.family": "DejaVu Sans"}):
            figure = Figure(figsize=(PLOT_WIDTH, max(6.5, 2.5 + .8 * len(page))), facecolor="white")
            ax = figure.subplots()
            _style_axis(ax, "Sample / technology / assembly mode; all recovered alternatives")
            display = []
            for row in page:
                p = row["protein"]
                label = row["label"]
                if row["assembly"] is not None:
                    label += "\nAssembly " + ("on" if row["assembly"] else "off")
                    if row["rank"]:
                        label += " / rank %d" % row["rank"]
                elif row["rank"]:
                    label += " %d" % row["rank"]
                if row["source"] == "Varcode":
                    note = "%d transcript(s)\n%s" % (len(row["transcript_ids"]),
                                                       "Reference + edit" if p else "No concrete prediction" if row["predictions"] else "No prediction supplied")
                elif p:
                    phases = ", ".join("%s/%d" % (s, phase) for s, phase in row["frames"])
                    note = "%d templates\nStrand/phase: %s" % (p["templates"], phases)
                    if row.get("input_quality", {}).get("primary_records_missing_qualities"):
                        note += "\nInput includes missing QUAL"
                else:
                    note = row["status"].replace("_", " ")
                    if row.get("input_quality", {}).get("primary_records_missing_qualities"):
                        note += "\nMissing QUAL in input"
                display.append((label, p, INK if row["source"] == "Varcode" else BLUE if row["assembly"] else GRAY, note))
            draw_protein_rows(ax, display, rectangle, differences)
            figure.suptitle(title, x=.035, y=.975, ha="left", fontsize=19, fontweight="bold")
            figure.subplots_adjust(left=.24, right=.79, bottom=.16, top=.85)
            figure.text(.24, .065, "Orange: nominated mutation. Magenta: differing residues at shared offsets; missing context is not a mismatch.", fontsize=10, color=GRAY)
            figure.text(.24, .03, "Page %d/%d. Varcode repeats for comparison. Exact shorter contexts grouped, not joined. All ranks/frames: evidence.json." %
                        (number, len(pages)), fontsize=9, color=GRAY)
        yield "protein-comparison-%02d" % number, figure


def save_protein_comparison(products, output_dir, dpi=PLOT_DPI, rows_per_page=PLOT_PROTEIN_ROWS_PER_PAGE):
    """Save all comparison pages, source-labelled evidence and transcript names."""
    from matplotlib.backends.backend_pdf import PdfPages
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=False)
    rows = comparison_rows(products)
    (directory / "evidence.json").write_text(json.dumps(dict(products=products, rows=rows), indent=2) + "\n")
    with PdfPages(directory / "all-figures.pdf") as pdf:
        for name, figure in protein_comparison_figures(products, rows_per_page):
            pdf.savefig(figure)
            for extension in ("png", "svg"):
                figure.savefig(directory / (name + "." + extension), dpi=dpi, facecolor="white")
            figure.clear()
    return directory
