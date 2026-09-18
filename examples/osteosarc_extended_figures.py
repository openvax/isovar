"""Consistent source/mode comparisons and sequence-resolved DLG5 evidence."""

import gzip
import json
from pathlib import Path

from isovar.fusion_visualization import _canvas
from isovar.protein_comparison import save_protein_comparison
from isovar.visualization import BLUE, GRAY, ORANGE, _side_note
from tests.data.osteosarc.expansion.inventory import digest
from tests.data.osteosarc.figure_comparisons.extended_footprints import CORPUS, load
from .osteosarc_footprint_figures import save_panels

LABELS = {
    "T1-ONT-dedup": "T1 / ONT dedup", "T1-ONT-tagged": "T1 / ONT tagged", "T1-short": "T1 / Illumina bulk",
    "T2-ONT-dedup": "T2 / ONT dedup", "T2-ONT-tagged": "T2 / ONT tagged", "T2-short": "T2 / Illumina bulk",
    "T0-BostonGene": "T0 BostonGene / Illumina", "T0-Personalis": "T0 Personalis / Illumina",
    "T1-Tempus": "T1 Tempus / Illumina", "T1-PacBio": "T1 / PacBio",
    "T3-ONT-dedup": "T3 / ONT dedup", "T3-ONT-tagged": "T3 / ONT tagged",
    "T3-scRNA": "T3 / Illumina scRNA", "T3-CD45neg": "T3 CD45neg / Illumina",
}


def load_dlg5():
    manifest = json.loads((CORPUS / "dlg5-manifest.json").read_text())
    path = CORPUS / manifest["file"]
    if digest(path) != manifest["sha256"]:
        raise ValueError("DLG5 evidence checksum mismatch")
    return json.loads(gzip.decompress(path.read_bytes()))


def table_panel(title, subtitle, columns, rows, note, footer):
    """A compact scientific evidence table with caveats outside the values."""
    figure, ax = _canvas(title, subtitle, 3 + .48 * len(rows))
    figure.subplots_adjust(left=.26, right=.79, top=.80, bottom=.18)
    for x, name in enumerate(columns, 1):
        ax.text(x, len(rows) + 1, name, ha="center", va="center", fontsize=11, weight="bold", color=GRAY)
    for y, (label, values) in zip(range(len(rows), 0, -1), rows):
        ax.text(-.035, y, label, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=11)
        for x, value in enumerate(values, 1):
            ax.text(x, y, str(value), ha="center", va="center", fontsize=14, color=BLUE)
        ax.axhline(y - .44, lw=.6, color="#eeeeee")
    ax.set(xlim=(.5, len(columns) + .5), ylim=(.4, len(rows) + 1.5))
    ax.axis("off")
    _side_note(ax, note)
    figure.text(.26, .075, footer, fontsize=10, color=GRAY)
    return figure


def indel_support(data):
    entries = data["indels"]
    products = [p["source"] for p in entries[0]["products"]]
    lookup = {(e["gene"], p["source"]): p for e in entries for p in e["products"]}
    rows = []
    for sid in products:
        values = []
        for e in entries:
            counts = lookup[e["gene"], sid]["default_counts"]
            values.append("%d / %d" % (counts["alt"], counts["ref"]))
        rows.append((LABELS[sid], values))
    return table_panel("Nine candidates | Extended RNA products", "Explicit small indels: alternate / reference template IDs",
                       [e["gene"] for e in entries], rows,
                       "Default allele filters\nProducts not pooled\n\nAbsent QUAL retained\nas unknown confidence\nStrict counts in JSON\n\nNot molecule counts",
                       "Other/conflicting observations remain in evidence.json. A zero alternate count is not biological absence.")


def deletion_support(entry):
    known = {(j["start"], j["end"]) for p in entry["products"] for j in p["junctions"]
             if j["operation"] == "N" and j["annotated"] and j["start"] <= entry["start"] and j["end"] >= entry["end"]}
    common = max(known, key=lambda k: sum(j["templates"] for p in entry["products"] for j in p["junctions"]
                                        if j["operation"] == "N" and (j["start"], j["end"]) == k)) if known else None
    rows = []
    for p in entry["products"]:
        skip = next((j["templates"] for j in p["junctions"] if j["operation"] == "N" and (j["start"], j["end"]) == common), 0)
        counts = p["breakpoint_matching_templates"]
        rows.append((LABELS[p["source"]], [p["aligned_templates"]["inside"], skip if common else "-",
                                          "%d / %d" % (counts["deletion"], counts["splice"])]))
    return table_panel(entry["gene"] + " | Extended RNA footprint", "Aligned sequence, ordinary splicing and DNA-endpoint gaps are separate",
                       ["Inside DNA interval", "Annotated skip", "Endpoint D / N"], rows,
                       "Primary / MAPQ >=20\nTemplate IDs\n\nD: alignment deletion\nN: skipped region\n\nNo allele-specific\ncausality inferred",
                       "Ordinary splicing is not evidence of the DNA deletion. No concrete mutant protein is assigned here.")


def fusion_support(entry):
    rows = [(LABELS[p["source"]], [p["complete_paths"], "%d / %d" % (
                sum(any(w["minimum_base_quality"] is not None for w in r["windows"]) for r in p["paths"]),
                sum(any(w["minimum_base_quality"] is None for w in r["windows"]) for r in p["paths"])), p["cell_umi_labels"]])
            for p in entry["products"]]
    return table_panel(entry["name"].replace("--", " / ") + " | RNA paths", "Observed supplementary paths; not paired mates or alternative placements",
                       ["Complete paths", "Q10 / QUAL unknown", "CB / UMI labels"], rows,
                       "Actual partner CIGARs\nMAPQ >=20\n\nTagged and dedup\nare processing relatives\nNever add their counts\n\nFrame unresolved",
                       "No validated coding sequence. Missing partner records can prevent a path call; zero is not absence of fusion RNA.")


def dlg5_panels(data):
    dna = [p for p in data["products"] if p["source"]["product"] == "DNA"]
    rows = [(p["source"]["id"], [p["signatures"]["purple_bnd"]["q20"], len(p["split_paths"]["observed_cigar"]),
                                  p["signatures"]["dragen_fields"]["q20"]]) for p in dna]
    yield "dna-junction", table_panel("DLG5 | Sequence-resolved DNA junction", "Two callers, one observed assembled haplotype; nominal fields are not the entire haplotype",
        ["BND sequence Q20", "Observed split paths", "Nominal DEL + insert"], rows,
        "Twenty-base flanks\nPrimary / MAPQ >=20\n\nColumns overlap\nDo not sum evidence\n\nOrganoid is culture\nNormals are not\nindependent replicates",
        "Purple/ESVEE's sequence matches DRAGEN's CONTIG and DNA reads. Nearby changes distinguish it from nominal DEL fields.")
    rna = [p for p in data["products"] if p["source"]["product"] not in ("DNA", "ONT-tagged")]
    rows = [(LABELS[p["source"]["id"]], ["%d / %d" % (p["start_codon"].get("reference_start", 0), p["start_codon"].get("quality_unresolved", 0)),
                p["footprints"]["purple_bnd"]["aligned_templates"]["inside"], p["signatures"]["purple_bnd"]["q20"]]) for p in rna]
    yield "rna-start-and-junction", table_panel("DLG5 | Retained start versus mutant junction", "The listed event removes the annotated DLG5-001 start region",
        ["Start Q20 / unresolved", "Inside DNA interval", "BND sequence Q20"], rows,
        "Genomic CAT = ATG\non DLG5 minus strand\n\nUnresolved: absent/low QUAL\nInside: aligned coverage\nDifferent measurements\n\nNo mutant protein\nframe justified",
        "Retained-start RNA does not exclude a subclonal DNA deletion. Missing mutant junction does not establish biological absence.")
    figure, ax = _canvas("DLG5 | What can be reconstructed?", "DNA sequence is concrete; mutant RNA and coding frame remain unresolved", 6)
    ax.axis("off")
    items = [(3, "DNA", "Explicit BND + reference flanks; observed tumor reads", BLUE),
             (2, "RNA", "No matching junction sequence in the queried products", GRAY),
             (1, "Protein", "Annotated start is deleted; no justified alternative CDS", ORANGE)]
    for y, label, text, color in items:
        ax.text(-.035, y, label, transform=ax.get_yaxis_transform(), ha="right", va="center", weight="bold", fontsize=14, color=color)
        ax.text(0, y, text, va="center", fontsize=13, color=color)
    ax.set(xlim=(0, 1), ylim=(.4, 3.6))
    _side_note(ax, "SV locus retrieval does\nnot require Varcode\nprotein prediction\n\nSplice changes are\nhypotheses until linked\nto the mutant allele")
    figure.text(.19, .07, "Do not substitute a genomic junction or an arbitrary open reading frame for an expressed mutant protein.", fontsize=11, color=GRAY)
    yield "interpretation", figure


def generate(output):
    data = load()
    root = Path(output) / "extended-rna"
    save_panels(root / "00-overview", [("indel-support", indel_support(data))])
    for entry in data["indels"]:
        products = [dict(p, label=LABELS[p["source"]]) for p in entry["products"]]
        save_protein_comparison(products, root / entry["gene"] / "protein-alternatives")
    for entry in data["deletions"]:
        save_panels(root / entry["gene"], [("rna-footprint", deletion_support(entry))])
    for entry in data["fusions"]:
        save_panels(root / entry["name"], [("path-support", fusion_support(entry))])
    dlg5 = load_dlg5()
    save_panels(root / "DLG5" / "dna-rna", dlg5_panels(dlg5))
    (root / "evidence.json").write_text(json.dumps(data, indent=2) + "\n")
    (root / "DLG5" / "dna-rna" / "evidence.json").write_text(json.dumps(dlg5, indent=2) + "\n")
    return root
