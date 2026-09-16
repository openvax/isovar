"""Auditable length, phase, haplotype and supplied-fusion gallery extension."""
import argparse
from collections import Counter, defaultdict
import gzip
from hashlib import sha256
import json
import logging
from pathlib import Path

import pysam

from isovar.fusion import fusion_from_dict, reconstruct_fusion
from isovar.fusion_visualization import GREEN, _canvas, save_fusion_figures
from isovar.visualization import BLUE, GRAY, ORANGE, _plot_imports, _protein_disagreements, _side_note
from tests.data.osteosarc.expansion.references import translate
from tests.osteosarc_protein_helpers import transcript_offset
from tests.data.osteosarc.figure_comparisons import nr2f2
from . import osteosarc_assembly_figures

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "tests/data/osteosarc/figure_comparisons/corpus"
FUSIONS = ROOT / "tests/data/fusions/corpus"
EXCLUDED = 4 | 256 | 512 | 1024 | 2048


def load_context():
    manifest = json.loads((CORPUS/"context-manifest.json").read_text())
    raw = (CORPUS/manifest["file"]).read_bytes()
    if sha256(raw).hexdigest() != manifest["sha256"]:
        raise ValueError("Context fixture digest mismatch")
    return json.loads(gzip.decompress(raw))


def original_reads(source):
    names = [str(i) for i in range(1,23)] + ["X","Y","M","MT"]
    names += ["chr"+n for n in names]
    header = pysam.AlignmentHeader.from_references(names,[300000000]*len(names))
    for sam in source["original_records"]:
        read = pysam.AlignedSegment.fromstring(sam,header)
        if not read.flag & EXCLUDED and read.mapping_quality >= 1 and read.query_sequence is not None:
            yield read


def verify_context(data):
    """Recount original observations and independently translate retained RNA."""
    for source in data["map2"]["sources"]:
        if "original_records" not in source:
            continue  # All 35 audited products are retained; five have raw records here.
        states = defaultdict(set)
        for read in original_reads(source):
            coordinates = {g:q for q,g in read.get_aligned_pairs(matches_only=True)}
            if 209694763 not in coordinates or 209694805 not in coordinates:
                continue
            lo,hi = coordinates[209694763],coordinates[209694805]+1
            if read.query_qualities is None or min(read.query_qualities[lo:hi]) < 20:
                continue
            sequence = read.query_sequence[lo:hi]
            label = next((k for k,v in data["map2"]["hypotheses"].items() if sequence == v),"other")
            key = (read.get_tag("RG") if read.has_tag("RG") else "",read.query_name)
            states[key].add(label)
        counts = dict(Counter(next(iter(s)) if len(s)==1 else "conflicting_template" for s in states.values()))
        assert counts == source["counts"]["q20"]
    for source in data["znf436"]:
        alt = {r["name"] for r in source["alt_observations"]}
        molecules = {(r.get_tag("CB"),r.get_tag("UB")) for r in original_reads(source) if r.query_name in alt}
        assert len(molecules) == 3
        for mode in source["modes"]:
            assert mode["validation"]["matches_expected"] and not mode["nonfocal_indels"]
            for witness in mode["witnesses"]:
                assert translate(witness["cdna"][witness["frame"]:])[0] == mode["protein"]
    assert max(r["prefix"] for r in data["znf436"][1]["alt_observations"]) + 1 + max(
        r["suffix"] for r in data["znf436"][1]["alt_observations"]) == 113

    for source in data["cd109"]:
        observations = [defaultdict(set),defaultdict(set)]
        joint = set()
        for read in original_reads(source):
            coords = {g:q for q,g in read.get_aligned_pairs(matches_only=True)}
            key = (read.get_tag("RG") if read.has_tag("RG") else "",read.query_name)
            states = []
            for i,(position,ref,alt) in enumerate(((73818404,"G","T"),(73812223,"C","T"))):
                q = coords.get(position)
                value = None
                if q is not None and read.query_qualities is not None and read.query_qualities[q]>=20:
                    base = read.query_sequence[q]
                    value = "ref" if base==ref else "alt" if base==alt else "other"
                    observations[i][key].add(value)
                states.append(value)
            if states == ["alt","alt"]:
                joint.add(key)
        for observed,label in zip(observations,("focal","linked")):
            counts = dict(Counter(next(iter(s)) if len(s)==1 else "conflict" for s in observed.values()))
            assert counts == source[label]
        assert len(joint) == source["same_alignment_joint"].get("alt/alt",0)
        assert len(set(observations[0]) & set(observations[1])) == source["template_keys_observing_both"]
    return cd109_protein(data)


def cd109_protein(data):
    """A 222-nt actual RNA window independently supports the two-edit protein."""
    model = data["models"]["ENST00000287097"]
    lo,hi = [transcript_offset(model["exons"],"+",p) for p in (73812224,73818405)]
    start = lo - (lo-model["cds_start"])%3 - 6
    end = hi + 3 - (hi-model["cds_start"])%3 + 6
    sequence = list(model["cdna"])
    sequence[hi] = "T"
    single = "".join(sequence[start:end])
    sequence[lo] = "T"
    compound = "".join(sequence[start:end])
    positions = [p for a,b in model["exons"] for p in range(a-1,b)]
    matches = []
    for read in original_reads(data["cd109"][0]):
        coords = {g:q for q,g in read.get_aligned_pairs(matches_only=True)}
        a,b = coords.get(positions[start]),coords.get(positions[end-1])
        if a is None or b is None or read.query_sequence[a:b+1] != compound:
            continue
        if min(read.query_qualities[a:b+1]) >= 20:
            matches.append(dict(read_id=read.query_name,sequence=read.query_sequence[a:b+1],
                                source_query_start=a,minimum_quality=min(read.query_qualities[a:b+1])))
    assert len(matches) == 3 and len(compound) == 222
    return dict(reference=translate(model["cdna"][start:end])[0],single=translate(single)[0],
                compound=translate(compound)[0],witnesses=matches,transcript=model["transcript_version"],
                mutation_offsets=[(lo-start)//3,(hi-start)//3],protein_start=(start-model["cds_start"])//3)


def protein_rows(title, subtitle, rows, offsets, note, footer=""):
    figure, ax = _canvas(title,subtitle,5.4)
    differences = _protein_disagreements(dict(amino_acids=p, mutation_start=-start, ends_with_stop_codon=False)
                                        for _,p,_,start in rows)
    for y,(label,protein,color,start) in enumerate(rows[::-1],1):
        for i,aa in enumerate(protein,start):
            ax.text(i+.5,y,aa,ha="center",va="center",fontsize=11 if len(protein)>55 else 14,
                    fontfamily="monospace",color=color)
            if i in differences:
                ax.plot([i+.12,i+.88],[y-.22,y-.22],color="#AA3377",lw=2.5)
        ax.text(-.035,y,label,transform=ax.get_yaxis_transform(),ha="right",va="center",fontsize=11)
    for offset in offsets:
        ax.axvspan(offset,offset+1,color=ORANGE,alpha=.12,zorder=-1)
    ax.set(xlim=(min(r[3] for r in rows)-1,max(r[3]+len(r[1]) for r in rows)+1),
           ylim=(.4,len(rows)+.7),xlabel="Aligned protein offset (amino acids)")
    _side_note(ax,note + "\n\nMagenta: tracks differ")
    figure.text(.19,.06,footer,fontsize=10,color=GRAY)
    return figure


def context_panels(data, cd109):
    long,short = data["znf436"]
    rows=[]
    for label,source,color in (("ONT / assembly on",long,BLUE),("Illumina / on",short,GREEN),("Illumina / off",short,GRAY)):
        mode=source["modes"][1 if label.endswith("off") else 0]
        rows.append((label,mode["protein"],color,-mode["mutation"][0]))
    yield "ZNF436-length", "protein", protein_rows("ZNF436   T1 matched timepoint", "Read length limits recovered protein context",
        rows,[0],"ONT: 49 aa\nIllumina on: 30 aa\nIllumina off: 29 aa\n\n3 alternate cell/UMI\ngroups per library",
        "Same patient/timepoint, different libraries; not a controlled platform experiment. All translations independently validated.")
    figure,ax=_canvas("ZNF436   T1 matched timepoint","Observed alternate-read spans",6.8)
    rows=[("ONT",a,BLUE) for a in long["alt_observations"]]+[("Illumina",a,GREEN) for a in short["alt_observations"]]
    for y,(label,row,color) in enumerate(rows[::-1],1):
        ax.plot([max(-80,-row["prefix"]),min(80,row["suffix"]+1)],[y,y],color=color,lw=4,solid_capstyle="butt")
        ax.text(-.03,y,label,transform=ax.get_yaxis_transform(),ha="right",va="center",fontsize=10)
    ax.axvline(0,color=ORANGE,ls="--",lw=1)
    ax.set(xlim=(-82,82),ylim=(.3,len(rows)+.7),xlabel="cDNA offset from alternate base (nt)")
    _side_note(ax,"3 ONT reads\n8 Illumina reads\n\nIllumina union: 113 nt\n49 aa require 147 nt\n\nONT spans are clipped\nto this viewing window")
    figure.text(.19,.06,"The short-read union is an optimistic bound, not an asserted phased assembly. Extra PCR copies do not add span.",fontsize=10,color=GRAY)
    yield "ZNF436-length","read-spans",figure

    figure,ax=_canvas("CD109   T2","Long RNA molecules directly connect both alleles",5.8)
    categories=[("ONT: alternate / alternate",21,BLUE), ("Illumina: focal alternate",22,GREEN),
                ("Illumina: linked alternate",19,GREEN), ("Illumina: observed at both",0,GRAY)]
    for y,(label,count,color) in enumerate(categories[::-1],1):
        ax.barh(y,count,color=color,height=.38)
        ax.text(-.03,y,label,transform=ax.get_yaxis_transform(),ha="right",va="center",fontsize=10)
        ax.text(count+.35,y,str(count),va="center",fontsize=12)
    ax.set(xlim=(0,27),ylim=(.3,4.8),xlabel="Distinct read/template IDs passing Q20 at the queried bases")
    _side_note(ax,"207 nt apart in cDNA\n6,181 bp in genome\n\n21 linked ONT UMIs\n0 linked Illumina UMIs\n\nShort reads: 90 nt")
    figure.text(.19,.06,"Separate short-read support at each site does not establish their phase. No single 25-mer contains both changed codons.",fontsize=10,color=GRAY)
    yield "CD109-phase","phase-support",figure
    yield "CD109-phase","protein",protein_rows("CD109   T2", "Observed two-edit protein context versus a single-edit prediction",
        [("Reference",cd109["reference"],GRAY,0),("Focal edit only",cd109["single"],GRAY,0),
         ("RNA: both edits",cd109["compound"],BLUE,0)],cd109["mutation_offsets"],
        "3 exact 222-nt reads\nQ20 at every base\n\n"+cd109["transcript"]+"\n(CD109-001)",
        "Local conditional translation, not a full-length Isovar protein. The second allele's germline/somatic origin is not established.")

    m=data["map2"]["proteins"][0]
    yield "MAP2-haplotypes","protein-hypotheses",protein_rows("MAP2   distinct alleles", "Reference-based predictions for explicit haplotypes",
        [("Listed 22-nt deletion",m["listed_22del"]["window"],GRAY,0),
         ("28-nt deletion alone",m["deletion_only"]["window"],GRAY,0),
         ("Observed compound",m["combined"]["window"],BLUE,0)],list(range(5,10)),
        "Predicted full lengths\n22-nt: 906 aa\n28-nt: 904 aa\nCompound: 904 aa\n\n"+m["transcript"]+"\n(MAP2-001)",
        "These are edited-reference predictions, not whole RNA protein reconstructions. Compound = 28-nt deletion + two substitutions.")
    figure,ax=_canvas("MAP2   anchored RNA and DNA evidence","Exact compound-haplotype support",6)
    selected=[("T1 DNA","d9200d085cb7774d"),("T2 DNA","707928e142bdb21a"),
              ("T1 bulk Illumina RNA","261f8be4ce5aa8c8"),("T2 bulk Illumina RNA","95a1e72e7ec10990"),
              ("T2 deduplicated ONT RNA","2bc1fc291308debb")]
    sources={s["source_id"]:s for s in data["map2"]["sources"]}
    for y,(label,sid) in enumerate(selected[::-1],1):
        count=sources[sid]["counts"]["q20"].get("combined",0)
        ax.barh(y,count,color=BLUE if "RNA" in label else GRAY,height=.4)
        ax.text(-.03,y,label,transform=ax.get_yaxis_transform(),ha="right",va="center",fontsize=10)
        ax.text(count+.3,y,str(count),va="center",fontsize=12)
    ax.set(xlim=(0,29),ylim=(.3,5.8),xlabel="Distinct templates with exact anchored sequence and Q20 at every retained base")
    _side_note(ax,"No exact isolated\n22-nt allele observed\nin the 35-product audit\n\nProcessed copies are\nnot summed together")
    figure.text(.19,.06,"Both genomic anchors must align. A large annotated splice skip is not counted as a deletion observation.",fontsize=10,color=GRAY)
    yield "MAP2-haplotypes","haplotype-support",figure


def nr2f2_panels(data):
    """Keep the DNA/RNA discrepancy and molecular-independence caveat visible."""
    figure, ax = _canvas("NR2F2 | RNA-supported, DNA-unconfirmed",
                         "Matched CeGaT T0 samples | chr15:96875577-96875579 | GRCh37", 5.0)
    columns = [0, 3.0, 5.0, 7.0]
    for x, label in zip(columns, ["Sample", "GTG retained", "GTG deleted", "Other"]):
        ax.text(x, 3.2, label, fontsize=12, weight="bold", color=GRAY)
    for y, key in zip([2.4, 1.65, 0.9], ["normal", "tumor", "rna"]):
        source = data["sources"][key]
        counts = source["audits"]["q20"]["counts"]["deletion"]
        ax.text(0, y, source["label"], fontsize=15, va="center")
        for x, label, color in zip(columns[1:], ["reference", "deletion", "other"], [BLUE, "#b02a78", GRAY]):
            ax.text(x + .5, y, str(counts.get(label, 0)), fontsize=22, color=color, va="center", ha="center")
        ax.axhline(y-.32, color="#e7e7e7", linewidth=.7)
    ax.set(xlim=(-.1, 8.3), ylim=(0, 3.7))
    ax.axis("off")
    _side_note(ax, "5 RNA templates\n1 endpoint family\n\nNot five proven\nindependent molecules.\n\nNo germline/somatic\norigin assigned.")
    figure.text(.08, .07, "Exact sequence between eight-base flanks; Q20 across the window; MAPQ >=20 (255 retained).\n"
                "Primary/QC-passing; duplicate-flagged records excluded; mates counted once. Other alleles remain visible.",
                fontsize=10, color=GRAY)
    yield "NR2F2-evidence", "dna-rna", figure

    figure, ax = _canvas("NR2F2 | Directly observed local haplotypes",
                         "Focal G>T and nearby deletion assessed together on one aligned segment", 5.4)
    labels = ["G + retained GTG", "T + retained GTG", "G + deleted GTG", "T + deleted GTG", "Other"]
    keys = ["reference/reference", "alternate/reference", "reference/deletion", "alternate/deletion"]
    for x, text in zip([0, 4.0, 5.7, 7.4], ["Same-segment haplotype", "Blood DNA", "Tumor DNA", "Tumor RNA"]):
        ax.text(x, 4.1, text, fontsize=11, weight="bold", color=GRAY, ha="left" if x == 0 else "center")
    for y, label in zip([3.35, 2.65, 1.95, 1.25, .55], labels):
        ax.text(0, y, label, fontsize=13, va="center")
        for x, sample in zip([4.0, 5.7, 7.4], ["normal", "tumor", "rna"]):
            counts = data["sources"][sample]["audits"]["q20"]["counts"]["joint"]
            i = labels.index(label)
            count = counts.get(keys[i], 0) if i < 4 else sum(v for k, v in counts.items() if k not in keys)
            ax.text(x, y, str(count), fontsize=19, va="center", ha="center",
                    color="#b02a78" if i == 3 else BLUE if i < 2 else GRAY)
        ax.axhline(y-.3, color="#e7e7e7", linewidth=.7)
    ax.set(xlim=(-.1, 8.3), ylim=(0, 4.6))
    ax.axis("off")
    _side_note(ax, "4 RNA templates\nphase T with deletion.\n\nSame endpoint family;\nPCR or alignment\nartifact not excluded.\n\nDNA supports T\nwithout the deletion.")
    figure.text(.08, .065, "One template counted once; only records passing both focal-base and full deletion-window Q20 checks.\n"
                "This is direct phase in observed RNA, not validation of a distinct biological deletion.", fontsize=10, color=GRAY)
    yield "NR2F2-evidence", "direct-haplotypes", figure


def generate(output_dir):
    data=load_context()
    cd109=verify_context(data)
    output=osteosarc_assembly_figures.generate(output_dir)
    _,rc_context,_=_plot_imports()
    from matplotlib.backends.backend_pdf import PdfPages
    from pypdf import PdfReader, PdfWriter

    custom=defaultdict(list)
    nr2f2_data = nr2f2.load()
    for source in nr2f2_data["sources"].values():
        assert nr2f2.recount(source, nr2f2_data["reference"]) == source["audits"]["q20"]
    for group,name,figure in context_panels(data,cd109):
        custom[group].append((name,figure))
    for group,name,figure in nr2f2_panels(nr2f2_data):
        custom[group].append((name,figure))
    with rc_context({"svg.fonttype":"none","pdf.fonttype":42}):
        for group,panels in custom.items():
            directory=output/group
            directory.mkdir()
            with PdfPages(directory/"all-figures.pdf") as pdf:
                for name,figure in panels:
                    figure.savefig(directory/(name+".png"),dpi=600,facecolor="white")
                    figure.savefig(directory/(name+".svg"),facecolor="white",metadata={"Date":None})
                    pdf.savefig(figure,facecolor="white")
                    figure.clear()
            evidence={"ZNF436-length":data["znf436"],"CD109-phase":dict(sources=data["cd109"],protein=cd109),
                      "MAP2-haplotypes":data["map2"], "NR2F2-evidence":nr2f2_data}[group]
            (directory/"evidence.json").write_text(json.dumps(evidence,indent=2)+"\n")
    fusion_entries = [(FUSIONS, e) for e in json.loads((FUSIONS/"manifest.json").read_text())]
    coding = FUSIONS.parent / "coding-corpus"
    fusion_entries += [(coding, e) for e in json.loads((coding/"manifest.json").read_text())]
    for corpus, entry in fusion_entries:
        if "input" not in entry:
            continue
        raw=(corpus/entry["input"]).read_bytes()
        assert sha256(raw).hexdigest()==entry["sha256"]
        supplied=json.loads(gzip.decompress(raw))
        fusion,refs,reads=fusion_from_dict(supplied)
        result=reconstruct_fusion(fusion,refs,reads)
        directory=save_fusion_figures(result,output/"fusion-rna",refs,supplied.get("reference_names"))
        (directory/"input.json").write_text(json.dumps(supplied,indent=2)+"\n")
    combined=PdfWriter()
    index=[]
    priority = ["DIAPH1-", "SLC25A12-", "TECPR1-", "ZNF436-length", "CD109-phase", "MAP2-haplotypes"]
    def order(path):
        name=str(path.relative_to(output))
        return (next((i for i,prefix in enumerate(priority) if name.startswith(prefix)),len(priority)),name)
    for path in sorted(output.rglob("all-figures.pdf"),key=order):
        pages=len(PdfReader(path).pages)
        index.append(dict(path=str(path.relative_to(output)),start_page=len(combined.pages)+1,pages=pages))
        combined.append(path,outline_item=str(path.parent.relative_to(output)))
    combined.add_metadata({"/Title":"Isovar: RNA reconstruction, haplotypes and fusion evidence"})
    combined.write(output/"isovar-all-figures.pdf")
    (output/"figure-index.json").write_text(json.dumps(index,indent=2)+"\n")
    with (output/"README.md").open("a") as handle:
        handle.write("\n\n## Extended evidence gallery\n\nCombined vector PDF: `isovar-all-figures.pdf`. "
            "Page index and bookmarks preserve individual examples; all panels also have separate 600-dpi PNG and SVG files. "
            "ZNF436 isolates a read-span limitation; CD109 demonstrates direct phase; MAP2 separates distinct haplotypes. "
            "ATP5MG-KMT2A (Sid) and BCR-ABL1 (external K562 control) have RNA-backed coding hypotheses; "
            "compatible noncoding/alternative-CDS annotations remain explicit. The other fusion windows "
            "remain unresolved. No uniquely expressed fusion protein, long-read-only fusion detection, "
            "or clinical suitability is claimed.\n\n")
        for item in index:
            handle.write("- Page %d: [%s](%s), %d pages\n" % (item["start_page"],item["path"],item["path"],item["pages"]))
    return output


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",default="figures/osteosarc")
    args=parser.parse_args()
    logging.disable(logging.INFO)
    print(generate(args.output_dir))


if __name__ == "__main__":
    main()
