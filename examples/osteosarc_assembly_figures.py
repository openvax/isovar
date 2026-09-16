"""Reproduce white-background figures from pinned, original osteosarc RNA.

Run from a repository checkout: python -m examples.osteosarc_assembly_figures
Only rendering adds files; input reads and scientific defaults are unchanged.
"""

import argparse
import json
import logging
from pathlib import Path
import tempfile

import pysam
from varcode import Variant

from isovar import ProteinSequenceCreator, ReadCollector
from isovar.default_parameters import PLOT_DPI
from isovar.variant_helpers import base0_interval_for_variant
from isovar.visualization import collect_visualization_data, save_variant_figures, timestamped_run_directory
from tests.data.osteosarc.expansion.inventory import digest
from tests.data.osteosarc.expansion.references import apply_variant, load_reference, reference_genome
from tests.data.osteosarc.expansion.runner import protein_check


CORPUS = Path(__file__).resolve().parents[1] / "tests/data/osteosarc/expansion/corpus"
CASE_IDS = (
    "09-DIAPH1-chr5-141526090-1b66c15da594a3ef",
    "34-SLC25A12-chr2-171813470-fbea0bf403520431",
    "36-TECPR1-chr7-98241127-8fda918141cedf55",
)
ADDITIONAL_CASE_IDS = (
    "10-DYNC1H1-chr14-101980529-53f498a544883d51",
    "11-DYNC1H1-chr14-102030200-1b66c15da594a3ef",
    "17-H1_2-chr6-26055824-1120a096937e29e1",
    "21-MAP2-chr2-209694768-2bc1fc291308debb",
    "25-NAV2-chr11-20080142-8fda918141cedf55",
    "28-NTF3-chr12-5494381-53f498a544883d51",
    "48-NR2F2-chr15-96332299-fbeb5d8a415e4d18",
)
PAIR_CORPUS = CORPUS.parents[1] / "figure_comparisons/corpus"
PAIR_CASE_IDS = ("PIP5K1A-T1-ONT", "PIP5K1A-T1-Illumina")


def case_caption(case, data):
    on, off = [m["protein"] for m in data["modes"]]
    gene = case["variant"]["gene"]
    if on is None:
        summary = ("No RNA-derived protein under the default coverage floor of two read objects. "
                   "%d alternate-supporting object(s) are available. Varcode supplies a prediction, "
                   "not evidence of mutant expression." % data["counts"]["alt"]["observations"])
    else:
        summary = ("Assembly on: %d aa and %d mutation-overlapping 25-mers; assembly off: %d aa and %d 25-mers. "
                   "Both outputs pass independent RNA translation, frame and interval checks. " % (
                       len(on["amino_acids"]), on["peptide_windows"],
                       len(off["amino_acids"]) if off else 0, off["peptide_windows"] if off else 0))
        if case["case_id"] in CASE_IDS:
            if not (off and len(on["amino_acids"]) > len(off["amino_acids"])
                    and off["amino_acids"] in on["amino_acids"] and not on["witness"]["spanning_observations"]):
                raise ValueError("Example is no longer simple context recovery: " + case["case_id"])
            summary += ("Both also match the independent reference-plus-edit expectation. The shorter sequence "
                        "is not mistranslated. No single read object spans the full reconstructed cDNA. ")
    if gene == "NAV2":
        summary += ("The current splice-aware output retains %d contributing transcript models, not a unique isoform. "
                    "Other Varcode predictions remain visible for comparison, not as isoform assignments. "
                    "The earlier seven-isoform V/L comparison predates splice-path restriction and is not the "
                    "current result. Black transcript models contribute to this RNA protein; gray models do not. "
                    "The 47-aa RNA sequence ends at W, before the reference predictions diverge into LR or VN; "
                    "those extra residues are predicted, not reconstructed. No junction remains inside this "
                    "retained cDNA window. Path compatibility may use alignment beyond the plotted window. "
                    % len(on["transcript_ids"]))
    elif gene == "NTF3":
        summary += ("For the nominated chr12:5494381 A>G call, Varcode predicts K56R (K69R on the other model). "
                    "RNA instead encodes serine: the adjacent change gives the source-linked compound AG>GT allele. "
                    "This is a translated RNA-versus-single-edit difference, not proof of a germline or somatic origin "
                    "for the second base. Both assembly modes recover it; assembly is not required here. ")
    elif gene == "H1-2":
        summary += ("H1-2 is named HIST1H1C in the pinned Ensembl 87 annotation. This is a 15-nt in-frame "
                    "deletion, not a structural variant. Varcode's repeat-normalized deletion label differs from "
                    "the source label; the full mutant sequence is independently checked. ")
    elif gene == "NR2F2":
        summary += ("This is the native GRCh37 CeGaT RNA product, not the GRCh38 locus with no alternate reads. "
                    "Three original templates carry the focal allele and a downstream 3-nt CIGAR deletion "
                    "at chr15:96875577-96875579. Both assembly modes retain it, whereas the single-edit "
                    "Varcode prediction does not. The in-frame deletion changes the downstream protein context; "
                    "its biological origin is not established. Magenta marks position-wise differences, "
                    "including the shifted sequence after that deletion. ")
    elif gene == "MAP2":
        summary += ("The nominal 22-nt call must not be conflated with the separately observed compound haplotype: "
                    "a 28-nt deletion plus two substitutions. Its local RNA sequence can align as 22D/2M/6D with "
                    "a mismatch. This single nominal alternate object is not independent confirmation of an "
                    "isolated 22-nt allele. See the separate MAP2 haplotype figures and anchored-sequence audit. ")
    if case["case_id"] in PAIR_CASE_IDS:
        summary += ("This is the complete acquired locus, not a selected read subset: " + case["case_id"] + ". "
                    "ONT and Illumina products are both catalogued T1, but are different processed libraries. "
                    "The ONT sequence passes independent reference-plus-edit validation; the Illumina product "
                    "has only one alternate read object. This demonstrates insufficient support in that product, "
                    "not that short reads intrinsically cannot reconstruct the sequence or that read length alone "
                    "explains the difference. ONT also reconstructs it without overlap assembly. ")
        return summary + "No cell/molecule-level independence or matched malignant-cell composition is established."
    return summary + ("Selected primary-only original RNA fixtures are not full-sample abundance estimates. "
                      "Varcode tracks assume only the nominated edit on each reference transcript.")


def nonfocal_witness_indels(protein):
    """Record retained alignment indels separately from the nominated allele."""
    indels = set()
    for translation in protein.translations:
        sequence = translation.untrimmed_variant_sequence
        focal = base0_interval_for_variant(translation.reference_context.variant)
        for read in sequence.reads:
            q0 = max(0, len(read.prefix) - len(sequence.prefix))
            q1 = min(len(read.sequence), len(read.prefix) + len(read.allele) + len(sequence.suffix))
            for left, right in zip(read.reference_blocks, read.reference_blocks[1:]):
                gap = (left[3], right[2])
                if (q0 < left[1] <= right[0] < q1 and gap not in read.splice_junctions
                        and (right[0] != left[1] or gap[0] != gap[1]) and gap != focal):
                    indels.add((gap[0], gap[1], right[0] - left[1]))
    return sorted(indels)


def generate(output_dir, case_ids=CASE_IDS + ADDITIONAL_CASE_IDS + PAIR_CASE_IDS):
    """Validate both modes and save reproducible original-read figures."""
    manifest = json.loads((CORPUS / "manifest.json").read_text())
    cases = {c["case_id"]: c for c in manifest["cases"]}
    if set(case_ids).intersection(PAIR_CASE_IDS):
        cases.update({c["case_id"]: c for c in json.loads((PAIR_CORPUS / "manifest.json").read_text())["cases"]})
    output = timestamped_run_directory(output_dir)
    summaries = []
    with tempfile.TemporaryDirectory(prefix="isovar-figure-reference-") as cache:
        reference_data = {}
        for case_id in case_ids:
            case = cases[case_id]
            reference_name = case.get("reference", "GRCh38")
            reference_dir = CORPUS / "references" / reference_name
            if reference_name not in reference_data:
                reference_manifest, models = load_reference(reference_dir)
                reference_data[reference_name] = (reference_manifest, models,
                    reference_genome(reference_dir, Path(cache) / reference_name))
            reference_manifest, models, genome = reference_data[reference_name]
            case_corpus = PAIR_CORPUS if case_id in PAIR_CASE_IDS else CORPUS
            r = case["variant"]
            bam_path = case_corpus / case["primary_bam"]
            for filename in (case["primary_bam"], case["primary_bam"] + ".bai"):
                if digest(case_corpus / filename) != case["files"][filename]:
                    raise ValueError("Fixture checksum mismatch: " + filename)
            variant = Variant(r["chrom"].removeprefix("chr"), r["pos"], r["ref"], r["alt"], ensembl=genome)
            expected = {tid: apply_variant(r, models[tid])
                        for tid in reference_manifest["variant_transcripts"][r["variant_id"]]}
            with pysam.AlignmentFile(bam_path) as bam:
                evidence = ReadCollector(use_secondary_alignments=False).read_evidence_for_variant(variant, bam)
            data = collect_visualization_data(variant, evidence, compare_assembly=True,
                                              transcript_id_whitelist=set(expected))
            data["genes"] = [r["gene"]]
            for prediction in data["reference_predictions"]:
                p = prediction["protein"]
                if p:
                    full = expected[prediction["transcript_id"]]["mutant_protein"]
                    offset = p["protein_start_1based"] - 1
                    if p["amino_acids"] != full[offset:offset + len(p["amino_acids"])]:
                        raise ValueError("Varcode baseline disagrees with independent reference edit: " + case_id)
            validation = []
            retained_indels = set()
            for assembly in (True, False):
                creator = ProteinSequenceCreator(variant_sequence_assembly=assembly)
                proteins = creator.sorted_protein_sequences_for_variant(
                    variant, evidence, transcript_id_whitelist=set(expected))
                if not proteins:
                    validation.append(dict(assembly=assembly, checks=[], status="no_protein"))
                    continue
                retained_indels.update(nonfocal_witness_indels(proteins[0]))
                checked = protein_check(proteins[0], expected, creator.protein_sequence_length,
                                        creator.protein_context_peptide_length)
                if checked["validation_status"] != "ok":
                    raise ValueError("Independent validation failed: " + case_id)
                # Record independent checks but omit the public source read names.
                validation.append(dict(assembly=assembly, checks=checked["checks"], status="ok"))
            on, off = [m["protein"] for m in data["modes"]]
            caption = case_caption(case, data)
            data["provenance"] = dict(
                case_id=case_id, source_id=case["source_id"], source_variant_url=r["source_url"],
                source_bam_url=case["source_url"],
                fixture=case["primary_bam"], fixture_sha256=digest(bam_path),
                corpus_manifest_sha256=digest(case_corpus / "manifest.json"),
                reference_manifest_sha256=digest(reference_dir / "manifest.json"),
                reference_files=reference_manifest["files"],
                original_allele=r, selection=("Complete acquired locus; primary-only derivative" if case_id in PAIR_CASE_IDS
                                             else "Pinned selected original reads; primary-only derivative"),
                excluded_flags=[256, 1024, 2048], merge_overlapping_fragments=True,
                nonfocal_witness_indels=len(retained_indels), nonfocal_alignment_indels=sorted(retained_indels),
                sample_label=case_id if case_id in PAIR_CASE_IDS else None,
                independent_validation=validation, caption=caption)
            directory = save_variant_figures(data, output / case_id if case_id in PAIR_CASE_IDS else output)
            (directory / "README.md").write_text(
                "# " + r["gene"] + " assembly comparison\n\n" + caption + "\n\n"
                + "\n\n".join("## " + label + "\n\n![" + label + "](" + name + ".png)\n\n"
                              + "[PNG](" + name + ".png) · [SVG](" + name + ".svg)"
                              for name, label in (("protein", "Protein context"), ("reads", "Read overlaps"),
                                                  ("coverage", "RNA coverage"), ("transcripts", "Transcripts")))
                + "\n\n[All panels PDF](all-figures.pdf) · [Overview](overview.svg) · [Figure notes](caption.md) · [Evidence](evidence.json)\n\n"
                + "Source variant: " + r["source_url"] + "\n\n"
                + "Source RNA product: `" + case["source_id"] + "`. Fixture: `" + case["primary_bam"] + "`.\n\n"
                + "See evidence.json for exact settings, transcript models, checksums and independent validation.\n")
            summaries.append(dict(case_id=case_id, directory=str(directory.relative_to(output)), caption=caption,
                                  on_length=len(on["amino_acids"]) if on else 0,
                                  off_length=len(off["amino_acids"]) if off else 0,
                                  on_windows=on["peptide_windows"] if on else 0,
                                  off_windows=off["peptide_windows"] if off else 0))
    (output / "manifest.json").write_text(json.dumps(dict(examples=summaries), indent=2) + "\n")
    (output / "README.md").write_text(
        "# Osteosarc assembly figures\n\n"
        "Reproduce: `python -m examples.osteosarc_assembly_figures --output-dir figures/osteosarc`\n\n"
        "White-background SVG (editable vector text) and " + str(PLOT_DPI) + "-dpi PNG. "
        "Each variant has separate protein, read-overlap, coverage and transcript figures, plus a 300-dpi overview preview, "
        "a multipage all-figures.pdf, captions and evidence.json. Source reads are public CC0 osteosarc.com data; see "
        "tests/data/osteosarc/README.md for source and selection details.\n\n"
        "Assembly recovers missing context in the first three examples; later cases show agreement, "
        "insufficient support, transcript ambiguity and RNA differences from a single-edit prediction. "
        "This does not establish a unique "
        "isoform, full-length protein, peptide presentation, or clinical suitability. The comparison "
        "changes only overlap assembly; mate merging stays enabled in both modes.\n\n"
        + "\n\n".join("## " + s["case_id"] + "\n\n" + s["caption"]
                      + "\n\n![Protein comparison](" + s["directory"] + "/protein.png)"
                      + "\n\n![Read overlaps](" + s["directory"] + "/reads.png)"
                      + "\n\n[All panels and captions](" + s["directory"] + "/README.md) · [Evidence]("
                      + s["directory"] + "/evidence.json)"
                    for s in summaries) + "\n")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="figures/osteosarc")
    args = parser.parse_args()
    logging.disable(logging.INFO)
    print(generate(args.output_dir))


if __name__ == "__main__":
    main()
