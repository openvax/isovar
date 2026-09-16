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


def check_no_nonfocal_witness_indels(protein):
    """Keep these demonstrations outside the known upstream-indel bug #265."""
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
                    raise ValueError("Non-focal retained indel: figure needs a separate frame audit (#265)")


def generate(output_dir):
    """Validate both modes and save a reproducible three-example figure batch."""
    manifest = json.loads((CORPUS / "manifest.json").read_text())
    cases = {c["case_id"]: c for c in manifest["cases"]}
    reference_dir = CORPUS / "references/GRCh38"
    reference_manifest, models = load_reference(reference_dir)
    output = timestamped_run_directory(output_dir)
    summaries = []
    with tempfile.TemporaryDirectory(prefix="isovar-figure-reference-") as cache:
        genome = reference_genome(reference_dir, cache)
        for case_id in CASE_IDS:
            case = cases[case_id]
            r = case["variant"]
            bam_path = CORPUS / case["primary_bam"]
            for filename in (case["primary_bam"], case["primary_bam"] + ".bai"):
                if digest(CORPUS / filename) != case["files"][filename]:
                    raise ValueError("Fixture checksum mismatch: " + filename)
            variant = Variant(r["chrom"].removeprefix("chr"), r["pos"], r["ref"], r["alt"], ensembl=genome)
            expected = {tid: apply_variant(r, models[tid])
                        for tid in reference_manifest["variant_transcripts"][r["variant_id"]]}
            with pysam.AlignmentFile(bam_path) as bam:
                evidence = ReadCollector(use_secondary_alignments=False).read_evidence_for_variant(variant, bam)
            data = collect_visualization_data(variant, evidence, compare_assembly=True,
                                              transcript_id_whitelist=set(expected))
            validation = []
            for assembly in (True, False):
                creator = ProteinSequenceCreator(variant_sequence_assembly=assembly)
                proteins = creator.sorted_protein_sequences_for_variant(
                    variant, evidence, transcript_id_whitelist=set(expected))
                if not proteins:
                    raise ValueError("Selected example no longer has a protein: " + case_id)
                check_no_nonfocal_witness_indels(proteins[0])
                checked = protein_check(proteins[0], expected, creator.protein_sequence_length,
                                        creator.protein_context_peptide_length)
                if checked["validation_status"] != "ok" or not checked["matches_expected"]:
                    raise ValueError("Independent validation failed: " + case_id)
                # Record independent checks but omit the public source read names.
                validation.append(dict(assembly=assembly, checks=checked["checks"]))
            on, off = [m["protein"] for m in data["modes"]]
            if not (len(on["amino_acids"]) > len(off["amino_acids"])
                    and off["amino_acids"] in on["amino_acids"]):
                raise ValueError("Example is no longer simple context recovery: " + case_id)
            caption = (
                "Assembly recovers %d aa rather than %d aa, providing %d rather than %d "
                "mutation-overlapping 25-mers. Both sequences pass independent translation "
                "and transcript-reference checks; the shorter output is not mistranslated. "
                "No single read object spans the full reconstructed cDNA. "
                "These are selected primary-only original RNA fixtures, not full-sample abundance estimates."
            ) % (len(on["amino_acids"]), len(off["amino_acids"]), on["peptide_windows"], off["peptide_windows"])
            if on["witness"]["spanning_observations"]:
                raise ValueError("A single observation now spans the assembled witness: " + case_id)
            data["provenance"] = dict(
                case_id=case_id, source_id=case["source_id"], source_variant_url=r["source_url"],
                source_bam_url=case["source_url"],
                fixture=case["primary_bam"], fixture_sha256=digest(bam_path),
                corpus_manifest_sha256=digest(CORPUS / "manifest.json"),
                reference_manifest_sha256=digest(reference_dir / "manifest.json"),
                reference_files=reference_manifest["files"],
                original_allele=r, selection="Pinned selected original reads; primary-only derivative",
                excluded_flags=[256, 1024, 2048], merge_overlapping_fragments=True,
                nonfocal_witness_indels=0,
                independent_validation=validation, caption=caption)
            directory = save_variant_figures(data, output)
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
            summaries.append(dict(case_id=case_id, directory=directory.name, caption=caption,
                                  on_length=len(on["amino_acids"]), off_length=len(off["amino_acids"]),
                                  on_windows=on["peptide_windows"], off_windows=off["peptide_windows"]))
    (output / "manifest.json").write_text(json.dumps(dict(examples=summaries), indent=2) + "\n")
    (output / "README.md").write_text(
        "# Osteosarc assembly figures\n\n"
        "Reproduce: `python -m examples.osteosarc_assembly_figures --output-dir figures/osteosarc`\n\n"
        "White-background SVG (editable vector text) and " + str(PLOT_DPI) + "-dpi PNG. "
        "Each variant has separate protein, read-overlap, coverage and transcript figures, plus a 300-dpi overview preview, "
        "a multipage all-figures.pdf, captions and evidence.json. Source reads are public CC0 osteosarc.com data; see "
        "tests/data/osteosarc/README.md for source and selection details.\n\n"
        "Assembly recovers missing local context in these examples. This does not establish a unique "
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
