"""Reconstruct the three strongest named Osteosarc fusion leads from original RNA.

python -m examples.osteosarc_sv_validation OUTPUT_DIRECTORY
Requires Varcode >=9.4.2 and installed Ensembl 95. Original RNA and Ensembl 87/115
RNA-model metadata are checked in; RNA acquisition/selection recipes live in
isovar.sid_data and tests/data/fusions/build_long_read.py. No local report,
sibling checkout, downloaded analysis script, or network is required.
"""
import argparse
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import tempfile

import pysam
from pyensembl import cached_release
import varcode
from varcode.sv_allele_parser import parse_symbolic_alt

from isovar import __version__
from isovar.sid_data import export_fixture
from isovar.sv_rna import reconstruct_sv_rna, sv_rna_input_from_dict
from isovar.sv_rna_comparison import compare_sv_rna_predictions

DATA = Path(__file__).resolve().parents[1] / "tests/data/fusions"


def predict(entry):
    row = entry["record"]
    variant = parse_symbolic_alt(row["chrom"], int(row["pos"]), row["ref"], row["alt"],
                                 info={"END": int(row["end"])}, genome=cached_release(95),
                                 convert_ucsc_contig_names=True)
    predictions, seen = [], set()
    for effect in variant.effects():
        candidates = [c.effect for c in getattr(effect, "candidates", ())] or [effect]
        for candidate in candidates:
            transcript = getattr(candidate, "transcript_id", None)
            partner = getattr(getattr(candidate, "partner_transcript", None), "id", None)
            sequence = getattr(candidate, "mutant_protein_sequence", None)
            key = type(candidate).__name__, transcript, partner, sequence
            if key in seen:
                continue
            seen.add(key)
            predictions.append(dict(prediction_id="model_%d" % len(predictions),
                effect_class=key[0], transcript_id=transcript, partner_transcript_id=partner,
                amino_acids=sequence, complete=False))
    return dict(event_id=entry["event_id"], reference_name="GRCh38", sample_id="T1",
                source=dict(annotator="varcode fast", version=varcode.__version__, ensembl_release=95,
                            vcf=entry), predictions=predictions)


def validate(entry, output, bam_path=None, rna_source=None, orientation="forward"):
    event = entry["event_id"]
    if event == "PARD3B--CDKN2B-AS1-CDKN2B" and bam_path is None:
        # Reuse the released original-read builder/audit rather than another
        # reconstruction of its source records, annotation or selection rules.
        from tests.data.fusions.audit_three_fusions import run_case
        with tempfile.TemporaryDirectory(prefix="isovar-sv-validation-") as temporary:
            result, _ = run_case(DATA / "three-fusions", event, "ONT-T1", orientation, Path(temporary))
        result["event_provenance"].update(dna_call=entry, alignment_scope="selected_fixture_records")
        predictions = predict(entry)
    else:
        reference_path = DATA / "corpus" / (event + "-T1.input.json.gz")
        data = json.loads(gzip.decompress(reference_path.read_bytes()))
        long_reads = json.loads((DATA / "long-read/manifest.json").read_text())
        selected = next((r for r in long_reads if r["event"] == event), None)
        fixture = ("fusions/long-read/" + selected["file"] if selected else
                   "fusions/corpus/" + reference_path.name + "#/original_records")
        source = rna_source or bam_path or (selected["url"] if selected else data["fusion"]["provenance"]["source"])
        fusion = dict(data["fusion"])
        if orientation == "reverse":
            fusion["donor"], fusion["acceptor"] = (
                dict(b, strand="-" if b["strand"] == "+" else "+")
                for b in (fusion["acceptor"], fusion["donor"]))
        inputs = sv_rna_input_from_dict(dict(event_id=event, reference_name="GRCh38", sample_id="T1",
            donor=fusion["donor"], acceptor=fusion["acceptor"], references=data["references"],
            event_provenance=dict(dna_call=entry, orientation=orientation, rna_fixture=fixture if bam_path is None else None,
                                 alignment_scope="selected_fixture_records" if bam_path is None else "supplied_alignment",
                                 reference_metadata_sha256=sha256(reference_path.read_bytes()).hexdigest())))
        predictions = predict(entry)
        with tempfile.TemporaryDirectory(prefix="isovar-sv-validation-") as temporary:
            sam, bam = Path(temporary) / "reads.sam", Path(temporary) / "reads.bam"
            if bam_path is None:
                export_fixture(fixture, sam)
                pysam.sort("-o", str(bam), str(sam))
                pysam.index(str(bam))
            with pysam.AlignmentFile(str(bam_path or bam)) as alignments:
                result = reconstruct_sv_rna(alignments, source=source, **inputs)
    comparison = compare_sv_rna_predictions(result, predictions)
    report = dict(isovar_version=__version__, result=result, prediction_comparison=comparison)
    (output / (event + "." + orientation + ".json")).write_text(json.dumps(report, indent=2) + "\n")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--bam", help="Use an indexed T1 RNA BAM (local path or public URL) instead of selected fixtures")
    parser.add_argument("--source", help="Original BAM source identity when --bam is a local acquisition")
    args = parser.parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=False)
    entries = json.loads((DATA / "validation-events.json").read_text())
    rows = []
    for entry in entries:
        for orientation in ("forward", "reverse"):
            report = validate(entry, args.output, args.bam, args.source, orientation)
            comparison = report["prediction_comparison"]
            classes = sorted({p["effect_class"] for p in comparison["predictions"]["predictions"]})
            hypotheses = comparison["hypotheses"] or [None]
            for h in hypotheses:
                support = h.get("orf", {}).get("rna_support", {}) if h else {}
                directions = support.get("fragment_query_orientations", {})
                rows.append(dict(event=entry["event_id"], orientation=orientation, dna_effect_classes="; ".join(classes),
                    rna_source=comparison["rna_source"],
                    rna_hypothesis_id=h["hypothesis_id"] if h else "",
                    rna_protein=h["amino_acids"] if h else "",
                    complete_candidate=h["complete_candidate"] if h else "",
                    full_interval_fragments=h["full_interval_fragments"] if h else "",
                    original_query_fragments=directions.get("original_query"),
                    reverse_complement_query_fragments=directions.get("reverse_complement"),
                    mixed_query_fragments=directions.get("mixed"),
                    uncertainty_flags="; ".join(h.get("orf", {}).get("uncertainty_flags", [])) if h else "",
                    comparison="; ".join(sorted({c["status"] for c in h["comparisons"]})) if h else comparison["status"],
                    limitations="; ".join(comparison["reconstruction_limitations"])))
            print(entry["event_id"], len(comparison["hypotheses"]), "RNA hypotheses", flush=True)
    with (args.output / "protein_comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
