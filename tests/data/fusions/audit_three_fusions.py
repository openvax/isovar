"""Reconstruct both orientations and independently audit original ORF witnesses.

This analysis helper is deliberately kept with the test-data construction code.
Native query orientation and poly(A) tails are RNA polarity evidence, not proof
of translation; reverse-complement ORFs are retained and labelled separately.
"""
import argparse
from collections import Counter, defaultdict
import gzip
from hashlib import sha256
import json
from pathlib import Path
import tempfile

import pysam

from isovar.sid_data import sam_digest
from isovar.sv_rna import reconstruct_sv_rna, sv_rna_input_from_dict
from tests.data.fusions.build_three_fusions import segment
from tests.data.osteosarc.expansion.references import translate
from tests.osteosarc_protein_helpers import reverse_complement

DATA = Path(__file__).with_name("three-fusions")


def pinned_json(directory, entry):
    raw = (directory / entry["file"]).read_bytes()
    if sha256(raw).hexdigest() != entry["sha256"]:
        raise ValueError("Fixture checksum mismatch: " + entry["file"])
    return json.loads(gzip.decompress(raw))


def run_case(directory, event, source, orientation, scratch):
    manifest = json.loads((directory / "manifest.json").read_text())
    entry, = [e for e in manifest["fixtures"] if e["event"] == event and e["source"] == source]
    data = pinned_json(directory, entry)
    models = pinned_json(directory, manifest["references"])[event]
    header = pysam.AlignmentHeader.from_dict(data["header"])
    records = []
    for checksum, line in sorted(data["records"].items()):
        if sam_digest(line) != checksum:
            raise ValueError("Original SAM checksum mismatch")
        records.append(pysam.AlignedSegment.fromstring(line, header))
    unsorted, path = scratch / "unsorted.bam", scratch / "reads.bam"
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for read in records:
            out.write(read)
    pysam.sort("--no-PG", "-o", str(path), str(unsorted))
    pysam.index(str(path))
    definition = manifest["events"][event]
    donor, acceptor = definition["donor"], definition["acceptor"]
    if orientation == "reverse":
        donor, acceptor = (dict(b, strand="-" if b["strand"] == "+" else "+")
                           for b in (acceptor, donor))
    inputs = sv_rna_input_from_dict(dict(
        event_id=event, reference_name="GRCh38", sample_id=source.split("-")[-1],
        donor=donor, acceptor=acceptor, references=models,
        event_provenance=dict(source=manifest["event_source"], orientation=orientation,
                              fixture=entry["file"])))
    with pysam.AlignmentFile(str(path)) as bam:
        result = reconstruct_sv_rna(
            bam, source=manifest["sources"][source]["url"], assemble=False,
            min_orf_amino_acids=8, max_records=50000, max_paths=200,
            max_extension_segments=500, **inputs)
    return result, records


def audit_candidates(result, records):
    groups = defaultdict(list)
    for read in records:
        groups[segment(read)].append(read)
    candidates = {}
    for path in result["paths"]:
        for candidate in path["exploratory_orfs"]["candidates"]:
            if not candidate["ends_with_stop_codon"] or not any(
                    path["junctions"][j]["relation"] == "breakpoint_junction"
                    for j in candidate["crossed_junctions"]):
                continue
            start, end = candidate["query_interval"]
            nt, aa = path["sequence"][start:end], candidate["amino_acids"]
            if translate(nt, table=1) != (aa, True) or len(nt) != 3 * (len(aa) + 1):
                raise ValueError("Independent translation disagrees")
            item = candidates.setdefault((nt, aa), dict(
                nucleotides=nt, amino_acids=aa, annotated_start=candidate["annotated_start"],
                initiation_observed=candidate["initiation_observed"],
                translation_observed=candidate["translation_observed"], witnesses=set()))
            for witness in candidate["full_interval_support"]["witnesses"]:
                identity = result["observations"][witness["observation"]]["identity"]
                item["witnesses"].add(tuple(identity))
    rows = []
    for item in candidates.values():
        nt = item["nucleotides"]
        evidence = []
        for identity in sorted(item.pop("witnesses")):
            primary = [r for r in groups[identity] if not r.is_secondary and not r.is_supplementary]
            row = dict(identity=list(identity), orientation="unavailable", minimum_quality=None,
                       terminal_A_run=None, cell=None, umi=None)
            if len(primary) == 1 and primary[0].query_sequence:
                read = primary[0]
                sequence = read.get_forward_sequence()
                native, opposite = sequence.count(nt), sequence.count(reverse_complement(nt))
                row["orientation"] = ("native" if native == 1 and not opposite else
                                      "opposite" if opposite == 1 and not native else "ambiguous")
                row["terminal_A_run"] = len(sequence) - len(sequence.rstrip("A"))
                row["cell"] = read.get_tag("CB") if read.has_tag("CB") else None
                row["umi"] = read.get_tag("UB") if read.has_tag("UB") else None
                if row["orientation"] == "native" and read.query_qualities is not None:
                    offset = sequence.index(nt)
                    row["minimum_quality"] = min(read.get_forward_qualities()[offset:offset + len(nt)])
            evidence.append(row)
        native = [row for row in evidence if row["orientation"] == "native"]
        item.update(witnesses=evidence, orientation_counts=dict(Counter(r["orientation"] for r in evidence)),
                    native_fragments=len(native),
                    native_q10_fragments=sum(r["minimum_quality"] is not None and r["minimum_quality"] >= 10 for r in native),
                    native_missing_quality_fragments=sum(r["minimum_quality"] is None for r in native),
                    native_terminal_a10_fragments=sum(r["terminal_A_run"] >= 10 for r in native),
                    native_cb_ub_labels=len({(tuple(r["identity"][:1]), r["cell"], r["umi"])
                                             for r in native if r["cell"] and r["umi"]}))
        rows.append(item)
    return sorted(rows, key=lambda c: (c["amino_acids"], c["nucleotides"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DATA)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((args.data / "manifest.json").read_text())
    rows = []
    with tempfile.TemporaryDirectory() as scratch:
        for entry in manifest["fixtures"]:
            for orientation in ("forward", "reverse"):
                result, records = run_case(args.data, entry["event"], entry["source"], orientation, Path(scratch))
                row = dict(event=entry["event"], source=entry["source"], orientation=orientation,
                           limitations=result["limitations"], candidates=audit_candidates(result, records))
                rows.append(row)
                print(row["event"], row["source"], orientation,
                      [(c["amino_acids"], c["native_fragments"]) for c in row["candidates"] if c["native_fragments"]], flush=True)
    args.output.write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()
