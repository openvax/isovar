"""Independently recount representative pilot witnesses from original SAM text.

Usage: python examples/check_sid_sv_pilot.py AUDIT REPORT OUTPUT.json
Requires Biopython for translation independently of Isovar's genetic-code code.
This checks sequence/count evidence, not RNA polarity, SV origin or translation.
SAM orientation and hard clipping follow https://samtools.github.io/hts-specs/SAMv1.pdf.
Translation uses NCBI standard genetic code 1.
"""

import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import re

from Bio.Seq import Seq


def load(path):
    with gzip.open(path, "rt") as handle:
        return json.load(handle)


def check_candidate(candidate, details, references):
    dna = candidate["nucleotide_sequence"]
    protein = str(Seq(dna[:len(dna) // 3 * 3]).translate(table=1))
    assert protein.rstrip("*") == candidate["amino_acids"]
    assert protein.endswith("*") == candidate["ends_with_stop_codon"]
    fragments, segments, cells, missing, original_records = set(), set(), set(), set(), {}
    witness_count = reference_count = 0
    for occurrence in candidate["occurrences"]:
        for reference in occurrence["reference_comparisons"]:
            sequence = references["models"][reference["transcript_id"]]["sequence"][reference["transcript_offset"]:]
            protein = str(Seq(sequence[:len(sequence) // 3 * 3]).translate(table=1))
            assert protein.split("*")[0] == reference["amino_acids"]
            reference_count += 1
        for witness in occurrence["witnesses"]:
            observation = details["observations"][witness["observation"]]
            bases, cell_labels, qualities_missing = {}, set(), False
            for record_id in observation["records"]:
                sam = details["original_records"][record_id]
                original_records[record_id] = sha256(sam.encode()).hexdigest()
                fields = sam.split("\t")
                flag = int(fields[1])
                cigar = re.findall(r"(\d+)([MIDNSHP=X])", fields[5])
                reverse = bool(flag & 16)
                sequence = str(Seq(fields[9]).reverse_complement()) if reverse else fields[9]
                hard = cigar[-1] if reverse else cigar[0]
                offset = int(hard[0]) if hard[1] == "H" else 0
                for query_position, base in enumerate(sequence, offset):
                    assert query_position not in bases or bases[query_position] == base
                    bases[query_position] = base
                tags = {v.split(":", 2)[0]: v.split(":", 2)[2] for v in fields[11:]}
                assert [tags.get("RG", ""), fields[0], flag & 192] == observation["identity"]
                if "CB" in tags:
                    cell_labels.add(tags["CB"])
                qualities_missing |= fields[10] == "*"
            start, end = witness["original_query_interval"]
            native_sequence = "".join(bases[q] for q in range(start, end))
            observed = (str(Seq(native_sequence).reverse_complement())
                        if observation["reverse_complement"] else native_sequence)
            assert observed == dna
            assert qualities_missing == observation["missing_qualities"]
            identity = tuple(observation["identity"])
            segments.add(identity)
            fragments.add(identity[:2])
            if qualities_missing:
                missing.add(identity)
            # These representatives have either no labels or unambiguous labels.
            # Do not invent a general conflict-resolution rule in this check.
            assert len(cell_labels) <= 1
            cells.update(cell_labels)
            witness_count += 1
    observed = dict(fragments=len(fragments), reads=len(segments), cells=len(cells),
                    missing_quality_reads=len(missing))
    assert observed == {key: candidate["rna_support"][key] for key in observed}
    return dict(source_id=candidate["source_id"], geometry_id=candidate["geometry_id"],
                orientation=candidate["orientation"], candidate_id=candidate["candidate_id"],
                amino_acids=candidate["amino_acids"], event_relations=candidate["event_relations"],
                independently_counted=observed, witness_intervals_checked=witness_count,
                reference_translations_checked=reference_count, original_record_hashes=original_records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("audit", type=Path)
    parser.add_argument("report", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    candidates = load(args.report / "candidate-orfs.json.gz")
    references = load(args.audit / "references.json.gz")
    checks = []
    for source_id in sorted({c["source_id"] for c in candidates}):
        candidate = max((c for c in candidates if c["source_id"] == source_id),
                        key=lambda c: c["rna_support"]["fragments"])
        details_path = args.audit / "reconstructions" / source_id / (
            candidate["geometry_id"] + "." + candidate["orientation"] + ".json.gz")
        checks.append(check_candidate(candidate, load(details_path), references))
    result = dict(check="highest_full_interval_support_per_pilot_product", checks=checks,
                  references_sha256=sha256((args.audit / "references.json.gz").read_bytes()).hexdigest(),
                  candidates_sha256=sha256((args.report / "candidate-orfs.json.gz").read_bytes()).hexdigest())
    with args.output.open("x") as handle:
        json.dump(result, handle, sort_keys=True, indent=2)
        handle.write("\n")
    print("Independently checked %d product representatives" % len(checks))


if __name__ == "__main__":
    main()
