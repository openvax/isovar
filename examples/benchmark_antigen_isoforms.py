#!/usr/bin/env python
"""Reproduce two synthetic transcript-mixture experiments with real Salmon.

Run from a repository checkout. Synthetic references deliberately exercise
different transcript lengths and an alternative cassette exon. They are not
measurements from patients or evidence that any gene is a cancer-testis antigen.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random
import time

import pysam

from isovar import reconstruct_antigen_isoforms, quantify_antigen_isoforms
from isovar.dna import reverse_complement_dna
from tests.test_antigen_isoforms import target_models
from tests.test_sv_rna import aligned, write_bam


def experiment(output, salmon, fraction, seed):
    output.mkdir()
    rng = random.Random(seed)
    targets, paths = target_models()
    sequences = [t.reference.sequence for t in targets]
    transcriptome = output / "background.fa"
    transcriptome.write_text("".join(">%s\n%s\n" % (t.reference.transcript_id, t.reference.sequence) for t in targets))
    first, second = output / "reads1.fastq", output / "reads2.fastq"
    records, truth_counts = [], [0, 0]
    weights = [fraction * (len(sequences[0]) - 200 + 1), (1 - fraction) * (len(sequences[1]) - 200 + 1)]
    with first.open("w") as f1, second.open("w") as f2:
        for i in range(2500):
            k = 0 if rng.random() < weights[0] / sum(weights) else 1
            truth_counts[k] += 1
            length = min(len(sequences[k]), max(160, round(rng.gauss(200, 15))))
            start = rng.randrange(len(sequences[k]) - length + 1)
            end = start + length
            name = "fragment%05d" % i
            r1, r2 = sequences[k][start:start + 100], sequences[k][end - 100:end]
            f1.write("@%s/1\n%s\n+\n%s\n" % (name, r1, "I" * 100))
            f2.write("@%s/2\n%s\n+\n%s\n" % (name, reverse_complement_dna(r2), "I" * 100))
            records.extend(aligned(name, r1, paths[k][start:start + 100], flag=65))
            records.extend(aligned(name, r2, paths[k][end - 100:end], flag=145))
    bam_path = write_bam(output / "reads.bam", records)
    started = time.monotonic()
    with pysam.AlignmentFile(bam_path) as bam:
        result = reconstruct_antigen_isoforms(bam, targets=targets, sample_id=output.name, source=str(bam_path),
                                             max_discovery_reads=6000, scratch_dir=output)
    reconstruction_seconds = time.monotonic() - started
    assert result["support"]["target_fragments"] == 2500
    assert result["support"]["unassigned_target_fragments"] == 0
    assert result["support"]["fragments_without_full_transcript_assignment"] == 0
    assert not result["discovery"]["truncated"]
    assert {r["sequence"] for r in result["isoforms"]} == set(sequences)
    started = time.monotonic()
    quantify_antigen_isoforms(result, transcriptome=transcriptome, reads1=first, reads2=second,
                             output_dir=output / "salmon", salmon=salmon)
    quantification_seconds = time.monotonic() - started
    gene, = result["quantification"]["genes"]
    isoform = next(row for row in result["isoforms"] if "T1" in row["annotated_transcript_ids"])
    estimate = gene["isoform_fractions"][isoform["id"]]
    assert abs(estimate - fraction) < 0.04, (fraction, estimate)
    dominant = next(r for r in result["isoforms"] if r["id"] == gene["dominant_isoform_id"])
    assert dominant["annotated_transcript_ids"] == ["T1" if fraction > 0.5 else "T2"]
    return dict(sample_id=output.name, expected_T1_fraction=fraction, estimated_T1_fraction=estimate,
                original_fragments=truth_counts, counted_fragments=result["support"]["target_fragments"],
                unassigned=result["support"]["unassigned_target_fragments"],
                dominant_transcripts=dominant["annotated_transcript_ids"],
                reconstruction_seconds=reconstruction_seconds, quantification_seconds=quantification_seconds,
                software=result["quantification"]["provenance"]["version"],
                inputs=result["quantification"]["provenance"]["inputs"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--salmon", required=True)
    args = parser.parse_args()
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)
    salmon = str(Path(args.salmon).resolve())
    result = dict(kind="synthetic_known_mixtures", seed=465,
                  script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  experiments=[experiment(output / name, salmon, fraction, seed)
                               for name, fraction, seed in [("sample-A", 0.8, 465), ("sample-B", 0.2, 466)]])
    (output / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
