"""Known RNA paths on the actual GDC combined reference, through the public CLI.

No network access: supply the GDC FASTA and a new experiment directory. The
reads are synthetic truth; this is not a viral survey of biological samples.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from time import perf_counter

import pysam
from pyensembl import Genome

from isovar import __version__
from isovar.antigen_reference import gdc_reference_dictionary


def run(genome_fasta, output_dir):
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)
    reference = Path(genome_fasta).resolve()
    header = pysam.AlignmentHeader.from_dict(dict(
        SQ=[dict(SN=r["name"], LN=r["length"], M5=r["md5"]) for r in gdc_reference_dictionary()["contigs"]],
        RG=[dict(ID="library", SM="synthetic")]))
    cache = output / "reference-cache"
    with Genome(reference_name="GRCh38", annotation_name="unannotated",
                genome_fasta_path_or_url=str(reference), cache_directory_path=str(cache)) as genome:
        hpv = genome.sequence("HPV16", 1001, 1400)
        ebv = genome.sequence("chrEBV", 90001, 90150) + genome.sequence("chrEBV", 90301, 90450)
        host = genome.sequence("chr1", 100001, 100100)
    allele = hpv[:200] + ("A" if hpv[200] != "A" else "C") + hpv[201:]
    records, viral_names = [], set()

    def add(name, contig, positions, sequence, flag=0, mapq=60):
        cigar, previous = [], None
        for position in positions:
            if previous is not None and position != previous + 1:
                cigar.append((3, position - previous - 1))
            if cigar and cigar[-1][0] == 0:
                cigar[-1] = 0, cigar[-1][1] + 1
            else:
                cigar.append((0, 1))
            previous = position
        record = pysam.AlignedSegment(header)
        record.query_name, record.flag = name, flag
        record.reference_id, record.reference_start = header.get_tid(contig), positions[0]
        record.cigartuples, record.query_sequence = cigar, sequence
        record.query_qualities, record.mapping_quality = [35] * len(sequence), mapq
        record.set_tag("RG", "library")
        records.append(record)
        if contig != "chr1":
            viral_names.add(name)

    hpv_positions = list(range(1000, 1400))
    for start in range(0, 301, 50):
        for replicate in range(10):
            add("hpv-%d-%d" % (start, replicate), "HPV16", hpv_positions[start:start + 100], hpv[start:start + 100],
                flag=16 if replicate % 2 else 0)
    for i in range(3):
        add("allele-%d" % i, "HPV16", hpv_positions, allele)
    ebv_positions = list(range(90000, 90150)) + list(range(90300, 90450))
    for i in range(3):
        add("ebv-full-%d" % i, "chrEBV", ebv_positions, ebv)
    for start in range(0, 201, 50):
        for i in range(2):
            add("ebv-short-%d-%d" % (start, i), "chrEBV", ebv_positions[start:start + 100], ebv[start:start + 100])
    add("host-or-virus", "HPV16", hpv_positions[:100], hpv[:100], mapq=0)
    add("host-or-virus", "chr1", list(range(120000, 120100)), hpv[:100], flag=256, mapq=0)
    for i in range(20):
        add("host-%d" % i, "chr1", list(range(100000, 100100)), host)
    bam_path = output / "reads.bam"
    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for record in records:
            bam.write(record)
    command = [sys.executable, "-m", "isovar", "viral", "--bam", str(bam_path),
               "--genome-fasta", str(reference), "--sample-id", "synthetic", "--output-dir", str(output / "result"),
               "--reference-cache", str(cache)]
    started = perf_counter()
    with (output / "cli.log").open("w") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    elapsed = perf_counter() - started
    result = json.loads((output / "result/viral.json").read_text())
    assert len(result["selected_contigs"]) == 200
    assert {r["sequence"] for r in result["sequences"]} == {hpv, allele, ebv}
    assert result["support"]["target_fragments"] == len(viral_names) == 87
    assert result["support"]["unassigned_target_fragments"] == 0
    assert sum(c["fragments"] for c in result["support"]["equivalence_classes"]) == 87
    assert result["input"]["records"] == 108
    ambiguous = next(w for w in result["support"]["witnesses"] if w["fragment"][1] == "host-or-virus")
    assert ambiguous["competing_placement_outside_targets"]
    assert not result["discovery"]["truncated"]
    with pysam.FastxFile(str(output / "result/sequences.fa")) as fasta:
        assert {r.sequence for r in fasta} == {hpv, allele, ebv}
    record = dict(isovar_version=__version__, pysam_version=pysam.__version__,
                  reference=result["reference"], input_records=108, viral_fragments=87,
                  selected_contigs=200, observed_sequences=3, unassigned_fragments=0,
                  reconstructed_sequence_sha256s=sorted(sha256(s.encode()).hexdigest() for s in (hpv, allele, ebv)),
                  cli_seconds=elapsed, command=command,
                  observation="Synthetic exact reads on the actual GDC reference; no GTF or mutation input",
                  script_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    (output / "experiment.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({k: v for k, v in record.items() if k not in {"reference", "command"}}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-fasta", required=True)
    parser.add_argument("--output-dir", required=True)
    options = parser.parse_args()
    run(options.genome_fasta, options.output_dir)
