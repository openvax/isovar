"""Check TAGET COLO829 polished representatives for the exact SV3 RNA join.

The 450 MB source FASTQ and tumor GFF are downloaded from Zenodo 10091914.
They are verified against the source hashes in the atlas README before use.
"""

import argparse
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import re

from Bio.Seq import Seq


HERE = Path(__file__).parent
COUNTS = HERE / "colo829_polished_isoform_counts.tsv.gz"
ASSESSMENT = HERE / "colo829_sv3_candidate_orfs.assessment.json"
REFERENCE = HERE.parents[1] / "tests/data/colo829_sv3/reference_slice.json"
OUTPUT_PREFIX = HERE / "colo829_sv3_polished_check"
SOURCE_URL = "https://zenodo.org/records/10091914"
FASTQ_SHA256 = "6b008d7610c03bdcc62e716bd365f4b448badc9c77554c58574068c7c9ab487b"
GFF_SHA256 = "506488c868150e35358514fd1b9528ca4070ed0f4297e50d42364ad6594f1b94"


def _digest(path):
    checksum = sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(chunk)
    return checksum.hexdigest()


def _motifs(assessment):
    join = assessment["unplaced_junction_bases"]
    motifs = dict(junction_plus=join, junction_reverse_complement=str(Seq(join).reverse_complement()))
    for candidate in assessment["candidates"]:
        for peptide, dna in candidate["junction_27nt_windows"].items():
            motifs["window_plus_" + peptide] = dna
            motifs["window_reverse_complement_" + peptide] = str(Seq(dna).reverse_complement())
    return motifs


def scan_sequences(fastq, motifs):
    """Count exact motifs on complete polished representative sequences."""
    hits = {key: [] for key in motifs}
    count = 0
    with Path(fastq).open() as stream:
        while header := stream.readline():
            sequence, plus, qualities = stream.readline(), stream.readline(), stream.readline()
            if not qualities or not header.startswith("@") or not plus.startswith("+"):
                raise ValueError("Incomplete or invalid polished FASTQ")
            sequence = sequence.strip().upper()
            if len(sequence) != len(qualities.strip()):
                raise ValueError("Polished FASTQ sequence and quality lengths differ")
            name = header[1:].split("|")[0].split()[0]
            count += 1
            for key, motif in motifs.items():
                if motif in sequence:
                    hits[key].append(name)
    return count, {key: sorted(set(names)) for key, names in hits.items()}


def _counts(path):
    with gzip.open(path, "rt") as stream:
        return {row["isoform_id"]: row for row in csv.DictReader(stream, delimiter="\t")}


def _exons(gff, names):
    exons = {name: [] for name in names}
    with Path(gff).open() as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            match = re.search(r'transcript_id "([^"]+)"', fields[8])
            if match and match.group(1) in names:
                exons[match.group(1)].append([fields[0], int(fields[3]), int(fields[4])])
    return exons


def assess(fastq, gff, counts=COUNTS, assessment_path=ASSESSMENT, reference_path=REFERENCE,
           validate_hashes=True):
    """Audit exact source-sequence support without treating motifs as isoforms."""
    fastq, gff = Path(fastq), Path(gff)
    fastq_hash, gff_hash = _digest(fastq), _digest(gff)
    if validate_hashes and (fastq_hash != FASTQ_SHA256 or gff_hash != GFF_SHA256):
        raise ValueError("TAGET FASTQ or tumor GFF differs from the pinned source")
    assessment = json.loads(Path(assessment_path).read_text())
    reference = json.loads(Path(reference_path).read_text())
    motifs = _motifs(assessment)
    count, hits = scan_sequences(fastq, motifs)
    cluster_counts = _counts(counts)
    if count != len(cluster_counts):
        raise ValueError("Polished representative count differs from the count table")
    genomic_span_count = sum(row["chrom"] == "chr1" and int(row["start"]) <= 224595102
                             and int(row["end"]) >= 224598323 for row in cluster_counts.values())
    hit_ids = sorted({name for names in hits.values() for name in names})
    if any(name not in cluster_counts for name in hit_ids):
        raise ValueError("A motif-positive representative lacks a count-table row")
    exons = _exons(gff, set(hit_ids))
    motif_rows = [dict(isoform_id=name,
                       tumor_full_length_count=int(cluster_counts[name]["tumor_full_length_count"]),
                       normal_full_length_count=int(cluster_counts[name]["normal_full_length_count"]),
                       genomic_span=[cluster_counts[name]["chrom"],
                                     int(cluster_counts[name]["start"]), int(cluster_counts[name]["end"])],
                       matched_motifs=[key for key, names in hits.items() if name in names],
                       tumor_gff_exons=exons[name]) for name in hit_ids]
    exact_windows = {key: names for key, names in hits.items() if key.startswith("window_") and names}
    status = ("candidate_exact_junction_window_present" if exact_windows else
              "no_exact_junction_window_in_polished_representatives")
    reverse = motifs["junction_reverse_complement"]
    positions = [reference["start"] + i + 1 for i in range(len(reference["sequence"]))
                 if reference["sequence"].startswith(reverse, i)]
    return dict(schema="isovar.colo829_sv3_polished_check.v1", event_id=assessment["event_id"],
                source_url=SOURCE_URL, source_fastq_sha256=fastq_hash, source_gff_sha256=gff_hash,
                representative_count=count, source_culture="TAGET_COLO829_separate_from_SMaHT_DNA_and_SRR8615617",
                dna_breakpoint_coordinate_spanning_representatives=genomic_span_count,
                exact_motif_hits=hits, motif_positive_representatives=motif_rows,
                reverse_motif_reference_positions_grch38_1based=positions,
                exact_polished_join_status=status, proteoform_call="not_established")


def _report(result):
    hits = result["exact_motif_hits"]
    window_count = sum(len(names) for key, names in hits.items() if key.startswith("window_"))
    lines = ["# COLO829 SV3: polished-isoform sequence cross-check", "",
             "**No exact polished SV3 join or resolved proteoform.** Across",
             f"{result['representative_count']:,} TAGET full-length representative sequences,",
             f"the forward 19-base junction motif occurs in {len(hits['junction_plus'])} and",
             f"the five strand-aware 27-base junction windows have {window_count} matches.", "",
             f"{result['dna_breakpoint_coordinate_spanning_representatives']} representatives have a genomic",
             "span across both nominated DNA breakend coordinates, but a genomic",
             "span does not establish an exon or RNA junction at either end.", "",
             "The reverse complement of the 19-base motif occurs in two",
             "representatives. It also occurs in the ordinary GRCh38 sequence at",
             "chr1:" + ", ".join(map(str, result["reverse_motif_reference_positions_grch38_1based"])) +
             ". A reverse-motif hit alone is not an SV join.", "",
             "| Representative | Tumor full-length count | Normal full-length count | "
             "Matched motif | Tumor GFF exons |",
             "|---|---:|---:|---|---|"]
    for row in result["motif_positive_representatives"]:
        exons = "; ".join(f"{contig}:{start}-{end}" for contig, start, end in row["tumor_gff_exons"])
        lines.append(f"| {row['isoform_id']} | {row['tumor_full_length_count']} | "
                     f"{row['normal_full_length_count']} | "
                     f"{', '.join(row['matched_motifs'])} | {exons} |")
    lines += ["", "Both reverse-motif representatives have a final exon in the",
              "downstream reference block and no exon at the nominated upstream",
              "SV3 breakend. Their tumor cluster counts cannot be assigned to",
              "the SV3 junction. The separate-culture polished data cannot",
              "exclude a rare or unpolished SV transcript.", "",
              f"Sources: [TAGET Zenodo record]({result['source_url']}) and the",
              "[source-pinned atlas](README.md#input-sha-256). Regenerate with",
              "python -m reports.hcc1395_colo829.build_colo829_sv3_polished_check",
              "--fastq /path/to/final.rep.fq --gff /path/to/final.COLO829T.gff.", ""]
    return "\n".join(lines)


def write_audit(fastq, gff, prefix=OUTPUT_PREFIX):
    """Write the sequence cross-check and a concise human-readable summary."""
    result = assess(fastq, gff)
    prefix = Path(prefix)
    json_path = Path(str(prefix) + ".json")
    json_path.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    Path(str(prefix) + ".md").write_text(_report(result))
    return json_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, default=OUTPUT_PREFIX)
    args = parser.parse_args()
    print(write_audit(args.fastq, args.gff, args.output_prefix))


if __name__ == "__main__":
    main()
