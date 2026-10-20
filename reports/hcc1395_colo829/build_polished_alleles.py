"""Audit COLO829 polished transcript bases at SMaHT somatic SNV sites.

This is sequence evidence from TAGET representative isoforms, not original
allele-spanning read support. A 0 requires a reference base at an aligned site;
an isoform with no coverage is never assigned 0.
"""

import argparse
from bisect import bisect_left
from collections import Counter, defaultdict
import csv
import gzip
from pathlib import Path

import pysam


FIELDS = [
    "variant_id", "chrom", "pos", "ref", "alt", "isoform_id", "allele",
    "observed_base", "mapping_quality", "tumor_full_length_count",
    "normal_full_length_count",
]
SUMMARY_FIELDS = [
    "variant_id", "chrom", "pos", "ref", "alt", "tumor_callability",
    "normal_callability", "tumor_alt_isoforms",
    "tumor_ref_isoforms", "tumor_other_isoforms", "normal_alt_isoforms",
    "normal_ref_isoforms", "normal_other_isoforms", "tumor_alt_fl_count",
    "tumor_ref_fl_count", "normal_alt_fl_count", "normal_ref_fl_count",
]


def read_variants(path, reference):
    """Return single-base PASS variants, verifying their reference allele."""
    by_chrom = defaultdict(dict)
    with gzip.open(path, "rt") as source:
        for line in source:
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alt, _, filt, *_ = line.rstrip("\n").split("\t")
            if filt != "PASS" or len(ref) != 1 or len(alt) != 1:
                continue
            pos = int(pos)
            if reference.fetch(chrom, pos - 1, pos).upper() != ref.upper():
                raise ValueError(f"Reference mismatch at {chrom}:{pos} {ref}>{alt}")
            if pos - 1 in by_chrom[chrom]:
                raise ValueError(f"Multiple alleles at {chrom}:{pos}")
            by_chrom[chrom][pos - 1] = (ref.upper(), alt.upper())
    return {chrom: (sorted(variants), variants) for chrom, variants in by_chrom.items()}


def read_counts(path):
    with path.open(newline="") as source:
        return {row["id"]: (int(row["COLO829T"]), int(row["COLO829BL"]))
                for row in csv.DictReader(source)}


def observations(alignment, variants, counts, min_mapq=20):
    """Yield only callable SNV/representative pairs from primary placements."""
    for read in alignment.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue
        chrom = alignment.get_reference_name(read.reference_id)
        if chrom not in variants:
            continue
        positions, alleles = variants[chrom]
        candidates = set()
        for start, end in read.get_blocks():
            left = bisect_left(positions, start)
            right = bisect_left(positions, end)
            candidates.update(positions[left:right])
        if not candidates:
            continue
        bases = {pos: read.query_sequence[index].upper()
                 for index, pos in read.get_aligned_pairs(matches_only=True)
                 if pos in candidates}
        isoform_id = read.query_name.split("|", 1)[0]
        tumor, normal = counts[isoform_id]
        for pos in sorted(candidates):
            base = bases.get(pos)
            if base is None:
                continue
            ref, alt = alleles[pos]
            allele = "1" if base == alt else "0" if base == ref else "?"
            yield {
                "variant_id": f"{chrom}:{pos + 1}:{ref}>{alt}", "chrom": chrom,
                "pos": pos + 1, "ref": ref, "alt": alt, "isoform_id": isoform_id,
                "allele": allele, "observed_base": base,
                "mapping_quality": read.mapping_quality,
                "tumor_full_length_count": tumor, "normal_full_length_count": normal,
            }


def write_outputs(rows, variants, output_dir):
    """Write sparse observations and one summary row for every source SNV."""
    output_dir.mkdir(parents=True, exist_ok=True)
    sums = defaultdict(Counter)
    seen = set()
    with gzip.open(output_dir / "colo829_polished_snv_alleles.tsv.gz", "wt", newline="") as target:
        writer = csv.DictWriter(target, FIELDS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            key = row["variant_id"], row["isoform_id"]
            if key in seen:
                raise ValueError(f"Duplicate placement: {key}")
            seen.add(key)
            writer.writerow(row)
            for sample, field in (("tumor", "tumor_full_length_count"),
                                  ("normal", "normal_full_length_count")):
                count = row[field]
                if count:
                    allele_name = {"1": "alt", "0": "ref", "?": "other"}[row["allele"]]
                    sums[row["variant_id"]][f"{sample}_{allele_name}_isoforms"] += 1
                    if allele_name in ("alt", "ref"):
                        sums[row["variant_id"]][f"{sample}_{allele_name}_fl_count"] += count
    with gzip.open(output_dir / "colo829_polished_snv_summary.tsv.gz", "wt", newline="") as target:
        writer = csv.DictWriter(target, SUMMARY_FIELDS, delimiter="\t")
        writer.writeheader()
        for chrom, (positions, alleles) in variants.items():
            for pos in positions:
                ref, alt = alleles[pos]
                variant_id = f"{chrom}:{pos + 1}:{ref}>{alt}"
                row = dict(variant_id=variant_id, chrom=chrom, pos=pos + 1,
                           ref=ref, alt=alt)
                for sample in ("tumor", "normal"):
                    callable_count = sum(sums[variant_id][f"{sample}_{allele}_isoforms"]
                                         for allele in ("alt", "ref", "other"))
                    row[f"{sample}_callability"] = "callable" if callable_count else "not_callable"
                row.update({field: sums[variant_id][field] for field in SUMMARY_FIELDS[7:]})
                writer.writerow(row)
    return len(seen), sum(len(positions) for positions, _ in variants.values())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    with pysam.FastaFile(str(args.reference)) as reference:
        variants = read_variants(args.vcf, reference)
    with pysam.AlignmentFile(str(args.bam), "rb") as alignment:
        totals = write_outputs(observations(alignment, variants, read_counts(args.counts)),
                               variants, args.output_dir)
    print(f"Callable isoform/SNV pairs: {totals[0]}; source SNVs: {totals[1]}")


if __name__ == "__main__":
    main()
