#!/usr/bin/env bash
# Align a complete, validated HCC1395 ONT RNA SRA run to GRCh38.
# Retain only original-read alignments at source small variants and SV ends.
set -euo pipefail

if [ "$#" -ne 6 ]; then
    echo "usage: $0 GRCh38-no-alt.fa hcc1395-regions.bed run.sra fastq-dump output-dir direct_rna|cdna" >&2
    exit 2
fi
reference=$1
regions=$2
sra=$3
fastq_dump=$4
output_dir=$5
assay=$6
run=$(basename "$sra" .sra)
case "$assay" in
    direct_rna) minimap_args=(-ax splice -uf -k14) ;;
    cdna) minimap_args=(-ax splice -k14) ;;
    *) echo "unknown assay: $assay" >&2; exit 2 ;;
esac
mkdir -p "$output_dir"

"$(dirname "$fastq_dump")/vdb-validate" "$sra" \
    > "$output_dir/$run.validate.log" 2>&1
"$fastq_dump" -Z --split-spot "$sra" \
    2> "$output_dir/$run.dump.log" \
    | minimap2 -t 6 "${minimap_args[@]}" "$reference" - \
        2> "$output_dir/$run.align.log" \
    | samtools view -b -L "$regions" - \
    | samtools sort -@ 4 -m 256M -T "$output_dir/$run.sort" \
        -o "$output_dir/$run.all_sites.grch38.bam" -
samtools index -@ 4 "$output_dir/$run.all_sites.grch38.bam"
samtools flagstat "$output_dir/$run.all_sites.grch38.bam" \
    > "$output_dir/$run.all_sites.flagstat.txt"
