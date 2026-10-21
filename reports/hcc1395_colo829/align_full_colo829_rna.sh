#!/usr/bin/env bash
# Align the complete, validated SRR8615617 SRA run to GRCh38.
# Only alignments at nominated small-variant loci and SV ends are retained.
set -euo pipefail

if [ "$#" -ne 5 ]; then
    echo "usage: $0 GRCh38-no-alt.fa atlas-regions.bed SRR8615617.sra fastq-dump output-dir" >&2
    exit 2
fi
reference=$1
regions=$2
sra=$3
fastq_dump=$4
output_dir=$5
mkdir -p "$output_dir"

"$(dirname "$fastq_dump")/vdb-validate" "$sra" \
    > "$output_dir/SRR8615617.validate.log" 2>&1
"$fastq_dump" -Z --split-spot "$sra" \
    2> "$output_dir/SRR8615617.dump.log" \
    | minimap2 -t 6 -ax sr --frag=yes "$reference" - \
        2> "$output_dir/SRR8615617.align.log" \
    | samtools view -b -L "$regions" - \
    | samtools sort -@ 4 -m 256M -T "$output_dir/sort" \
        -o "$output_dir/SRR8615617.all_sites.grch38.bam" -
samtools index -@ 4 "$output_dir/SRR8615617.all_sites.grch38.bam"
samtools flagstat "$output_dir/SRR8615617.all_sites.grch38.bam" \
    > "$output_dir/SRR8615617.all_sites.flagstat.txt"
