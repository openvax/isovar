"""Download public source regions for rebuild.py (manual, never run in CI)."""

import argparse
from hashlib import md5
from pathlib import Path
import subprocess
from urllib.request import urlretrieve

import pysam


GIAB = ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
        "NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark")
ONT = ("https://s3.amazonaws.com/nanopore-human-wgs/rna/bamFiles/"
       "NA12878-DirectRNA.pass.dedup.NoU.fastq.hg38.minimap2.sorted.bam")
PACBIO = "https://jbrowse.org/demos/cancer_sv/K562_isoseq.bam"


def fetch(destination):
    # Refuse to overwrite an existing source collection.
    destination.mkdir(parents=True, exist_ok=False)
    for resource, filename in ((8466, "star.bam"), (8467, "star.bam.bai")):
        urlretrieve(f"https://experimenthub.bioconductor.org/fetch/{resource}", destination / filename)
    assert md5((destination / "star.bam").read_bytes()).hexdigest() == "dc15bcd56e25e8ad09937fe64594438a"
    for source, region, filename in (
        (ONT, "chr6:43767094-43788458", "ont_chr6.sam"),
        (PACBIO, "chr22:23286000-23293000", "pacbio_k562.sam"),
    ):
        subprocess.run(["samtools", "view", "--no-PG", "-h", "-o",
                        str(destination / filename), source, region], check=True,
                       cwd=destination)
    index = destination / "giab.vcf.gz.tbi"
    urlretrieve(GIAB + ".vcf.gz.tbi", index)
    with pysam.TabixFile(GIAB + ".vcf.gz", index=str(index)) as vcf:
        for chrom, start, end in (("chr4", 0, 1000000), ("chr6", 43767093, 43788458)):
            records = list(vcf.fetch(chrom, start, end))
            (destination / (chrom + ".vcf")).write_text(
                "\n".join(list(vcf.header) + records) + "\n")
    urlretrieve(GIAB + ".bed", destination / "HG001.bed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path, help="New directory; existing directories are refused")
    fetch(parser.parse_args().destination.resolve())
