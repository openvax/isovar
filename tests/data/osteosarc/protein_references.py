"""Pin original Ensembl 87 records for six transcript-specific protein checks.

Run with the directory containing the original release-87 GTF, cDNA and
peptide FASTA archives and a NEW output directory. No Isovar/Varcode-derived
sequences enter these reference assets. CI reads the checked-in subset offline.
"""

import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import re


TRANSCRIPTS = {
    "PIP5K1A": "ENST00000368888",
    "MAP2": "ENST00000360351",
    "H1-2": "ENST00000343677",
    "EXOC4": "ENST00000253861",
    "GTF3C5": "ENST00000372108",
    "DYNC1H1": "ENST00000360184",
}
FILES = {
    "reference.gtf.gz": ("gtf", "Homo_sapiens.GRCh38.87.gtf.gz"),
    "reference.cdna.fa.gz": ("fasta", "Homo_sapiens.GRCh38.cdna.all.fa.gz"),
    "reference.pep.fa.gz": ("fasta", "Homo_sapiens.GRCh38.pep.all.fa.gz"),
}


def fasta_records(path):
    """Yield original FASTA headers and unwrapped sequences."""
    with gzip.open(path, "rt") as handle:
        header, sequence = None, []
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header, sequence = line.rstrip(), []
            else:
                sequence.append(line.strip())
        if header is not None:
            yield header, "".join(sequence)


def rebuild(source, output):
    output.mkdir(exist_ok=False)
    ids = set(TRANSCRIPTS.values())
    manifest = {"assembly": "GRCh38", "ensembl_release": 87,
                "transcripts": TRANSCRIPTS, "files": {}}
    genes = set()
    # Identify the six parent genes from the original cDNA headers.
    cdna = list((h, s) for h, s in fasta_records(source / FILES["reference.cdna.fa.gz"][1])
                if h[1:].split(".")[0] in ids)
    for header, _ in cdna:
        genes.add(re.search(r"gene:(ENSG\d+)", header)[1])
    assert len(cdna) == len(ids)
    for name, (kind, original) in FILES.items():
        path = source / original
        if kind == "gtf":
            with gzip.open(path, "rt") as handle:
                lines = []
                for line in handle:
                    if line.startswith("#"):
                        lines.append(line)
                        continue
                    transcript = re.search(r'transcript_id "([^"]+)"', line)
                    gene = re.search(r'gene_id "([^"]+)"', line)
                    if ((transcript and transcript[1] in ids)
                            or (line.split("\t")[2] == "gene" and gene[1] in genes)):
                        lines.append(line)
            data = "".join(lines).encode()
            url = "https://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/" + original
        else:
            records = cdna if ".cdna." in name else [
                (h, s) for h, s in fasta_records(path)
                if re.search(r"transcript:(ENST\d+)", h)[1] in ids]
            assert len(records) == len(ids)
            data = "".join(h + "\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60))
                           + "\n" for h, s in records).encode()
            url = ("https://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/"
                   + ("cdna/" if ".cdna." in name else "pep/") + original)
        with (output / name).open("wb") as handle:
            with gzip.GzipFile(filename="", mode="wb", fileobj=handle, mtime=0) as archive:
                archive.write(data)
        manifest["files"][name] = {
            "source_url": url, "source_sha256": sha256(path.read_bytes()).hexdigest(),
            "subset_sha256": sha256((output / name).read_bytes()).hexdigest(),
            "uncompressed_sha256": sha256(data).hexdigest(),
        }
    (output / "protein_reference_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rebuild(args.source, args.output)
