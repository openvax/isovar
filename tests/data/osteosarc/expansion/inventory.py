"""Snapshot public source metadata for the all-vaccine RNA audit (#218).

Acquisition is explicit and never runs in CI. Existing snapshots are checked
against their receipt, not silently refreshed. Parsing functions are offline.
"""

import argparse
from collections import defaultdict
import csv
from hashlib import sha256
from html.parser import HTMLParser
import io
import json
from pathlib import Path
import re
import subprocess


BUCKET = "https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/"
METADATA = {
    "variants.html": "https://osteosarc.com/variants/",
    "bams.json": "https://osteosarc.com/bams/bams.json",
    "vafs.tsv": "https://osteosarc.com/variants/variant_vafs_long.tsv",
    "vafs-columns.tsv": "https://osteosarc.com/variants/variant_vafs_long.columns.tsv",
    "data.html": "https://osteosarc.com/data/",
    "bucket_listing.json": "https://osteosarc.com/bucket_listing.json",
    "license.yaml": "https://raw.githubusercontent.com/awslabs/open-data-registry/main/datasets/sid-osteosarc.yaml",
}
RNA_CATEGORIES = {"RNA", "scRNA ONT", "PacBio", "Tumor scRNA", "Blood scRNA"}


def digest(path):
    """Hash files incrementally, including larger original regional sources."""
    value = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def write_json(path, value):
    """Create a deterministic generated artifact without overwriting work."""
    with Path(path).open("x") as handle:
        handle.write(json.dumps(value, indent=2, sort_keys=True) + "\n")


def fetch_snapshot(url, path):
    """Bounded HTTPS fetch with a fail-closed, independently hashed receipt."""
    path = Path(path)
    receipt = path.with_name(path.name + ".receipt.json")
    if path.exists():
        metadata = json.loads(receipt.read_text())
        if metadata["url"] != url or metadata["sha256"] != digest(path):
            raise ValueError(f"Snapshot drift: {path}")
        return metadata
    if receipt.exists():
        raise ValueError(f"Receipt without its snapshot: {path}")
    partial = path.with_name(path.name + ".partial")
    if partial.exists():
        raise ValueError(f"Unfinished download needs inspection: {partial}")
    subprocess.run([
        "curl", "--fail", "--location", "--silent", "--show-error",
        "--connect-timeout", "15", "--max-time", "120",
        "--retry", "3", "--retry-max-time", "180", "--output", str(partial), url,
    ], check=True, timeout=240)
    metadata = {"url": url, "sha256": digest(partial), "bytes": partial.stat().st_size}
    partial.rename(path)
    write_json(receipt, metadata)
    return metadata


class VariantIndex(HTMLParser):
    """Read only identity/membership cells, not tooltip or clinical claims."""

    def __init__(self):
        super().__init__()
        self.rows = []
        self.current = None
        self.cells = []
        self.cell = None

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "tr" and "data-vaccines" in attrs:
            self.current = {"vaccine_count": int(attrs["data-vaccines"])}
            self.cells = []
        elif self.current is not None:
            if tag == "td":
                self.cell = []
            elif tag == "a" and attrs.get("href", "").startswith("/variant/"):
                self.current["variant_id"] = attrs["href"].rstrip("/").split("/")[-1]

    def handle_data(self, data):
        if self.cell is not None:
            self.cell.append(data)

    def handle_endtag(self, tag):
        if tag == "td" and self.current is not None and self.cell is not None:
            self.cells.append(" ".join("".join(self.cell).split()))
            self.cell = None
        elif tag == "tr" and self.current is not None:
            if "variant_id" not in self.current or len(self.cells) != 9:
                raise ValueError("Variant index schema changed")
            self.rows.append(dict(self.current, gene=self.cells[0],
                                  source_location=self.cells[1], source_protein_label=self.cells[2]))
            self.current = None


def vaccine_variants(index_html, vafs_text):
    """Join exact variant IDs to unique original alleles; reject ambiguity."""
    index = VariantIndex()
    index.feed(index_html)
    selected = [r for r in index.rows if r["vaccine_count"] > 0]
    if not selected or len({v["variant_id"] for v in index.rows}) != len(index.rows):
        raise ValueError("Missing or duplicate variant identities in index")
    alleles = defaultdict(set)
    for row in csv.DictReader(io.StringIO(vafs_text), delimiter="\t"):
        alleles[row["variant_id"]].add((row["chrom"], int(row["pos"]), row["ref"], row["alt"]))
    result = []
    for record in selected:
        matches = alleles[record["variant_id"]]
        if len(matches) != 1:
            raise ValueError(f"Missing or ambiguous genomic allele for {record['variant_id']}: {matches}")
        chrom, pos, ref, alt = next(iter(matches))
        if (record["source_location"] != f"{chrom}:{pos}" or pos < 1
                or not re.fullmatch("[ACGT]+", ref) or not re.fullmatch("[ACGT]+", alt)):
            raise ValueError(f"Invalid/inconsistent genomic identity: {record}")
        result.append(dict(record, assembly="GRCh38", chrom=chrom, pos=pos, ref=ref, alt=alt,
                           source_url=f"https://osteosarc.com/variant/{record['variant_id']}/"))
    return result


def catalogue_sources(catalogue, vafs_text):
    """Preserve upstream assertions and conflicts instead of guessing samples."""
    labels = defaultdict(set)
    for row in csv.DictReader(io.StringIO(vafs_text), delimiter="\t"):
        labels[row["bam_file"]].add((row["sample_label"], row["timepoint"], row["sample_date"], row["assay_type"]))
    result = []
    for category in catalogue["categories"]:
        if category["name"] not in RNA_CATEGORIES:
            continue
        for record in category["bams"]:
            key = record["url"]
            result.append(dict(
                record, source_id=sha256(key.encode()).hexdigest()[:16], key=key,
                url=catalogue["baseUrl"] + key, category=category["name"],
                catalogue_genome_assertion=catalogue["genome"],
                count_export_claims=[dict(label=x[0], timepoint=x[1], date=x[2], assay=x[3])
                                     for x in sorted(labels[Path(key).name])],
                attribution=("pooled_library_unresolved" if category["name"] == "Blood scRNA"
                             else "source_label_only"),
                discovery="bam_catalogue",
            ))
    if len({r["key"] for r in result}) != len(result):
        raise ValueError("Duplicate catalogue alignment URLs")
    return result


def alignment_scope(key):
    """Classify inventory paths conservatively; headers decide coordinate space."""
    name = key.lower()
    if not name.endswith((".bam", ".cram")):
        return None
    if ("/vdj" in name or "tcr" in name or "bcr" in name
            or name.endswith(("/all_contig.bam", "/concat_ref.bam", "/consensus.bam"))
            or "/cider/" in name):
        return "excluded_receptor_or_targeted_alignment"
    if any(x in name for x in ("/dna/", "wes", "wgs", "/esvee/", "/teal/", "sarek", "vendor/natera/")):
        return "excluded_dna_alignment"
    if "totranscriptome" in name or name.endswith(".transcript.bam"):
        return "excluded_transcript_coordinate_alignment"
    if name.startswith("vendor/cegat/"):
        return ("rna_candidate" if key == "vendor/cegat/P116686_2_S000048/P116686_3.bam"
                else "excluded_dna_alignment")
    if (name.startswith(("rna-seq/", "ont/", "pacbio/", "hudson_lab/", "ucsf/", "t1/"))
            or "/rna/" in name or "rnaseq" in name
            or name.startswith(("kamil/blood/", "kamil/tumor/"))):
        return "rna_candidate"
    return "unresolved_assay"


def reconcile_bucket(catalogue_inventory, listing):
    """Retain every BAM/CRAM disposition and every distinct candidate product."""
    catalogue = {r["key"]: r for r in catalogue_inventory["catalogue_alignments"]}
    files = {r[0]: r for r in listing["files"]}
    if len(files) != len(listing["files"]):
        raise ValueError("Duplicate bucket keys")
    records = []
    for key in sorted(set(files) | set(catalogue)):
        scope = alignment_scope(key)
        if scope is None:
            continue
        if key in catalogue and scope != "rna_candidate":
            raise ValueError(f"Catalogue/path assay conflict: {key}")
        name = key.lower()
        row = dict(catalogue.get(key, {}))
        row.update(source_id=sha256(key.encode()).hexdigest()[:16], key=key, url=BUCKET + key,
                   scope=scope, listed_in_bucket=key in files,
                   bucket_bytes=files[key][1] if key in files else None,
                   bucket_mtime=files[key][2] if key in files else None,
                   listed_indexes=[x for x in (key + ".bai", key[:-4] + ".bai", key + ".csi", key + ".crai")
                                   if x in files])
        if key not in catalogue:
            row.update(discovery="bucket_listing", name=key,
                       tissue=("blood" if "hudson_lab/" in name or "kamil/blood/" in name else "source_path_only"),
                       attribution=("unassigned_library" if "unassigned_alignments" in name else
                                    "pooled_library_unresolved" if "hudson_lab/" in name or "kamil/blood/" in name
                                    else "source_path_only"))
        records.append(row)
    return dict(schema_version=1, bucket_generated_at=listing["generated_at"],
                metadata=catalogue_inventory["sources"], variants=catalogue_inventory["variants"],
                alignments=records)


def snapshot(destination):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    metadata = {}
    for name, url in METADATA.items():
        metadata[name] = fetch_snapshot(url, destination / name)
        print("snapshot", name, metadata[name]["bytes"], flush=True)
    variants = vaccine_variants((destination / "variants.html").read_text(),
                                (destination / "vafs.tsv").read_text())
    sources = catalogue_sources(json.loads((destination / "bams.json").read_text()),
                                (destination / "vafs.tsv").read_text())
    path = destination / "catalogue_inventory.json"
    value = dict(schema_version=1, sources=metadata, variants=variants, catalogue_alignments=sources)
    if path.exists():
        if json.loads(path.read_text()) != value:
            raise ValueError("Inventory changed; choose a new output directory")
    else:
        write_json(path, value)
    print(f"Pinned {len(variants)} vaccine variants and {len(sources)} catalogue RNA alignments", flush=True)
    reconciled = reconcile_bucket(value, json.loads((destination / "bucket_listing.json").read_text()))
    output = destination / "inventory.json"
    if output.exists():
        if json.loads(output.read_text()) != reconciled:
            raise ValueError("Reconciled inventory changed; choose a new output path")
    else:
        write_json(output, reconciled)
    from collections import Counter
    print("Alignment dispositions", dict(Counter(r["scope"] for r in reconciled["alignments"])), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    snapshot(parser.parse_args().destination)
