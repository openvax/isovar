"""Pin vaccine membership, source claims and library-processing metadata."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from hashlib import sha256
from html.parser import HTMLParser
import io
import json
from pathlib import Path
import re
import sys
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from tests.data.osteosarc.expansion.inventory import BUCKET, digest, fetch_snapshot, write_json  # noqa: E402


class PageSections(HTMLParser):
    """Retain only genomic/protein definitions and vaccine inclusion text."""

    def __init__(self):
        super().__init__()
        self.heading = None
        self.current = None
        self.sections = {}
        self.ignore = 0

    def handle_starttag(self, tag, attrs):
        if tag in ("script", "style"):
            self.ignore += 1
        if tag == "h2":
            self.heading = []

    def handle_endtag(self, tag):
        if tag in ("script", "style"):
            self.ignore -= 1
        if tag == "h2":
            self.current = " ".join("".join(self.heading).split())
            self.sections[self.current] = []
            self.heading = None

    def handle_data(self, value):
        if self.ignore:
            return
        if self.heading is not None:
            self.heading.append(value)
        elif self.current in ("Genomic context", "Protein context", "Vaccines targeting this variant"):
            self.sections[self.current].append(value)

    def result(self):
        return {k: " ".join(" ".join(v).split()) for k, v in self.sections.items() if v}


def variant_details(record, directory):
    path = directory / (record["variant_id"] + ".html")
    receipt = fetch_snapshot(record["source_url"], path)
    parser = PageSections()
    parser.feed(path.read_text())
    sections = parser.result()
    if "Vaccines targeting this variant" not in sections:
        raise ValueError("Vaccine membership section missing from selected variant")
    return dict(variant_id=record["variant_id"], receipt=receipt, source_sections=sections,
                source_refseq_ids=sorted(set(re.findall(r"\bNM_\d+(?:\.\d+)?\b", path.read_text()))))


def cellranger_config(text):
    """Preserve explicit sample/library tables; do not infer donor identities."""
    sections, current = {}, None
    for row in csv.reader(io.StringIO(text)):
        if not row or not row[0].strip() or row[0].startswith("#"):
            continue
        if len(row) == 1 and row[0].startswith("[") and row[0].endswith("]"):
            current = row[0][1:-1]
            sections[current] = []
        elif current is not None:
            # Analysis service paths/token-file locations are not sample
            # identity evidence and are not needed in the checked-in report.
            if "token" not in row[0].lower():
                sections[current].append(row)
    return sections


def acquire_config(record, directory):
    key = record[0]
    path = directory / (sha256(key.encode()).hexdigest()[:16] + ".csv")
    receipt = fetch_snapshot(BUCKET + quote(key, safe="/"), path)
    return dict(key=key, receipt=receipt, sections=cellranger_config(path.read_text()))


def prepare(destination):
    destination = Path(destination)
    directory = destination / "source-metadata"
    directory.mkdir(exist_ok=True)
    inventory = json.loads((destination / "inventory-live.json").read_text())
    auxiliary = json.loads((destination / "auxiliary-metadata-inventory.json").read_text())
    configs = [r for r in auxiliary if r[0].endswith("/config.csv")
               and r[0].startswith(("hudson_lab/", "kamil/", "T1/")) and r[1] < 100000]
    with ThreadPoolExecutor(max_workers=4) as executor:
        details = list(executor.map(lambda v: variant_details(v, directory), inventory["variants"]))
        libraries = list(executor.map(lambda r: acquire_config(r, directory), configs))
    checksums = fetch_snapshot(BUCKET + "md5sums.txt", directory / "md5sums.txt")
    # These are published source MD5 claims, not local verification that the
    # full remote objects were downloaded. Multipart ETags are not MD5s.
    md5s = {}
    for line in (directory / "md5sums.txt").read_text().splitlines():
        fields = line.split(maxsplit=1)
        if len(fields) == 2 and re.fullmatch("[a-fA-F0-9]{32}", fields[0]):
            md5s[fields[1].lstrip("*").removeprefix("./")] = fields[0].lower()
    value = dict(source_inventory_sha256=digest(destination / "inventory-live.json"),
                 variants=details, library_configs=libraries, published_checksum_receipt=checksums,
                 alignment_md5_claims={s["source_id"]: md5s[s["key"]] for s in inventory["alignments"] if s["key"] in md5s})
    output = destination / "source-metadata.json"
    if output.exists():
        if json.loads(output.read_text()) != value:
            raise ValueError("Source metadata drift; choose a new snapshot")
    else:
        write_json(output, value)
    print("Pinned variant pages", len(details), "library configs", len(libraries), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    prepare(parser.parse_args().destination)
