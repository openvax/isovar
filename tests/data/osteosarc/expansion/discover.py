"""Reconcile the dated website snapshot against live public RNA-related prefixes."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from hashlib import sha256
import json
from pathlib import Path
import sys
from urllib.parse import urlencode
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from tests.data.osteosarc.expansion.inventory import (  # noqa: E402
    BUCKET, fetch_snapshot, reconcile_bucket, write_json,
)


# Published RNA paths, plus analysis/spatial prefixes that may contain
# additional genomic RNA products. Imaging/reference-only roots are outside
# scope; their previously listed BAM dispositions remain in the inventory.
PREFIXES = (
    "ONT/", "T1/", "ctat_lr_fusion/", "elucidate/", "genomics/",
    "hms_spatial/", "hudson_lab/", "kamil/", "pacbio/", "rna-seq/",
    "scratch/", "ucsf/", "vendor/", "genomics_reprocessing/", "neoantigen_prediction/", "web/", "webpage/",
)
NS = {"s": "http://s3.amazonaws.com/doc/2006-03-01/"}


def list_prefix(prefix, directory):
    """Paginate explicitly; a truncated or failed listing is never complete."""
    records, receipts, token = [], [], None
    seen_tokens = set()
    for page in range(1000):
        query = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
        if token is not None:
            query["continuation-token"] = token
        name = sha256(prefix.encode()).hexdigest()[:16] + f"-{page:04d}.xml"
        path = Path(directory) / name
        receipts.append(fetch_snapshot(BUCKET + "?" + urlencode(query), path))
        root = ET.fromstring(path.read_bytes())
        if root.findtext("s:Prefix", namespaces=NS) != prefix:
            raise ValueError("S3 listing prefix mismatch")
        for item in root.findall("s:Contents", NS):
            key = item.findtext("s:Key", namespaces=NS)
            if not key.startswith(prefix):
                raise ValueError("S3 key escaped requested prefix")
            modified = item.findtext("s:LastModified", namespaces=NS)
            records.append([key, int(item.findtext("s:Size", namespaces=NS)),
                            int(datetime.fromisoformat(modified.replace("Z", "+00:00")).timestamp())])
        truncated = root.findtext("s:IsTruncated", namespaces=NS)
        if truncated == "false":
            print("live inventory", prefix, len(records), "objects", flush=True)
            return records, receipts
        token = root.findtext("s:NextContinuationToken", namespaces=NS)
        if truncated != "true" or not token or token in seen_tokens:
            raise ValueError("Invalid/repeated S3 pagination token")
        seen_tokens.add(token)
    raise ValueError("S3 prefix exceeded the bounded page limit")


def discover(destination, suffix=""):
    destination = Path(destination)
    directory = destination / "live-listings"
    directory.mkdir(exist_ok=True)
    root_receipt = fetch_snapshot(BUCKET + "?list-type=2&delimiter=%2F&max-keys=1000", directory / "root.xml")
    root = ET.fromstring((directory / "root.xml").read_bytes())
    if root.findtext("s:IsTruncated", namespaces=NS) != "false":
        raise ValueError("Root listing is incomplete")
    roots = [x.findtext("s:Prefix", namespaces=NS) for x in root.findall("s:CommonPrefixes", NS)]
    snapshot = json.loads((destination / "bucket_listing.json").read_text())
    kept = [r for r in snapshot["files"] if not r[0].startswith(PREFIXES)]
    receipts = [root_receipt]
    with ThreadPoolExecutor(max_workers=4) as executor:
        for records, pages in executor.map(lambda p: list_prefix(p, directory), PREFIXES):
            kept.extend(records)
            receipts.extend(pages)
    snapshot["files"] = kept
    snapshot["generated_at"] = "live S3 prefix snapshots; see checksummed per-page receipts"
    base = json.loads((destination / "catalogue_inventory.json").read_text())
    reconciled = reconcile_bucket(base, snapshot)
    reconciled["live_discovery"] = dict(prefixes=PREFIXES, root_prefixes=roots, receipts=receipts,
                                       prior_bucket_generated_at=json.loads((destination / "bucket_listing.json").read_text())["generated_at"])
    write_json(destination / f"inventory-live{suffix}.json", reconciled)
    # Keep relevant auxiliary metadata keys for attribution/reference work,
    # without retaining imaging filenames in the checked-in report.
    metadata = [r for r in kept if r[0].startswith(PREFIXES)
                and any(token in r[0].lower() for token in ("config", "sample", "donor", "demultiplex", "manifest"))
                and r[0].endswith((".csv", ".json", ".tsv", ".yaml", ".yml"))]
    write_json(destination / f"auxiliary-metadata-inventory{suffix}.json", metadata)
    print("Live inventory complete", len(reconciled["alignments"]), "alignment products", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--suffix", default="")
    args = parser.parse_args()
    if any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_" for c in args.suffix):
        parser.error("Invalid snapshot suffix")
    discover(args.destination, args.suffix)
