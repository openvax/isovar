"""Regenerate the three validation alleles from checksum-pinned public VCFs.

python -m tests.data.fusions.build_validation_events --snapshot NAME --output events.json
Acquisition uses the shared Osteosarc cache. --source-directory DIR instead
reads original compressed VCFs named by their URL basename, entirely offline.
The checked-in manifest is the selection recipe (source, SHA256, record ID).
"""
import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
from urllib.parse import urlsplit

from isovar.sid_data import open_dataset

RECIPE = Path(__file__).with_name("validation-events.json")


def regenerate(recipe, acquire):
    result = []
    for entry in recipe:
        old = entry["record"]
        path = acquire(old["source_url"])
        if sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            raise ValueError("Source checksum mismatch: " + old["source_url"])
        with gzip.open(path, "rt") as handle:
            matches = [line.rstrip("\n").split("\t") for line in handle
                       if not line.startswith("#") and line.split("\t")[2] == old["record_id"]]
        if len(matches) != 1:
            raise ValueError("Expected one record: " + old["record_id"])
        chrom, pos, identifier, ref, alt, _, _, info = matches[0][:8]
        if "," in alt:
            raise ValueError("Selection requires an explicit ALT index")
        fields = dict(item.split("=", 1) for item in info.split(";") if "=" in item)
        record = dict(old, chrom=chrom, pos=pos, record_id=identifier, ref=ref,
                      alt=alt, info=info, end=fields.get("END", pos))
        result.append(dict(entry, record=record))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot")
    parser.add_argument("--cache")
    parser.add_argument("--source-directory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.source_directory:
        acquire = lambda url: args.source_directory / Path(urlsplit(url).path).name
    else:
        dataset = open_dataset(args.snapshot, args.cache)
        acquire = lambda url: dataset.download(dataset.asset(url))
    output = regenerate(json.loads(RECIPE.read_text()), acquire)
    with args.output.open("x") as handle:
        json.dump(output, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
