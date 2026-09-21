"""Content identity for checksum-pinned offline annotation datasets."""

import gzip
from hashlib import sha256
import json
from pathlib import Path


REFERENCE_FILES = ("reference.gtf.gz", "reference.cdna.fa.gz", "reference.pep.fa.gz")


def reference_dataset_identity(directory, checksums):
    """Verify pinned archives and fingerprint their uncompressed contents.

    Parameters
    ----------
    directory : path-like
        Directory containing the GTF, transcript and protein reference files.
    checksums : mapping of str to str
        Expected SHA256 of each compressed asset, keyed by filename.

    Returns
    -------
    str
        Dataset name independent of directory, gzip metadata and descriptive
        manifest fields. File roles participate in the fingerprint.
    """
    contents = {}
    for name in REFERENCE_FILES:
        if name not in checksums:
            raise ValueError(f"Missing reference asset checksum: {name}")
        data = (Path(directory) / name).read_bytes()
        if sha256(data).hexdigest() != checksums[name]:
            raise ValueError(f"Reference asset checksum mismatch: {name}")
        contents[name] = sha256(gzip.decompress(data)).hexdigest()
    payload = json.dumps(contents, sort_keys=True, separators=(",", ":")).encode()
    return "reference-sha256-" + sha256(payload).hexdigest()
