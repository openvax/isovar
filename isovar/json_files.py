"""Read and write JSON files, gzip-compressed when the name ends in ``.gz``."""

import gzip
import json
import os
from pathlib import Path
import tempfile


def read_json(path):
    """Load JSON from ``path``, decompressing a ``.gz`` file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(Path(path).expanduser(), "rt") as stream:
        return json.load(stream)


def write_json(path, value, *, indent=2, sort_keys=False):
    """Write ``value`` to ``path`` atomically, gzip-compressed for a ``.gz`` name.

    A failed write leaves any existing file unchanged.
    """
    path = Path(path).expanduser()
    text = json.dumps(value, indent=indent, sort_keys=sort_keys, separators=None if indent else (",", ":")) + "\n"
    descriptor, temporary = tempfile.mkstemp(dir=path.parent, prefix="." + path.name, suffix=".tmp")
    try:
        with os.fdopen(descriptor, "wb") as stream:
            if str(path).endswith(".gz"):
                with gzip.GzipFile(filename="", mode="wb", fileobj=stream, mtime=0) as compressed:
                    compressed.write(text.encode())
            else:
                stream.write(text.encode())
        os.chmod(temporary, 0o666 & ~_umask())
        os.replace(temporary, path)
    except BaseException:
        os.unlink(temporary)
        raise


def json_stem(path):
    """``path`` without a ``.json`` or ``.json.gz`` suffix, for files written beside it."""
    path = Path(path)
    for suffix in (".json.gz", ".json"):
        if path.name.endswith(suffix):
            return path.with_name(path.name[:-len(suffix)])
    return path


def _umask():
    mask = os.umask(0)
    os.umask(mask)
    return mask
