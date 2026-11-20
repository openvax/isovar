"""Read and write JSON files, gzip-compressed when the name ends in ``.gz``."""

import gzip
import json
import os
from pathlib import Path
import tempfile


def read_json(path, *, loader=json.load):
    """Load JSON from ``path``, decompressing a ``.gz`` file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(Path(path).expanduser(), "rt") as stream:
        return loader(stream)


def write_json(path, value, *, indent=2, sort_keys=False, canonical=False):
    """Write ``value`` to ``path`` atomically, gzip-compressed for a ``.gz`` name.

    A failed write leaves any existing file unchanged. When ``canonical`` is
    true, write compact, key-sorted, strict byte-stable JSON.
    """
    path = Path(path).expanduser()
    if canonical:
        payload = (
            json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False) + "\n"
        ).encode("utf-8")
    else:
        payload = (
            json.dumps(value, indent=indent, sort_keys=sort_keys, separators=None if indent else (",", ":")) + "\n"
        ).encode("utf-8")
    descriptor, temporary = tempfile.mkstemp(dir=path.parent, prefix="." + path.name, suffix=".tmp")
    try:
        with os.fdopen(descriptor, "wb") as stream:
            if str(path).endswith(".gz"):
                with gzip.GzipFile(filename="", mode="wb", fileobj=stream, mtime=0) as compressed:
                    compressed.write(payload)
            else:
                stream.write(payload)
        os.chmod(temporary, 0o666 & ~_umask())
        os.replace(temporary, path)
    except BaseException:
        os.unlink(temporary)
        raise


def publish_json(path, value, *, compressed=False):
    """Atomically publish canonical JSON, refusing to overwrite existing content."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = (
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False) + "\n"
    ).encode("utf-8")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
            temporary = Path(handle.name)
            if compressed:
                with gzip.GzipFile(filename="", mode="wb", fileobj=handle, mtime=0) as gz:
                    gz.write(payload)
            else:
                handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


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
