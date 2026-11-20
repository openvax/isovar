"""Stream file identities without loading large original inputs into memory."""

from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path


def file_sha256(path: str | Path) -> str:
    """Return the SHA-256 digest of exact file bytes, including compression."""
    digest = sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def software_identity() -> dict:
    """Record the reconstruction engine and its native dependencies."""
    from . import __version__

    root = Path(__file__).parent
    files = {p.relative_to(root).as_posix(): file_sha256(p) for p in sorted(root.rglob("*.py"))}
    engine_sha256 = sha256(
        json.dumps(files, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode("utf-8")
    ).hexdigest()
    return dict(
        isovar_version=__version__,
        engine_sha256=engine_sha256,
        dependencies={name: version(name) for name in ("pysam", "edlib")},
    )
