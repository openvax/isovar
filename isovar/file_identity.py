"""Stream file identities without loading large original inputs into memory."""

from hashlib import sha256
from pathlib import Path


def file_sha256(path: str | Path) -> str:
    """Return the SHA-256 digest of exact file bytes, including compression."""
    digest = sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()
