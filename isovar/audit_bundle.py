"""Portable, checksum-verified files for offline scientific audits.

Bundles contain data and structured metadata, never executable replay commands.
The manifest format and the public API are documented in ``docs/audit-bundles.md``.
"""

from copy import deepcopy
from collections.abc import Mapping
from hashlib import sha256
from importlib.metadata import version
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile
import tempfile


FORMAT = "isovar-audit-bundle"
SCHEMA_VERSION = 1
MANIFEST = "manifest.json"
__all__ = ["AuditBundle", "AuditBundleError"]


class AuditBundleError(ValueError):
    """An audit bundle, input or replay checkpoint failed validation."""


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def identity(value):
    return sha256(canonical(value)).hexdigest()


def digest(path):
    result = sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def software_identity():
    """Record the reconstruction engine and its native dependencies."""
    from . import __version__

    root = Path(__file__).parent
    files = {p.relative_to(root).as_posix(): digest(p) for p in sorted(root.rglob("*.py"))}
    return dict(isovar_version=__version__, engine_sha256=identity(files),
                dependencies={name: version(name) for name in ("pysam", "edlib")})


def _name(value):
    if not isinstance(value, str) or not value or any(c in value for c in ("\\", ":", "\0")):
        raise AuditBundleError("File names must be nonempty relative POSIX paths")
    path = PurePosixPath(value)
    if path.is_absolute() or any(p in (".", "..") for p in value.split("/")) or path.as_posix() != value:
        raise AuditBundleError("Invalid relative file name: " + value)
    return value


def _validate(manifest):
    if not isinstance(manifest, dict) or manifest.get("format") != FORMAT:
        raise AuditBundleError("Not an Isovar audit bundle")
    if type(manifest.get("schema_version")) is not int or manifest["schema_version"] != SCHEMA_VERSION:
        raise AuditBundleError("Unsupported audit bundle schema version")
    if not isinstance(manifest.get("files"), dict) or not isinstance(manifest.get("metadata"), dict):
        raise AuditBundleError("Bundle files and metadata must be objects")
    if not isinstance(manifest.get("software"), dict):
        raise AuditBundleError("Missing bundle software identity")
    for name, pin in manifest["files"].items():
        _name(name)
        if not isinstance(pin, dict) or type(pin.get("size")) is not int or pin["size"] < 0:
            raise AuditBundleError("Invalid file size: " + name)
        if not isinstance(pin.get("sha256"), str) or not re.fullmatch("[0-9a-f]{64}", pin["sha256"]):
            raise AuditBundleError("Invalid SHA-256: " + name)
    unsigned = {key: value for key, value in manifest.items() if key != "bundle_id"}
    if manifest.get("bundle_id") != identity(unsigned):
        raise AuditBundleError("Bundle manifest identity mismatch")


def publish_json(path, value):
    """Atomically publish JSON, refusing to overwrite existing content."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(canonical(value) + b"\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


class AuditBundle:
    """An immutable manifest and its locally stored audit files.

    Use :meth:`create`, :meth:`open` or :meth:`unpack` to obtain an instance.
    Manifest and file names stay unchanged when the directory is moved.
    """

    def __init__(self, directory, manifest):
        _validate(manifest)
        self.directory = Path(directory).resolve()
        self._manifest = deepcopy(manifest)

    @property
    def manifest(self):
        """Return a detached copy of the versioned, JSON-serializable manifest."""
        return deepcopy(self._manifest)

    @property
    def bundle_id(self):
        """Content identity of the manifest, including every input pin."""
        return self._manifest["bundle_id"]

    @classmethod
    def create(cls, destination, files, metadata=None, *, expected_sha256=None):
        """Copy explicitly selected files and publish a complete manifest.

        Parameters
        ----------
        destination : path-like
            New bundle directory. Existing destinations are never modified.
        files : mapping of str to path-like
            Logical relative names mapped to local source files. No directory
            traversal, network retrieval or implicit whole-library inclusion.
        metadata : dict, optional
            JSON-serializable scientific provenance and structured replay plans.
        expected_sha256 : mapping of str to str, optional
            Expected checksums for logical file names. Supplied pins are checked
            during copying, before a complete bundle can be published.

        Returns
        -------
        AuditBundle
            The published bundle, with its copied files verified.
        """
        if not isinstance(files, Mapping):
            raise AuditBundleError("Files must be a mapping of logical names to local paths")
        files = {_name(name): Path(path) for name, path in files.items()}
        expected_sha256 = {} if expected_sha256 is None else dict(expected_sha256)
        if set(expected_sha256) - set(files):
            raise AuditBundleError("Checksum names must refer to declared input files")
        metadata = {} if metadata is None else deepcopy(metadata)
        if not isinstance(metadata, dict):
            raise AuditBundleError("Metadata must be an object")
        canonical(metadata)
        destination = Path(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.mkdir()  # Reserve ownership before writing or cleaning up.
        try:
            pins = {}
            for name, source in sorted(files.items()):
                if not source.is_file():
                    raise AuditBundleError("Missing input file: " + str(source))
                checksum = digest(source)
                if name in expected_sha256 and checksum != expected_sha256[name]:
                    raise AuditBundleError("Input checksum mismatch: " + name)
                target = destination / "files" / name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, target)
                if digest(target) != checksum:
                    raise AuditBundleError("Input changed during copy: " + str(source))
                pins[name] = dict(sha256=checksum, size=target.stat().st_size)
            manifest = dict(format=FORMAT, schema_version=SCHEMA_VERSION, files=pins,
                            metadata=metadata, software=software_identity())
            manifest["bundle_id"] = identity(manifest)
            publish_json(destination / MANIFEST, manifest)
            return cls.open(destination)
        except BaseException:
            shutil.rmtree(destination)
            raise

    @classmethod
    def open(cls, directory, verify=True):
        """Open a local bundle; verify all payload files by default."""
        directory = Path(directory)
        try:
            manifest = json.loads((directory / MANIFEST).read_text())
        except (OSError, ValueError) as error:
            raise AuditBundleError("Cannot read bundle manifest: " + str(directory)) from error
        bundle = cls(directory, manifest)
        if verify:
            bundle.verify()
        return bundle

    def path(self, name):
        """Resolve a declared logical file name within this bundle."""
        _name(name)
        if name not in self._manifest["files"]:
            raise AuditBundleError("File is not declared in the bundle: " + name)
        path = self.directory / "files" / name
        if not path.resolve().is_relative_to(self.directory / "files"):
            raise AuditBundleError("Bundle file escapes its payload directory: " + name)
        return path

    def verify(self):
        """Verify manifest identity and every file's size and SHA-256.

        Returns
        -------
        dict
            Bundle identity, number of verified files and total payload bytes.
        """
        _validate(self._manifest)
        try:
            on_disk = json.loads((self.directory / MANIFEST).read_text())
        except (OSError, ValueError) as error:
            raise AuditBundleError("Missing or invalid on-disk bundle manifest") from error
        if on_disk != self._manifest:
            raise AuditBundleError("On-disk bundle manifest changed after opening")
        total = 0
        for name, pin in self._manifest["files"].items():
            path = self.path(name)
            if not path.is_file() or path.stat().st_size != pin["size"] or digest(path) != pin["sha256"]:
                raise AuditBundleError("Missing or changed bundle file: " + name)
            total += pin["size"]
        return dict(bundle_id=self.bundle_id, files=len(self._manifest["files"]), bytes=total)

    def pack(self, destination):
        """Write a verified tar.gz backup, refusing an existing archive."""
        self.verify()
        destination = Path(destination)
        if destination.resolve().is_relative_to(self.directory):
            raise AuditBundleError("Archive output must be outside the immutable bundle")
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as handle:
                temporary = Path(handle.name)
            with tarfile.open(temporary, "w:gz", dereference=True) as archive:
                archive.add(self.directory / MANIFEST, arcname=MANIFEST, recursive=False)
                for name in sorted(self._manifest["files"]):
                    archive.add(self.path(name), arcname="files/" + name, recursive=False)
            # Verify the source again: changing inputs must not produce a backup
            # advertised as complete. Restoring also verifies archived bytes.
            self.verify()
            os.link(temporary, destination)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
        return destination

    @classmethod
    def unpack(cls, archive, destination):
        """Restore and verify a tar.gz backup into a new directory."""
        destination = Path(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.mkdir()
        try:
            with tarfile.open(archive, "r:gz") as source:
                members = source.getmembers()
                names = [m.name for m in members]
                if len(names) != len(set(names)) or any(not m.isfile() for m in members):
                    raise AuditBundleError("Archive has duplicate or nonregular entries")
                if MANIFEST not in names or source.getmember(MANIFEST).size > 16 * 1024 * 1024:
                    raise AuditBundleError("Missing or oversized archive manifest")
                manifest = json.load(source.extractfile(MANIFEST))
                _validate(manifest)
                expected = {MANIFEST, *("files/" + name for name in manifest["files"])}
                if set(names) != expected:
                    raise AuditBundleError("Archive does not contain exactly its declared files")
                for name, pin in manifest["files"].items():
                    member = source.getmember("files/" + name)
                    if member.size != pin["size"]:
                        raise AuditBundleError("Archive file size mismatch: " + name)
                    path = destination / "files" / name
                    path.parent.mkdir(parents=True, exist_ok=True)
                    with source.extractfile(member) as handle, path.open("xb") as target:
                        shutil.copyfileobj(handle, target, length=1024 * 1024)
                publish_json(destination / MANIFEST, manifest)
            return cls.open(destination)
        except BaseException as error:
            shutil.rmtree(destination)
            if isinstance(error, (tarfile.TarError, json.JSONDecodeError)):
                raise AuditBundleError("Cannot read audit archive: " + str(error)) from error
            raise
