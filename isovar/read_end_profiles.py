"""Explicit, versioned library profiles for shared technical-sequence inference."""

from dataclasses import asdict, dataclass
from functools import cached_property
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Literal, Mapping

import pysam

from ._file_identity import file_sha256
from .default_parameters import (
    READ_END_WINDOW, MIN_ADAPTER_OVERLAP, MAX_ADAPTER_ERROR_RATE,
    MIN_POLY_A_LENGTH, MAX_POLY_A_ERROR_RATE,
)


TechnicalRole = Literal['adapter', 'primer', 'linker']


def _text(value: object, name: str, optional: bool = False) -> None:
    if optional and value is None:
        return
    if not isinstance(value, str) or not value.strip():
        raise ValueError(name + ' must be a nonempty string')


@dataclass(frozen=True)
class ReadProcessingStep:
    """Documented upstream processing, without inferring a trimming outcome."""

    tool: str
    operation: str
    version: str | None = None
    source: str | None = None

    def __post_init__(self) -> None:
        _text(self.tool, 'Processing tool')
        _text(self.operation, 'Processing operation')
        _text(self.version, 'Processing version', optional=True)
        _text(self.source, 'Processing source', optional=True)


@dataclass(frozen=True)
class ReadLibraryProvenance:
    """Supplied library facts; unknown values stay None, never guessed from PL.

    Sources may identify a sample sheet, acquisition receipt or manufacturer
    documentation. Processing records describe what was run; absence of a visible
    adapter or tail does not establish biological absence or read completeness.
    """

    platform: str | None = None
    instrument: str | None = None
    preparation: str | None = None
    kit: str | None = None
    kit_revision: str | None = None
    library_state: str | None = None
    processing: tuple[ReadProcessingStep, ...] = ()
    sources: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        for name in ('platform', 'instrument', 'preparation', 'kit', 'kit_revision', 'library_state'):
            _text(getattr(self, name), name, optional=True)
        if isinstance(self.processing, (str, bytes)) or isinstance(self.sources, (str, bytes)):
            raise ValueError('Processing and sources must be sequences, not a single string')
        object.__setattr__(self, 'processing', tuple(self.processing))
        object.__setattr__(self, 'sources', tuple(self.sources))
        if any(not isinstance(step, ReadProcessingStep) for step in self.processing):
            raise ValueError('Library processing must contain ReadProcessingStep values')
        for source in self.sources:
            _text(source, 'Library provenance source')

    @classmethod
    def from_dict(cls, values: Mapping[str, Any]) -> 'ReadLibraryProvenance':
        fields = dict(values)
        fields['processing'] = tuple(ReadProcessingStep(**step) for step in fields.get('processing', ()))
        return cls(**fields)


@dataclass(frozen=True)
class Adapter:
    """Explicit technical motif in sequencing orientation, not BAM orientation.

    ``end`` is left or right; terminal partial suffixes/prefixes respectively
    are allowed. ``mate`` is 0 (either/unpaired), 1 or 2. Sequence must be literal
    A/C/G/T; variable barcode positions are not informative adapter sequence.
    ``role`` distinguishes adapter, primer and linker candidates without changing
    matching thresholds. ``source`` records the actual sequence provenance.
    """

    name: str
    sequence: str
    end: Literal['left', 'right'] = 'right'
    mate: int = 0
    min_overlap: int = MIN_ADAPTER_OVERLAP
    max_error_rate: float = MAX_ADAPTER_ERROR_RATE
    role: TechnicalRole = 'adapter'
    source: str | None = None

    def __post_init__(self) -> None:
        _text(self.name, 'Adapter name')
        _text(self.sequence, 'Adapter sequence')
        _text(self.source, 'Adapter source', optional=True)
        object.__setattr__(self, 'sequence', self.sequence.upper())
        if (self.end not in ('left', 'right') or self.role not in ('adapter', 'primer', 'linker')
                or type(self.mate) is not int or self.mate not in (0, 1, 2)
                or type(self.min_overlap) is not int
                or not 4 <= self.min_overlap <= len(self.sequence)
                or set(self.sequence) - set('ACGT')
                or not isinstance(self.max_error_rate, (int, float))
                or not 0 <= self.max_error_rate < 0.5):
            raise ValueError('Invalid adapter sequence, end, mate, role or matching threshold')


@dataclass(frozen=True)
class ReadEndProfile:
    """Versioned library configuration; unknown chemistry needs no guessed kit.

    Only explicit sequences are searched. ``both_orientations`` supports either
    cDNA orientation; it is not inferred from an instrument. Provenance is supplied
    evidence, not a classifier. Poly-A/T candidates do not establish full length.
    """

    name: str = 'unknown'
    version: str = '1'
    adapters: tuple[Adapter, ...] = ()
    both_orientations: bool = False
    end_window: int = READ_END_WINDOW
    min_poly_a_length: int = MIN_POLY_A_LENGTH
    max_poly_a_error_rate: float = MAX_POLY_A_ERROR_RATE
    provenance: ReadLibraryProvenance = ReadLibraryProvenance()

    def __post_init__(self) -> None:
        _text(self.name, 'Profile name')
        _text(self.version, 'Profile version')
        object.__setattr__(self, 'adapters', tuple(self.adapters))
        if (type(self.min_poly_a_length) is not int or type(self.end_window) is not int
                or type(self.both_orientations) is not bool
                or not 4 <= self.min_poly_a_length <= self.end_window
                or not isinstance(self.max_poly_a_error_rate, (int, float))
                or not 0 <= self.max_poly_a_error_rate < 0.5
                or not isinstance(self.provenance, ReadLibraryProvenance)
                or any(not isinstance(adapter, Adapter) for adapter in self.adapters)
                or any(len(adapter.sequence) > self.end_window for adapter in self.adapters)):
            raise ValueError('Invalid read-end profile, provenance or matching thresholds')

    @classmethod
    def from_dict(cls, values: Mapping[str, Any]) -> 'ReadEndProfile':
        fields = dict(values)
        fields['adapters'] = tuple(Adapter(**adapter) for adapter in fields.get('adapters', ()))
        if 'provenance' in fields:
            fields['provenance'] = ReadLibraryProvenance.from_dict(fields['provenance'])
        return cls(**fields)

    @cached_property
    def fingerprint(self) -> str:
        """Identify the full configuration, including supplied sequence provenance."""
        return sha256(json.dumps(asdict(self), sort_keys=True).encode()).hexdigest()


UNKNOWN_PROFILE = ReadEndProfile()


def read_end_profile_from_fasta(filename: str | Path, *, sha256: str,
                                profile: Mapping[str, Any], adapter_defaults: Mapping[str, Any] | None = None
                                ) -> ReadEndProfile:
    """Build a profile from a checksum-pinned FASTA and explicit shared settings.

    ``profile`` supplies the ordinary name/version, orientation, provenance and
    tail settings. ``adapter_defaults`` supplies end, mate, role and matching
    thresholds for every FASTA entry. Names and sequences come from FASTA; no kit,
    mate or orientation is guessed from a header. Empty, duplicate-ID and FASTQ
    inputs are rejected. Use JSON adapters when individual entries differ.
    """
    filename = Path(filename).expanduser()
    if file_sha256(filename) != sha256:
        raise ValueError('Adapter FASTA checksum mismatch')
    if 'adapters' in profile:
        raise ValueError('FASTA profile cannot also supply inline adapters')
    defaults = dict(adapter_defaults or {})
    if {'name', 'sequence'} & defaults.keys():
        raise ValueError('Adapter FASTA supplies names and sequences')
    defaults.setdefault('source', f'sha256:{sha256}')
    adapters = []
    names = set()
    with pysam.FastxFile(str(filename)) as records:
        for record in records:
            if record.quality is not None or record.name in names:
                raise ValueError('Adapter FASTA requires unique IDs and no FASTQ qualities')
            names.add(record.name)
            adapters.append(dict(name=record.name, sequence=record.sequence, **defaults))
    if not adapters:
        raise ValueError('Adapter FASTA is empty')
    if file_sha256(filename) != sha256:
        raise ValueError('Adapter FASTA changed during reading')
    return ReadEndProfile.from_dict(dict(profile, adapters=adapters))


def _profile_from_json(values: Mapping[str, Any], directory: Path | None) -> ReadEndProfile:
    fields = dict(values)
    fasta = fields.pop('adapter_fasta', None)
    if fasta is None:
        return ReadEndProfile.from_dict(fields)
    if not isinstance(fasta, Mapping) or set(fasta) - {'path', 'sha256', 'defaults'}:
        raise ValueError('adapter_fasta accepts path, sha256 and defaults only')
    if not {'path', 'sha256'} <= fasta.keys():
        raise ValueError('adapter_fasta requires path and sha256')
    if directory is None:
        raise ValueError('Inline profiles cannot reference FASTA paths; use a pinned profile file_id')
    return read_end_profile_from_fasta(directory / fasta['path'], sha256=fasta['sha256'],
                                       profile=fields, adapter_defaults=fasta.get('defaults'))


def read_end_profiles_from_dict(values: Mapping[str, Any], *, directory: Path | None = None
                                ) -> ReadEndProfile | dict[str, ReadEndProfile]:
    """Parse inline profile JSON; external FASTA requires an explicit base directory."""
    if 'read_groups' in values:
        if set(values) != {'read_groups'}:
            raise ValueError('read_groups cannot be combined with a global profile')
        return {group: _profile_from_json(profile, directory)
                for group, profile in values['read_groups'].items()}
    return _profile_from_json(values, directory)


def read_end_profiles_from_json(filename: str | Path) -> ReadEndProfile | dict[str, ReadEndProfile]:
    """Load one explicit profile or a per-RG map; FASTA paths are JSON-relative.

    Missing read groups use the unknown profile with no cross-library fallback.
    An empty-string group identifies records without an RG tag.
    """
    filename = Path(filename).expanduser()
    return read_end_profiles_from_dict(json.loads(filename.read_text()), directory=filename.parent)


def read_end_profile_file_id(value: Any) -> str | None:
    """Identify an explicit manifest profile reference, rejecting mixed declarations."""
    if not isinstance(value, Mapping) or 'file_id' not in value:
        return None
    if set(value) != {'file_id'}:
        raise ValueError('A profile file_id cannot be combined with inline settings')
    _text(value['file_id'], 'Profile file ID')
    return value['file_id']


def resolve_read_end_profile(value: Any, files: Mapping[str, Path]) -> ReadEndProfile | dict[str, ReadEndProfile] | None:
    """Use the same profiles for typed callers, inline JSON and pinned file IDs.

    The caller verifies manifest file digests; referenced FASTA files are checked
    by the profile loader. No directory or nearest-kit fallback is inferred.
    """
    if value is None or isinstance(value, ReadEndProfile):
        return value
    identifier = read_end_profile_file_id(value)
    if identifier is not None:
        if identifier not in files:
            raise ValueError('Read-end profile refers to undeclared file ID: ' + identifier)
        return read_end_profiles_from_json(files[identifier])
    if not isinstance(value, Mapping):
        raise ValueError('Read-end profile must be a typed profile, inline mapping or file_id')
    if all(isinstance(profile, ReadEndProfile) for profile in value.values()):
        return dict(value)
    return read_end_profiles_from_dict(value)
