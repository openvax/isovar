"""Interpret nominated SV source rows without silently repairing coordinates.

Varcode owns breakend geometry. This layer preserves source provenance, separates
PCR validation from caller support and applies only explicitly selected, pinned
source conventions. A resolved adjacency is a nomination, not RNA support.
"""

from collections import Counter
from dataclasses import asdict, dataclass
import gzip
from importlib.metadata import version
from pathlib import Path
import re
from typing import Any, Literal

from pyensembl import Genome
from varcode import ContigResolver
from varcode.structural_variant import StructuralVariant
from varcode.sv_allele_parser import breakend_mate, parse_symbolic_alt

from ._file_identity import file_sha256
from . import __version__


SourceProfile = Literal['standard', 'colo829_paulin_paz_v2', 'hcc1395_svscope_2026']
SourceState = Literal['resolved_adjacency', 'unresolved', 'invalid']
SOURCE_PROFILES: dict[SourceProfile, str | None] = {
    'standard': None,
    'colo829_paulin_paz_v2': '8ad075bf43aa71d89c591396569e813959ada44a30f171982c2d299411290bac',
    'hcc1395_svscope_2026': '892222606e3951f20a8608b28ff32e6dc320257d573c9b832108156e9a361424',
}
_CALLER_COORDINATES = re.compile(
    r'^(?:ONT|PB)_TDScope\.(?P<kind>DEL|DUP|INV|BND)\.'
    r'(?P<contig>[^:]+):(?P<position>\d+)-(?P<mate_contig>[^:]+):(?P<mate_position>\d+)(?:_\d+)?$')


def _info(column: str) -> dict[str, str | bool]:
    """Preserve INFO values literally; reject ambiguous repeated keys."""
    fields: dict[str, str | bool] = {}
    for entry in () if column == '.' else column.split(';'):
        key, separator, value = entry.partition('=')
        if not key or key in fields:
            raise ValueError('Empty or repeated INFO key: ' + key)
        fields[key] = value if separator else True
    return fields


@dataclass(frozen=True)
class CallerCoordinateClaim:
    """Coordinates encoded in a caller ID, with unverified coordinate convention.

    These are source claims only. They are never promoted to resolved endpoints
    or used to infer missing retained sides, sequence or RNA support.
    """

    evidence_id: str
    kind: str
    contig: str
    position: int
    mate_contig: str
    mate_position: int


@dataclass(frozen=True)
class SVSourceRecord:
    """One original row and its explicit, source-specific interpretation."""

    line_number: int
    original_row: str
    record_id: str
    state: SourceState
    reasons: tuple[str, ...]
    interpreted_alt: str | None
    transformations: tuple[str, ...]
    pcr_validation: bool | None
    caller_evidence: tuple[str, ...]
    coordinate_claims: tuple[CallerCoordinateClaim, ...]
    variant: StructuralVariant | None

    def for_annotation(self, genome: Genome, resolver: ContigResolver) -> StructuralVariant:
        """Resolve both source contigs against one explicit annotation vocabulary.

        Reuse a resolver built from ``genome.contigs()`` for the whole inventory.
        Missing or ambiguous contigs are errors, not evidence of an intergenic
        event. The original source variant and reported coordinates stay intact.
        """
        if self.variant is None:
            raise ValueError('Unresolved source records cannot be annotated')
        if (genome.reference_name != self.variant.reference_name
                or resolver.reference_name != genome.reference_name):
            raise ValueError('SV source annotation assembly mismatch')
        fields = self.variant.to_dict()
        for key in ('contig', 'mate_contig'):
            if fields[key] is None:
                continue
            resolution = resolver.resolve(fields[key])
            if resolution.resolved_name is None:
                raise ValueError('Unresolved annotation contig: ' + resolution.diagnostic)
            fields[key] = resolution.resolved_name
        return StructuralVariant(**dict(fields, genome=genome, normalize_contig_names=False,
                                         convert_ucsc_contig_names=False))

    def to_dict(self) -> dict[str, Any]:
        """Return a detached JSON record, including oriented source coordinates."""
        return dict(line_number=self.line_number, original_row=self.original_row,
            record_id=self.record_id, state=self.state, reasons=list(self.reasons),
            interpreted_alt=self.interpreted_alt, transformations=list(self.transformations),
            pcr_validation=self.pcr_validation, caller_evidence=list(self.caller_evidence),
            coordinate_claims=[dict(asdict(claim), coordinate_convention='source_id_unverified')
                               for claim in self.coordinate_claims],
            confidence_offsets={name: list(offsets) if offsets is not None else None
                                for name, offsets in (('local', self.variant.ci_start), ('mate', self.variant.ci_end))}
                               if self.variant else None,
            junctions=[[dict(contig=end.contig, position_1based=end.position, keeps=end.keeps)
                        for end in pair] for pair in self.variant.junctions] if self.variant else [])


def _metadata(info: dict[str, str | bool]) -> tuple[bool | None, tuple[str, ...], tuple[CallerCoordinateClaim, ...]]:
    pcr = info.get('PCR_Validation')
    if pcr not in (None, 'True', 'False'):
        raise ValueError('PCR_Validation must be literal True or False when supplied')
    support = info.get('SUPP_VEC', info.get('SUPPORT', ''))
    if not isinstance(support, str):
        raise ValueError('Caller support must retain its source string')
    evidence = tuple(support.split(',')) if support else ()
    claims = []
    for identifier in evidence:
        if match := _CALLER_COORDINATES.fullmatch(identifier):
            fields = match.groupdict()
            claims.append(CallerCoordinateClaim(identifier, fields['kind'], fields['contig'],
                int(fields['position']), fields['mate_contig'], int(fields['mate_position'])))
    return None if pcr is None else pcr == 'True', evidence, tuple(claims)


def _interpret_row(line: str, number: int, profile: SourceProfile, reference_name: str) -> SVSourceRecord:
    columns = line.split('\t')
    record_id = columns[2] if len(columns) > 2 else '.'
    pcr: bool | None = None
    evidence: tuple[str, ...] = ()
    claims: tuple[CallerCoordinateClaim, ...] = ()
    interpreted_alt: str | None = None
    transformations: tuple[str, ...] = ()
    variant: StructuralVariant | None = None
    state: SourceState = 'unresolved'
    reasons: tuple[str, ...] = ()
    try:
        if len(columns) < 8 or not columns[0] or int(columns[1]) < 1:
            raise ValueError('SV source rows require eight fields and a positive VCF POS')
        info = _info(columns[7])
        pcr, evidence, claims = _metadata(info)
        interpreted_alt = columns[4]
        if profile == 'colo829_paulin_paz_v2' and interpreted_alt == '<BND>':
            embedded = info.get('END')
            mate = breakend_mate(embedded) if isinstance(embedded, str) else None
            if not isinstance(embedded, str) or mate is None:
                raise ValueError('Pinned COLO829 BND lacks the documented END breakend encoding')
            partner = mate[0]
            if info.get('CHR2') != partner or info.get('SVTYPE') != 'BND':
                raise ValueError('COLO829 embedded breakend conflicts with CHR2 or SVTYPE')
            interpreted_alt = embedded
            info.pop('END')
            transformations = ('colo829_v2_literal_breakend_from_END',)
        if profile == 'hcc1395_svscope_2026' and set(columns[3] + interpreted_alt) <= set('ACGTN'):
            # The author's writer uses a normal-sample consensus as REF and a
            # discovery-window start as POS. It is not an ordinary anchored VCF
            # allele; moving POS or replacing REF would invent a new definition.
            reasons = ('source_consensus_alleles_require_reference_anchoring',)
        else:
            variant = parse_symbolic_alt(columns[0], int(columns[1]), columns[3], interpreted_alt,
                                          info, genome=reference_name, normalize_contig_names=False,
                                          convert_ucsc_contig_names=False)
            if variant is None:
                reasons = ('not_a_resolved_structural_allele',)
            elif not variant.junctions or any(end.keeps not in ('left', 'right')
                                               for pair in variant.junctions for end in pair):
                reasons = ('oriented_adjacency_or_inserted_sequence_unresolved',)
                variant = None
            else:
                state = 'resolved_adjacency'
    except (TypeError, ValueError) as error:
        state, reasons, variant = 'invalid', (str(error),), None
    return SVSourceRecord(number, line, record_id, state, reasons, interpreted_alt,
                          transformations, pcr, evidence, claims, variant)


def read_sv_source(path: str | Path, *, sha256: str, reference_name: str,
                   profile: SourceProfile = 'standard') -> tuple[SVSourceRecord, ...]:
    """Validate a source digest and interpret every original data row.

    Parameters
    ----------
    path : str or pathlib.Path
        Unmodified source VCF, optionally gzip compressed.
    sha256 : str
        Required digest of the source file bytes, including compression.
    reference_name : str
        Declared assembly. Named source profiles require GRCh38.
    profile : str, optional
        Explicit source convention; defaults to strict standard interpretation.
        Named profiles accept only their documented original file checksum.

    Returns
    -------
    tuple of SVSourceRecord
        Every original row, including invalid and unresolved definitions. PCR
        absence is unknown, not a negative validation result. Input genotypes
        and caller support remain in the original row; no event is RNA-supported
        merely because its source adjacency was interpretable.
        Use ``record.for_annotation(genome, resolver)`` to resolve names against
        an annotation vocabulary before effect prediction.
    """
    path = Path(path)
    if profile not in SOURCE_PROFILES:
        raise ValueError('Unknown SV source profile')
    if not reference_name:
        raise ValueError('SV source requires an explicit reference assembly')
    if len(sha256) != 64 or file_sha256(path) != sha256:
        raise ValueError('SV source checksum mismatch')
    pinned = SOURCE_PROFILES[profile]
    if pinned is not None and (sha256 != pinned or reference_name != 'GRCh38'):
        raise ValueError('Named SV source profiles require their pinned GRCh38 original')
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as stream:
        records = tuple(_interpret_row(line.rstrip('\r\n'), number, profile, reference_name)
                        for number, line in enumerate(stream, 1) if line.strip() and not line.startswith('#'))
    if file_sha256(path) != sha256:
        raise ValueError('SV source changed during interpretation')
    return records


def sv_source_report(records: tuple[SVSourceRecord, ...], *, sha256: str, reference_name: str,
                     profile: SourceProfile) -> dict[str, Any]:
    """Summarize row interpretation without treating it as RNA or protein evidence."""
    return dict(schema='isovar.sv_source.v1', source_sha256=sha256, reference_name=reference_name,
                producer_version=__version__, varcode_version=version('varcode'),
                interpretation_sha256=file_sha256(__file__),
                profile=profile, counts=dict(Counter(record.state for record in records)),
                pcr_validated_records=sum(record.pcr_validation is True for record in records),
                records=[record.to_dict() for record in records])
