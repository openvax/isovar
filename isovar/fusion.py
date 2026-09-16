"""Translate sequence-resolved fusion transcripts, never symbolic SV alleles.

Coordinates are plain 0-based, half-open offsets. Input sequence and mappings
must be supplied by an RNA assembly/caller; this module does not discover
fusions, fill sequence from reference, or infer phase across unobserved gaps.
"""

from dataclasses import asdict, dataclass, field
from hashlib import sha256
from typing import Optional, Tuple

from .default_parameters import FUSION_PEPTIDE_LENGTHS, MIN_FUSION_FRAGMENTS
from .genetic_code import standard_genetic_code


def _integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError("%s must be a nonnegative integer" % name)


def _dna(sequence):
    if not sequence or set(sequence) - set("ACGT"):
        raise ValueError("Fusion sequences must contain explicit uppercase A/C/G/T bases")


@dataclass(frozen=True)
class FusionBreakpoint:
    """Oriented partner boundary on the named genomic reference."""

    contig: str
    position: int
    strand: str

    def __post_init__(self):
        _integer(self.position, "breakpoint position")
        if not self.contig or self.strand not in ("+", "-"):
            raise ValueError("A breakpoint needs a contig and +/- strand")


@dataclass(frozen=True)
class FusionBlock:
    """One gapless cDNA-to-genome alignment block, in transcript orientation."""

    query_start: int
    query_end: int
    contig: str
    reference_start: int
    reference_end: int
    strand: str

    def __post_init__(self):
        for name in ("query_start", "query_end", "reference_start", "reference_end"):
            _integer(getattr(self, name), name)
        if (self.query_end <= self.query_start or
                self.reference_end - self.reference_start != self.query_end - self.query_start):
            raise ValueError("Fusion blocks must be nonempty and gapless")
        if not self.contig or self.strand not in ("+", "-"):
            raise ValueError("A block needs a contig and +/- strand")


@dataclass(frozen=True)
class FusionTranscript:
    """Supplied RNA sequence with an explicit donor/insert/acceptor boundary.

    ``junction_start:junction_end`` contains untemplated inserted bases, or is
    empty for a direct join. Blocks cover the two retained partner sequences;
    internal alignment gaps may be represented but cannot establish a frame
    unless an exact, collinear annotated donor match can be validated.
    Provenance must identify the sample, caller/assembler, version, parameters,
    input source and assembly/contig identifier; it is retained verbatim.
    """

    event_id: str
    reference_name: str
    sequence: str
    junction_start: int
    junction_end: int
    donor: FusionBreakpoint
    acceptor: FusionBreakpoint
    blocks: Tuple[FusionBlock, ...]
    provenance: dict

    def __post_init__(self):
        _dna(self.sequence)
        for name in ("junction_start", "junction_end"):
            _integer(getattr(self, name), name)
        if not 0 < self.junction_start <= self.junction_end < len(self.sequence):
            raise ValueError("Fusion junction needs observed sequence on both sides")
        if not self.event_id or not self.reference_name:
            raise ValueError("Fusion event and reference identity are required")
        required = {"sample_id", "method", "version", "parameters", "source", "contig_id"}
        if not required.issubset(self.provenance) or any(
                not self.provenance[k] for k in required - {"parameters"}):
            raise ValueError("Fusion caller/assembly provenance is incomplete")
        if not isinstance(self.provenance["parameters"], dict):
            raise ValueError("Provenance parameters must be a dictionary")
        blocks = tuple(self.blocks)
        object.__setattr__(self, "blocks", blocks)
        left, right = [], []
        previous = 0
        for block in blocks:
            if block.query_start < previous or block.query_end > len(self.sequence):
                raise ValueError("Fusion blocks overlap or lie outside the sequence")
            previous = block.query_end
            if block.query_end <= self.junction_start:
                partner, side = self.donor, left
            elif block.query_start >= self.junction_end:
                partner, side = self.acceptor, right
            else:
                raise ValueError("A partner block crosses the fusion junction")
            if (block.contig, block.strand) != (partner.contig, partner.strand):
                raise ValueError("Fusion block disagrees with partner orientation")
            side.append(block)
        if not left or not right:
            raise ValueError("Both fusion partners need aligned sequence")
        if left[-1].query_end != self.junction_start or right[0].query_start != self.junction_end:
            raise ValueError("Fusion mappings do not reach the declared junction")
        donor_edge = left[-1].reference_end if self.donor.strand == "+" else left[-1].reference_start
        acceptor_edge = right[0].reference_start if self.acceptor.strand == "+" else right[0].reference_end
        if donor_edge != self.donor.position or acceptor_edge != self.acceptor.position:
            raise ValueError("Fusion mapping does not agree with the event breakpoints")


@dataclass(frozen=True)
class FusionReference:
    """Versioned reference transcript, with exon/CDS intervals in interbase units.

    Exons are genomic intervals, in ascending genomic order. The cDNA and CDS
    offsets are in transcript orientation. CDS end includes the stop codon.
    A missing CDS remains noncoding/unknown; a longest ORF is never substituted.
    """

    transcript_id: str
    reference_name: str
    annotation: str
    contig: str
    strand: str
    exons: Tuple[Tuple[int, int], ...]
    sequence: str
    cds_start: Optional[int] = None
    cds_end: Optional[int] = None

    def __post_init__(self):
        _dna(self.sequence)
        if not all((self.transcript_id, self.reference_name, self.annotation, self.contig)):
            raise ValueError("Reference transcript/annotation identity is required")
        if self.strand not in ("+", "-"):
            raise ValueError("Reference strand must be +/-")
        exons = tuple(tuple(e) for e in self.exons)
        object.__setattr__(self, "exons", exons)
        end = 0
        for start, stop in exons:
            _integer(start, "exon start")
            _integer(stop, "exon end")
            if start < end or stop <= start:
                raise ValueError("Reference exons must be ordered and nonoverlapping")
            end = stop
        if sum(b - a for a, b in exons) != len(self.sequence):
            raise ValueError("Transcript sequence and exon lengths disagree")
        if (self.cds_start is None) != (self.cds_end is None):
            raise ValueError("Both CDS boundaries or neither must be supplied")
        if self.cds_start is not None:
            _integer(self.cds_start, "CDS start")
            _integer(self.cds_end, "CDS end")
            coding = self.sequence[self.cds_start:self.cds_end]
            if (not 0 <= self.cds_start < self.cds_end <= len(self.sequence) or len(coding) % 3 or
                    coding[:3] not in standard_genetic_code.start_codons or
                    coding[-3:] not in standard_genetic_code.stop_codons or
                    any(coding[i:i + 3] in standard_genetic_code.stop_codons for i in range(0, len(coding) - 3, 3))):
                raise ValueError("Reference CDS must have a justified complete start/frame/stop")

    def spliced_offset(self, position):
        """Map one genomic base (0-based) to a transcript base."""
        offset = 0
        for start, end in (self.exons if self.strand == "+" else self.exons[::-1]):
            if start <= position < end:
                return offset + (position - start if self.strand == "+" else end - 1 - position)
            offset += end - start
        return None


@dataclass(frozen=True)
class FusionRead:
    """An observed contiguous RNA substring and its original alignment mapping.

    Blocks use offsets into this read's ``sequence``. A read can be a clipped
    window, with ``source_query_start`` recording its offset in the original
    transcript-oriented read. Fragment IDs must identify physical templates
    within a sample/library (or validated cell/UMI groups), not processed BAMs.
    Supplementary records for a read must be resolved by the input adapter,
    never submitted as independent reads or joined across alternative mappings.
    """

    sample_id: str
    library_id: str
    fragment_id: str
    read_id: str
    source: str
    source_query_start: int
    cdna_start: int
    sequence: str
    blocks: Tuple[FusionBlock, ...] = field(default_factory=tuple)

    def __post_init__(self):
        _dna(self.sequence)
        for name in ("source_query_start", "cdna_start"):
            _integer(getattr(self, name), name)
        if not all((self.sample_id, self.library_id, self.fragment_id, self.read_id, self.source)):
            raise ValueError("RNA evidence needs sample/library/read/fragment/source identities")
        object.__setattr__(self, "blocks", tuple(self.blocks))


def _coordinates(blocks):
    result = {}
    for block in blocks:
        positions = range(block.reference_start, block.reference_end)
        if block.strand == "-":
            positions = reversed(positions)
        for query, position in zip(range(block.query_start, block.query_end), positions):
            value = (block.contig, position, block.strand)
            if query in result:
                raise ValueError("Overlapping or alternative mappings in one RNA observation")
            result[query] = value
    return result


def _match_reference(fusion, reference, side):
    partner = fusion.donor if side == "donor" else fusion.acceptor
    lo, hi = (0, fusion.junction_start) if side == "donor" else (fusion.junction_end, len(fusion.sequence))
    if (reference.reference_name, reference.contig, reference.strand) != (
            fusion.reference_name, partner.contig, partner.strand):
        return None
    coordinates = _coordinates([b for b in fusion.blocks if lo <= b.query_start and b.query_end <= hi])
    offsets = [reference.spliced_offset(coordinates[q][1]) if q in coordinates else None for q in range(lo, hi)]
    if not offsets or offsets[0] is None or offsets != list(range(offsets[0], offsets[0] + hi - lo)):
        return None
    start = offsets[0]
    if fusion.sequence[lo:hi] != reference.sequence[start:start + hi - lo]:
        return None
    return start


def reconstruct_fusion(fusion, references=(), reads=(), peptide_lengths=FUSION_PEPTIDE_LENGTHS,
                       min_fragments=MIN_FUSION_FRAGMENTS):
    """Validate supplied fusion RNA and return a JSON-serializable evidence result.

    Frames are transferred only from exact, collinear annotated donor matches.
    All compatible reference models are retained. No reference bases are added;
    novel acceptor sequence is translated in the donor frame, not its own ORF.
    A partial 5' CDS is explicitly conditional on the annotated upstream frame.
    Missing support, missing frame and alternative frames remain distinct.

    Returns
    -------
    dict
        RNA sequence, junction/mapping/provenance, evidence counts and all
        supported translation hypotheses. Only ``status == 'translated'`` has
        one resolved translation; ``ambiguous`` must not be silently ranked.
    """
    _integer(min_fragments, "min_fragments")
    if min_fragments < 1:
        raise ValueError("At least one direct supporting fragment is required")
    lengths = sorted(set(peptide_lengths))
    for length in lengths:
        _integer(length, "peptide length")
        if length < 1:
            raise ValueError("Peptide lengths must be positive")
    if not lengths:
        raise ValueError("At least one peptide length is required")
    references = tuple(references)
    ids = [(r.annotation, r.transcript_id) for r in references]
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate reference transcript identity")
    coordinates = _coordinates(fusion.blocks)
    observations, fragments, direct, seen = [], set(), set(), {}
    for read in reads:
        if read.sample_id != fusion.provenance["sample_id"]:
            raise ValueError("RNA evidence belongs to another sample")
        start, end = read.cdna_start, read.cdna_start + len(read.sequence)
        if end > len(fusion.sequence) or fusion.sequence[start:end] != read.sequence:
            raise ValueError("RNA observation does not match the supplied fusion sequence")
        mapping = _coordinates(read.blocks)
        if any(q >= len(read.sequence) or coordinates.get(q + start) != p for q, p in mapping.items()):
            raise ValueError("RNA observation has a conflicting partner mapping")
        key = (read.sample_id, read.library_id, read.read_id)
        # Copies from another processed file do not add evidence. Disagreeing
        # placements of one read are a conflict, not a larger inferred witness.
        signature = (read.fragment_id, start, read.source_query_start, read.sequence, tuple(read.blocks))
        if key in seen:
            if seen[key] != signature:
                raise ValueError("Conflicting observations for the same RNA read")
            continue
        seen[key] = signature
        fragment = (read.sample_id, read.library_id, read.fragment_id)
        fragments.add(fragment)
        spans = (start < fusion.junction_start <= fusion.junction_end < end and
                 fusion.junction_start - 1 - start in mapping and fusion.junction_end - start in mapping)
        if spans:
            direct.add(fragment)
        observations.append(dict(asdict(read), directly_spans_junction=spans))
    matches = {side: [(ref, offset) for ref in references
                     if (offset := _match_reference(fusion, ref, side)) is not None]
               for side in ("donor", "acceptor")}
    result = dict(schema_version=1, event_id=fusion.event_id, reference_name=fusion.reference_name,
                  cdna_sequence=fusion.sequence, sequence_sha256=sha256(fusion.sequence.encode()).hexdigest(),
                  junction_interval=[fusion.junction_start, fusion.junction_end],
                  donor=asdict(fusion.donor), acceptor=asdict(fusion.acceptor), blocks=[asdict(b) for b in fusion.blocks],
                  provenance=fusion.provenance, status="unresolved_frame", reasons=[], translations=[],
                  compatible_transcripts={side: [r.transcript_id for r, _ in rows] for side, rows in matches.items()},
                  reference_annotations=sorted({r.annotation for r in references}),
                  evidence=dict(reads=len(observations), fragments=len(fragments),
                                directly_spanning_fragments=len(direct), observations=observations),
                  parameters=dict(peptide_lengths=lengths, min_fragments=min_fragments))
    if len(direct) < min_fragments:
        result.update(status="insufficient_support", reasons=["insufficient_direct_junction_fragments"])
        return result
    groups = {}
    for reference, offset in matches["donor"]:
        if reference.cds_start is None:
            result["reasons"].append("donor_CDS_unavailable:" + reference.transcript_id)
            continue
        projected_start = reference.cds_start - offset
        if projected_start >= fusion.junction_start:
            result["reasons"].append("junction_before_donor_CDS:" + reference.transcript_id)
            continue
        if offset + fusion.junction_start > reference.cds_end - 3:
            result["reasons"].append("junction_after_donor_CDS:" + reference.transcript_id)
            continue
        start = projected_start if projected_start >= 0 else projected_start % 3
        protein, stop = standard_genetic_code.translate(fusion.sequence[start:], first_codon_is_start=projected_start >= 0)
        coding_end = start + len(protein) * 3
        if coding_end <= fusion.junction_end:
            result["reasons"].append("no_translated_acceptor_context:" + reference.transcript_id)
            continue
        key = (start, protein, stop, projected_start >= 0)
        if key not in groups:
            junction = [fusion.junction_start - start, fusion.junction_end - start]
            peptides = []
            for length in lengths:
                for i in range(len(protein) - length + 1):
                    boundaries = [b for b in sorted(set(junction)) if 3 * i < b < 3 * (i + length)]
                    if boundaries:
                        peptides.append(dict(sequence=protein[i:i + length], protein_interval=[i, i + length],
                                             junction_boundaries_in_cds=boundaries))
            groups[key] = dict(amino_acids=protein, translation_start=start,
                               cds_start=(start if projected_start >= 0 else None),
                               complete_5prime=projected_start >= 0, ends_with_stop_codon=stop,
                               trailing_partial_codon_bases=(0 if stop else (len(fusion.sequence) - start) % 3),
                               junction_in_translated_cds=junction, junction_peptides=peptides,
                               donor_transcript_ids=[], frame_evidence=[], acceptor_frames=[])
            for acceptor, acceptor_offset in matches["acceptor"]:
                if acceptor.cds_start is not None and acceptor.cds_start <= acceptor_offset < acceptor.cds_end - 3:
                    in_frame = (fusion.junction_end - start) % 3 == (acceptor_offset - acceptor.cds_start) % 3
                    groups[key]["acceptor_frames"].append(dict(transcript_id=acceptor.transcript_id, in_frame=in_frame))
        groups[key]["donor_transcript_ids"].append(reference.transcript_id)
        groups[key]["frame_evidence"].append(dict(transcript_id=reference.transcript_id, annotation=reference.annotation,
            donor_reference_offset=offset, reference_cds_start=reference.cds_start,
            basis="exact_collinear_donor_match", upstream_frame_assumed=projected_start < 0))
    result["translations"] = list(groups.values())
    if groups:
        # Noncoding alternatives are also real uncertainty, not permission to
        # pick the one donor model that happens to produce a peptide.
        result["status"] = "ambiguous" if len(groups) > 1 or result["reasons"] else "translated"
    elif not matches["donor"]:
        result["reasons"].append("no_exact_collinear_annotated_donor")
    return result


def fusion_from_dict(data):
    """Decode the explicit JSON input used by ``isovar fusion``."""
    transcript = dict(data["fusion"])
    transcript["donor"] = FusionBreakpoint(**transcript["donor"])
    transcript["acceptor"] = FusionBreakpoint(**transcript["acceptor"])
    transcript["blocks"] = tuple(FusionBlock(**b) for b in transcript["blocks"])
    reads = []
    for entry in data.get("reads", []):
        entry = dict(entry)
        entry["blocks"] = tuple(FusionBlock(**b) for b in entry.get("blocks", []))
        reads.append(FusionRead(**entry))
    return FusionTranscript(**transcript), tuple(FusionReference(**r) for r in data.get("references", [])), tuple(reads)
