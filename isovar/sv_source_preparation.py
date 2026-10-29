"""Prepare validated SV nominations for the existing RNA audit workflow."""

from dataclasses import asdict, dataclass
from typing import Any

from pyensembl import Genome
from varcode import ContigResolver
from varcode.structural_variant import Breakend

from .fusion import FusionBreakpoint, FusionReference
from .sv_rna_reconciliation import hypotheses_from_varcode
from .sv_sources import SVSourceRecord


@dataclass(frozen=True)
class SVSourceView:
    """One oriented nomination, in the original source contig vocabulary."""

    view_id: str
    donor: FusionBreakpoint
    acceptor: FusionBreakpoint
    regions: tuple[tuple[str, int, int], ...]

    def to_dict(self, reference_name: str) -> dict[str, Any]:
        """Return the audit's existing input schema, without RNA evidence."""
        return dict(reference_name=reference_name, donor=asdict(self.donor),
                    acceptor=asdict(self.acceptor), regions=[list(region) for region in self.regions])


@dataclass(frozen=True)
class PreparedSVSourceRecord:
    """Reference hypotheses and nomination views; never observed RNA or protein."""

    reference_name: str
    annotation: str
    references: tuple[FusionReference, ...]
    hypotheses: tuple[dict[str, Any], ...]
    views: tuple[SVSourceView, ...]
    limitations: tuple[str, ...]
    gene_names: tuple[str, ...]
    nominal_igv: tuple[str, ...]
    annotation_flank: int

    def reference_payload(self) -> dict[str, Any]:
        """Serialize reference models for an ordinary SV audit manifest."""
        return dict(models=[asdict(model) for model in self.references],
                    limitations=list(self.limitations), annotation=self.annotation,
                    reference_name=self.reference_name,
                    selection=f'all transcripts of genes overlapping a {self.annotation_flank}bp breakpoint window')


def _breakpoint(end: Breakend, donor: bool) -> FusionBreakpoint:
    """Convert a retained 1-based base to the corresponding interbase boundary."""
    left = end.keeps == 'left'
    return FusionBreakpoint(end.contig, end.position if left else end.position - 1,
                            ('+' if left else '-') if donor else ('-' if left else '+'))


def _reference_model(transcript: Any, source_contig: str, reference_name: str,
                     annotation: str) -> tuple[FusionReference | None, str | None]:
    sequence = transcript.sequence
    if not sequence or set(sequence.upper()) - set('ACGT'):
        return None, transcript.id + ': missing_or_ambiguous_cdna'
    fields = dict(transcript_id=transcript.id, reference_name=reference_name, annotation=annotation,
                  contig=source_contig, strand=transcript.strand,
                  exons=sorted((exon.start - 1, exon.end) for exon in transcript.exons), sequence=sequence.upper())
    complete = (transcript.contains_start_codon and transcript.contains_stop_codon
                and transcript.start_codon_complete and transcript.stop_codon_complete)
    if complete:
        try:
            return FusionReference(**fields, cds_start=transcript.first_start_codon_spliced_offset,
                                   cds_end=transcript.last_stop_codon_spliced_offset + 1), None
        except ValueError as error:
            return FusionReference(**fields), transcript.id + ': CDS validation: ' + str(error)
    return FusionReference(**fields), None


def prepare_sv_source_record(record: SVSourceRecord, genome: Genome,
                             resolver: ContigResolver, *, annotation: str,
                             annotation_flank: int = 10000, region_flank: int = 500) -> PreparedSVSourceRecord:
    """Prepare one resolved source through the shared reference/hypothesis flow.

    Parameters
    ----------
    record : SVSourceRecord
        A resolved adjacency from ``read_sv_source``. Unresolved inputs raise
        ValueError; retain them explicitly as unassessed manifest events.
    genome : pyensembl.Genome
        Explicit annotation context; sequence downloads are the caller's concern.
    resolver : varcode.ContigResolver
        Reuse one resolver built from this genome's contigs across the inventory.
    annotation : str
        Reproducible annotation label, for example ``Ensembl 87``.
    annotation_flank, region_flank : int
        Nonnegative bases on either side of the nominated retained base, for
        transcript selection and regional RNA acquisition, respectively.

    Returns
    -------
    PreparedSVSourceRecord
        Both orientations of every adjacency, ordinary reference models and all
        Varcode hypotheses. Missing annotation contigs remain explicit limitations
        and do not remove RNA nominations on custom, viral or decoy contigs.
        Source aliases are preserved; no liftover, RNA support or translation is
        inferred. Breakpoints remain nominal; confidence intervals are retained
        in the accompanying original source interpretation.
    """
    variant = record.variant
    if variant is None:
        raise ValueError('Preparing an SV source requires a resolved adjacency')
    if (genome.reference_name != variant.reference_name
            or resolver.reference_name != genome.reference_name):
        raise ValueError('SV source preparation assembly mismatch')
    if not annotation:
        raise ValueError('SV source preparation requires an annotation label')
    for value in (annotation_flank, region_flank):
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError('SV source flanks must be nonnegative integers')
    ends = tuple(dict.fromkeys(end for pair in variant.junctions for end in pair))
    resolved = {end.contig: resolver.resolve(end.contig) for end in ends}
    limitations = [f'annotation_contig_unresolved:{name}: {resolution.diagnostic}'
                   for name, resolution in resolved.items() if resolution.resolved_name is None]
    transcripts: dict[str, tuple[Any, str]] = {}
    gene_names: set[str] = set()
    for end in ends:
        contig = resolved[end.contig].resolved_name
        if contig is None:
            continue
        for gene in genome.genes_at_locus(contig, max(1, end.position - annotation_flank),
                                         end.position + annotation_flank):
            for transcript in gene.transcripts:
                if transcript.id in transcripts and transcripts[transcript.id][1] != end.contig:
                    raise ValueError('One annotation transcript maps to conflicting source contig aliases')
                transcripts[transcript.id] = (transcript, end.contig)
        gene_names.update(gene.name for gene in genome.genes_at_locus(
            contig, max(1, end.position - 1), end.position + 1) if gene.name)
    references = []
    for identifier in sorted(transcripts):
        transcript, source_contig = transcripts[identifier]
        model, limitation = _reference_model(transcript, source_contig, genome.reference_name, annotation)
        if model is not None:
            references.append(model)
        if limitation is not None:
            limitations.append(limitation)
    if not references:
        limitations.append('No usable reference transcripts in breakpoint windows; exon-anchored clip seeding unavailable')
    hypotheses = (() if any(value.resolved_name is None for value in resolved.values()) else
                  tuple(hypotheses_from_varcode(record.for_annotation(genome, resolver).effects())))
    views = []
    for index, (first, second) in enumerate(variant.junctions, 1):
        for orientation, (left, right) in (('forward', (first, second)), ('reverse', (second, first))):
            regions = tuple((end.contig, max(0, end.position - 1 - region_flank),
                             end.position + region_flank) for end in (left, right))
            views.append(SVSourceView(f'{index}-{orientation}', _breakpoint(left, True),
                                      _breakpoint(right, False), regions))
    return PreparedSVSourceRecord(genome.reference_name, annotation, tuple(references), hypotheses,
                                  tuple(views), tuple(limitations), tuple(sorted(gene_names)),
                                  tuple(dict.fromkeys(f'{end.contig}:{end.position}-{end.position}' for end in ends)),
                                  annotation_flank)
