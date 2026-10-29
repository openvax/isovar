"""Explicit alignment references and PyEnsembl annotations for antigen RNA.

GDC's combined FASTA is a recommended alignment profile, not an Ensembl
release. No function here downloads data or changes existing genome defaults.
"""

from dataclasses import dataclass
from hashlib import md5, sha256
from importlib.resources import files
import json
from pathlib import Path

from pyensembl import Genome
from pyensembl.download_cache import cache_root
from varcode import ContigResolver

from .contig_resolution import UnresolvedContigError, require_resolved_contig

from .audit_bundle import digest
from .dna import reverse_complement_dna
from .fusion import FusionReference

GDC_REFERENCE = "GRCh38.d1.vd1"


def gdc_reference_dictionary():
    """Return GDC's published contig lengths/MD5s and acquisition provenance."""
    return json.loads(files("isovar").joinpath("data/references/gdc-grch38.d1.vd1.json").read_text())


def gdc_viral_panel():
    """Return GDC's explicitly listed viral contigs and historical source IDs."""
    return json.loads(files("isovar").joinpath("data/references/gdc-grch38.d1.vd1-viruses.json").read_text())


def _annotation_cache(gtf, transcript_fasta=None, protein_fasta=None, *, directory=None):
    """Reuse identical local annotations without same-basename index collisions."""
    identities = [digest(Path(path).expanduser()) if path is not None else None
                  for path in (gtf, transcript_fasta, protein_fasta)]
    key = sha256(json.dumps(identities).encode()).hexdigest()
    base = Path(directory).expanduser() if directory is not None else Path(cache_root()) / "isovar"
    return str(base / key)


def gdc_genome(*, genome_fasta, gtf, annotation_version, cache_directory=None,
               annotation_name="GENCODE", transcript_fasta=None, protein_fasta=None):
    """Create a custom PyEnsembl genome using local GDC DNA and human GTF.

    ``annotation_version`` is the supplied GENCODE release, independently of
    the GDC alignment profile. Call ``index()`` before querying transcripts.
    The annotation need not describe every viral/decoy contig in the FASTA.
    ``cache_directory`` is a cache root; annotation content selects its child
    directory, so samples can reuse indexes without reusing stale source data.
    """
    return Genome(
        reference_name="GRCh38", annotation_name=annotation_name, annotation_version=annotation_version,
        gtf_path_or_url=str(Path(gtf).expanduser().resolve()),
        genome_fasta_path_or_url=str(Path(genome_fasta).expanduser().resolve()),
        transcript_fasta_paths_or_urls=[] if transcript_fasta is None else [str(Path(transcript_fasta).expanduser().resolve())],
        protein_fasta_paths_or_urls=[] if protein_fasta is None else [str(Path(protein_fasta).expanduser().resolve())],
        cache_directory_path=_annotation_cache(gtf, transcript_fasta, protein_fasta, directory=cache_directory))


def contig_in_reference(contig, names):
    """Resolve an exact name or unambiguous alias through the shared vocabulary.

    ``names`` accepts a reused ContigResolver or an iterable of reference names.
    Missing and ambiguous contigs raise the common UnresolvedContigError;
    resolving a name does not establish sequence identity or perform liftover.
    """
    resolver = names if isinstance(names, ContigResolver) else ContigResolver(names)
    return require_resolved_contig(resolver.resolve(contig))


def validate_alignment_reference(bam, genome, *, profile=None, contigs=()):
    """Check BAM/FASTA lengths and MD5s without guessing the alignment assembly.

    A named GDC profile requires its complete sequence dictionary in the BAM
    and FASTA. Target contig DNA is checked against GDC's published MD5, and
    BAM M5 values are checked wherever supplied. For custom references, missing
    BAM M5 remains explicitly unverified; equal lengths do not prove identity.
    """
    fasta = genome.fasta
    if fasta is None:
        raise ValueError("A local genomic FASTA is required for isoform reconstruction")
    header = bam.header.to_dict()
    rows = header.get("SQ", [])
    lengths = {row["SN"]: row["LN"] for row in rows}
    if len(lengths) != len(rows):
        raise ValueError("Duplicate BAM reference contig names")
    expected = None
    if profile is not None:
        if profile != GDC_REFERENCE:
            raise ValueError("Unknown alignment profile %r; omit profile for a custom reference" % profile)
        if genome.reference_name != "GRCh38":
            raise ValueError("The GDC profile requires GRCh38 annotation")
        expected = {row["name"]: row for row in gdc_reference_dictionary()["contigs"]}
        if lengths != {name: row["length"] for name, row in expected.items()}:
            raise ValueError("BAM sequence dictionary does not match GDC GRCh38.d1.vd1")
        if set(fasta.keys()) != set(expected):
            raise ValueError("FASTA contigs do not match GDC GRCh38.d1.vd1")
    for row in rows:
        name = row["SN"]
        if name not in fasta or len(fasta[name]) != row["LN"]:
            raise ValueError("BAM/FASTA contig or length mismatch: %s" % name)
        if expected and row.get("M5") and row["M5"].lower() != expected[name]["md5"]:
            raise ValueError("BAM/GDC sequence checksum mismatch: %s" % name)
    checked, unverified = [], []
    resolver = ContigResolver(lengths, reference_name=genome.reference_name)
    by_name = {row["SN"]: row for row in rows}
    for contig in sorted(set(contigs)):
        name = contig_in_reference(contig, resolver)
        checksum = md5()
        for start in range(0, lengths[name], 1024 * 1024):
            checksum.update(str(fasta[name][start:min(lengths[name], start + 1024 * 1024)]).upper().encode("ascii"))
        actual = checksum.hexdigest()
        required = expected[name]["md5"] if expected else by_name[name].get("M5")
        if required and actual != required.lower():
            raise ValueError("Reference sequence checksum mismatch: %s" % name)
        (checked if required else unverified).append(name)
    return dict(profile=profile, assembly=genome.reference_name,
                annotation=dict(name=genome.annotation_name, version=genome.annotation_version),
                checked_sequence_contigs=checked, bam_sequence_identity_unverified=unverified,
                viral_annotation_inferred=False)


@dataclass(frozen=True)
class IsoformTarget:
    """One annotated competitor, with a gene identity and explicit coding status."""

    gene_id: str
    gene_name: str
    reference: FusionReference
    coding_status: str = "complete"

    def __post_init__(self):
        if not self.gene_id or not self.gene_name:
            raise ValueError("Isoform targets need a gene ID and name")
        if self.reference.cds_start is None and self.coding_status == "complete":
            object.__setattr__(self, "coding_status", "no_complete_annotated_CDS")


def isoform_targets_from_genome(genome, genes, *, contig_names):
    """Load every annotated transcript of explicit gene IDs or names.

    Reference cDNA is extracted from the selected genome's exons, avoiding a
    second potentially incompatible transcript FASTA. Missing sequence and
    unsupported coding models are returned as explicit ``unavailable`` rows.
    """
    selected, unavailable = {}, []
    resolver = ContigResolver(contig_names, reference_name=genome.reference_name)
    for name in genes:
        try:
            found = [genome.gene_by_id(name)]
        except ValueError:
            found = genome.genes_by_name(name)
        if not found:
            raise ValueError("Gene %r is absent from the supplied annotation" % name)
        for gene in found:
            selected[gene.id] = gene
    if not selected:
        raise ValueError("At least one explicit antigen target is required")
    targets = []
    for gene in sorted(selected.values(), key=lambda g: g.id):
        for transcript in sorted(gene.transcripts, key=lambda t: t.id):
            try:
                resolution = resolver.resolve(transcript.contig)
                contig = require_resolved_contig(resolution)
                exons = tuple(sorted((a - 1, b) for a, b in transcript.exon_intervals))
                sequence = "".join(genome.sequence(contig, a + 1, b) for a, b in exons)
                if transcript.strand == "-":
                    sequence = reverse_complement_dna(sequence)
                kwargs = dict(transcript_id=transcript.versioned_transcript_id, reference_name=genome.reference_name,
                              annotation="%s:%s" % (genome.annotation_name, genome.annotation_version),
                              contig=contig, strand=transcript.strand, exons=exons, sequence=sequence)
                reference = FusionReference(**kwargs)
                status = "no_complete_annotated_CDS"
                if resolution.normalized_name == 'MT':
                    status = "unsupported_mitochondrial_genetic_code"
                elif transcript.start_codon_complete and transcript.stop_codon_complete:
                    try:
                        reference = FusionReference(
                            **kwargs, cds_start=transcript.first_start_codon_spliced_offset,
                            cds_end=transcript.last_stop_codon_spliced_offset + 1)
                        status = "complete"
                    except ValueError as error:
                        status = "unsupported_CDS: " + str(error)
                targets.append(IsoformTarget(gene.id, gene.name or gene.id, reference, status))
            except UnresolvedContigError as error:
                unavailable.append(dict(gene_id=gene.id, transcript_id=transcript.id,
                                        reason=str(error), contig_resolution=error.to_dict()))
            except ValueError as error:
                unavailable.append(dict(gene_id=gene.id, transcript_id=transcript.id, reason=str(error)))
    return targets, unavailable


def reference_file_identity(genome_fasta, gtf):
    """Hash the supplied bytes; historical identities never follow live releases."""
    return {kind: dict(path=str(Path(path).resolve()), sha256=digest(path))
            for kind, path in (("genome_fasta", genome_fasta), ("gtf", gtf))}
