"""Versioned reference evidence for human upstream open reading frames.

Genomic coordinates are 1-based inclusive on the explicitly named assembly.
These are reference observations, not evidence for a patient's mutant protein.
"""

from copy import deepcopy
from dataclasses import asdict, dataclass
import gzip
import hashlib
from importlib.resources import files
import json
from pathlib import Path
from typing import Optional


EVIDENCE_METHODS = ("ribosome", "shotgun", "hla", "ms_unclassified", "ms")


def _positive_integer(value, name):
    if type(value) is not int or value < 1:
        raise ValueError("%s must be a positive integer" % name)


def _contig(value):
    value = str(value)
    if value.startswith("chr"):
        value = value[3:]
    return "MT" if value in ("M", "MT") else value


@dataclass(frozen=True)
class UORFPeptide:
    """A reported MS peptide, with every compatible protein interval.

    ``protein_intervals`` are 1-based inclusive. ``mapping_count`` and
    spectrum counts have the scope of the named source, not all human proteins
    or independent biological specimens. Missing counts remain unknown.
    """

    method: str
    sequence: str
    protein_intervals: tuple
    source: str
    accession: str
    mapping_count: Optional[int] = None
    spectrum_count: Optional[int] = None
    datasets: tuple = ()
    sequence_match: str = "exact"
    attribution: str = "source_orf"
    quality: str = "reported"
    review_ratings: tuple = ()
    spectrum_ids: tuple = ()

    @property
    def unique_in_source(self):
        """Whether the source reports one protein mapping (not global novelty)."""
        return self.mapping_count == 1

    @property
    def is_supporting(self):
        """Reported or convincingly reviewed, excluding low/mixed reviews."""
        return self.quality in ("reported", "reviewed_support")


@dataclass(frozen=True)
class UORFRecord:
    """One reported, spliced reference ORF, not a gene-wide evidence label."""

    orf_id: str
    gene_id: str
    gene_name: str
    transcript_id: str
    contig: str
    strand: str
    coding_blocks: tuple
    protein_sequence: str
    biotype: str
    ribosome_studies: tuple
    peptides: tuple
    genome_build: str = "GRCh38"

    def to_dict(self):
        """Return a serializable copy of this record."""
        return asdict(self)

    def has_evidence(self, method, *, require_unique=False, require_reviewed=False):
        """Select reference observations by assay, without clinical promotion."""
        if method not in EVIDENCE_METHODS:
            raise ValueError("Unknown evidence method: %s" % method)
        if method == "ribosome":
            if require_unique or require_reviewed:
                raise ValueError("MS filters do not apply to ribosome evidence")
            return bool(self.ribosome_studies)
        return any(
            (method == "ms" or p.method == method)
            and p.is_supporting
            and (not require_unique or p.unique_in_source)
            and (not require_reviewed or p.quality == "reviewed_support")
            for p in self.peptides)

    def map_genomic_position(self, position, *, genome_build):
        """Map one genomic base into the spliced ORF and reference MS coverage.

        Parameters
        ----------
        position : int
            A 1-based genomic position on this record's contig.
        genome_build : str
            Must equal the record's assembly. No implicit liftover is done.

        Returns
        -------
        dict or None
            Coding nucleotide offset (0-based), residue position (1-based),
            reference residue, start/stop/body role and covering reference
            peptides. An intronic or outside position returns None. This does
            not compute a mutant sequence or demonstrate mutant translation.
        """
        _positive_integer(position, "position")
        if genome_build != self.genome_build:
            raise ValueError("Expected genome build %s, got %s" % (self.genome_build, genome_build))
        offset = 0
        blocks = self.coding_blocks if self.strand == "+" else self.coding_blocks[::-1]
        for start, end in blocks:
            if start <= position <= end:
                offset += position - start if self.strand == "+" else end - position
                aa = offset // 3 + 1
                stop = aa > len(self.protein_sequence)
                covered = [p for p in self.peptides
                           if any(a <= aa <= b for a, b in p.protein_intervals)]
                return dict(
                    orf_id=self.orf_id, genome_build=self.genome_build,
                    coding_nucleotide_offset=offset,
                    amino_acid_position=aa,
                    reference_amino_acid="*" if stop else self.protein_sequence[aa - 1],
                    role="stop_codon" if stop else "start_codon" if offset < 3 else "body",
                    reference_peptides=[asdict(p) for p in covered],
                    reference_ms_support=any(p.is_supporting for p in covered),
                    mutant_translation_evidence="not_assessed")
            offset += end - start + 1
        return None


class UORFCatalog:
    """Offline, assembly-specific uORF records and source provenance."""

    def __init__(self, records, metadata):
        self.records = tuple(records)
        self._metadata = deepcopy(metadata)
        self.genome_build = metadata["genome_build"]
        self._by_id = {r.orf_id: r for r in self.records}
        if len(self._by_id) != len(self.records):
            raise ValueError("Duplicate ORF identifiers")
        for r in self.records:
            if r.genome_build != self.genome_build or r.strand not in ("+", "-"):
                raise ValueError("Invalid assembly or strand for %s" % r.orf_id)
            if r.biotype not in ("uORF", "uoORF"):
                raise ValueError("Invalid upstream ORF/evidence for %s" % r.orf_id)
            if (not r.coding_blocks or not r.protein_sequence
                    or any(a < 1 or b < a for a, b in r.coding_blocks)
                    or any(b >= c for (_, b), (c, _) in zip(r.coding_blocks, r.coding_blocks[1:]))
                    or sum(b - a + 1 for a, b in r.coding_blocks) != 3 * (len(r.protein_sequence) + 1)):
                raise ValueError("Invalid spliced coding blocks for %s" % r.orf_id)
            for p in r.peptides:
                if p.method not in ("shotgun", "hla", "ms_unclassified") or not p.protein_intervals:
                    raise ValueError("Invalid peptide evidence for %s" % r.orf_id)
                if (p.quality not in ("reported", "reviewed_support", "reviewed_low_quality", "reviewed_mixed")
                        or p.attribution not in ("source_orf", "sequence_compatible")
                        or p.sequence_match not in ("exact", "I_L_equivalent")
                        or not p.source or not p.accession):
                    raise ValueError("Invalid peptide provenance for %s" % r.orf_id)
                for count in (p.mapping_count, p.spectrum_count):
                    if count is not None:
                        _positive_integer(count, "evidence count")
                for curator, rating in p.review_ratings:
                    if not curator or type(rating) is not int or not 1 <= rating <= 5:
                        raise ValueError("Invalid spectrum review rating")
                if p.quality != "reported" and not p.review_ratings:
                    raise ValueError("Reviewed evidence requires review ratings")
                for a, b in p.protein_intervals:
                    if not 1 <= a <= b <= len(r.protein_sequence):
                        raise ValueError("Peptide interval outside %s" % r.orf_id)
                    if r.protein_sequence[a - 1:b].replace("I", "L") != p.sequence.replace("I", "L"):
                        raise ValueError("Peptide sequence does not match %s" % r.orf_id)

    @property
    def metadata(self):
        """Return a copy of dataset version, sources, scope and provenance."""
        return deepcopy(self._metadata)

    def get(self, orf_id):
        """Return an exact ORF ID; unknown IDs raise KeyError."""
        return self._by_id[orf_id]

    def query(self, *, gene=None, transcript_id=None, contig=None, start=None,
              end=None, genome_build=None, evidence=None, require_unique=False,
              require_reviewed=False, biotype=None):
        """Filter records; genomic intervals are 1-based inclusive and spliced.

        A locus query requires an explicit assembly. Gene and transcript lookup
        is a discovery aid, not evidence transfer between transcript isoforms.
        Transcript identifiers match exactly as supplied in the source; these
        are unversioned identifiers. Empty results do not prove no translation.
        ``require_unique`` means unique within the evidence source's mapping
        scope, not unique in every normal human proteome.
        """
        locus = any(x is not None for x in (contig, start, end))
        if locus:
            if contig is None or start is None or end is None or genome_build is None:
                raise ValueError("Locus queries require contig, start, end and genome_build")
            _positive_integer(start, "start")
            _positive_integer(end, "end")
            if end < start:
                raise ValueError("end must be >= start")
            contig = _contig(contig)
        if genome_build is not None and genome_build != self.genome_build:
            raise ValueError("Expected genome build %s, got %s" % (self.genome_build, genome_build))
        if evidence is not None and evidence not in EVIDENCE_METHODS:
            raise ValueError("Unknown evidence method: %s" % evidence)
        if require_unique and evidence not in EVIDENCE_METHODS[1:]:
            raise ValueError("require_unique applies to MS evidence only")
        if require_reviewed and evidence not in EVIDENCE_METHODS[1:]:
            raise ValueError("require_reviewed applies to MS evidence only")
        if biotype is not None and biotype not in ("uORF", "uoORF"):
            raise ValueError("Unknown upstream ORF biotype: %s" % biotype)
        return tuple(r for r in self.records
                     if (gene is None or gene in (r.gene_name, r.gene_id))
                     and (transcript_id is None or r.transcript_id == transcript_id)
                     and (not locus or _contig(r.contig) == contig
                          and any(a <= end and b >= start for a, b in r.coding_blocks))
                     and (biotype is None or r.biotype == biotype)
                     and (evidence is None or r.has_evidence(
                         evidence, require_unique=require_unique, require_reviewed=require_reviewed)))


def load_uorf_catalog(path=None, *, expected_sha256=None):
    """Load the packaged reference offline, or a compatible JSON/JSON.gz file.

    No genome downloads, spreadsheet dependencies or LLM are required. The
    catalogue stores protein sequences, not genomic REF alleles; callers must
    validate alleles against their own assembly before reconstructing mutants.
    """
    raw = (files("isovar").joinpath("data/uorf-evidence/catalog.json.gz").read_bytes()
           if path is None else Path(path).read_bytes())
    if path is None:
        manifest = json.loads(files("isovar").joinpath("data/uorf-evidence/manifest.json").read_text())
        if hashlib.sha256(raw).hexdigest() != manifest["catalog_sha256"]:
            raise ValueError("Packaged uORF catalogue checksum mismatch")
    if expected_sha256 is not None and hashlib.sha256(raw).hexdigest() != expected_sha256:
        raise ValueError("uORF catalogue checksum mismatch")
    if raw.startswith(b"\x1f\x8b"):
        raw = gzip.decompress(raw)
    data = json.loads(raw)
    if data.get("schema") != "isovar.uorf_evidence.v1":
        raise ValueError("Unsupported uORF evidence schema")
    if data["metadata"].get("coordinates") != "1-based-inclusive":
        raise ValueError("Unsupported uORF coordinate convention")
    records = []
    for row in data["records"]:
        row = dict(row)
        row["coding_blocks"] = tuple(tuple(x) for x in row["coding_blocks"])
        row["ribosome_studies"] = tuple(row["ribosome_studies"])
        peptides = []
        for p in row["peptides"]:
            p = dict(p)
            p["protein_intervals"] = tuple(tuple(x) for x in p["protein_intervals"])
            p["datasets"] = tuple(p.get("datasets", ()))
            p["review_ratings"] = tuple(tuple(x) for x in p.get("review_ratings", ()))
            p["spectrum_ids"] = tuple(p.get("spectrum_ids", ()))
            peptides.append(UORFPeptide(**p))
        row["peptides"] = tuple(peptides)
        records.append(UORFRecord(**row))
    return UORFCatalog(records, data["metadata"])
