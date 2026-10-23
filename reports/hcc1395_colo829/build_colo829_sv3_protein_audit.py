"""Export COLO829 SV3 protein hypotheses and assess their proteoform evidence.

The input BAM is the source-derived selected-site fixture. Its genomic RNA
alignment includes the SV breakends, but not the neighboring CNIH3 exons.
"""

import argparse
from collections import defaultdict
import json
from pathlib import Path

from Bio.Seq import Seq
from orfdata import load_catalog
from pyensembl import EnsemblRelease
import pysam

from isovar.fusion import FusionBreakpoint
from isovar.read_collector import ReadCollector
from isovar.sv_rna import reconstruct_sv_rna
from isovar.sv_rna_orf_export import export_sv_rna_orfs, write_sv_rna_orfs

from .build_sv_junction_audit import COLO_BAM, EVENTS


HERE = Path(__file__).parent
REFERENCE_SLICE = HERE.parents[1] / "tests/data/colo829_sv3/reference_slice.json"
OUTPUT_PREFIX = HERE / "colo829_sv3_candidate_orfs"
EVENT = next(event for event in EVENTS if event[1] == "COLO829_SOMATIC_SV3_DEL")


def _reference_slice(path):
    reference = json.loads(Path(path).read_text())
    if (reference.get("schema") != "isovar.reference_slice.v1"
            or reference.get("contig") != "chr1"
            or reference.get("start") != 224594494
            or reference.get("end") != 224598926
            or len(reference.get("sequence", "")) != reference["end"] - reference["start"]):
        raise ValueError("Expected the pinned COLO829 SV3 GRCh38 reference slice")
    return reference


def _reconstruction(bam):
    _, call_id, _, donor, acceptor, _, _ = EVENT
    with pysam.AlignmentFile(str(bam), "rb") as alignment:
        return reconstruct_sv_rna(
            alignment, event_id=call_id, reference_name="GRCh38",
            sample_id="COLO829", source="SRA:SRR8615617",
            donor=FusionBreakpoint(*donor), acceptor=FusionBreakpoint(*acceptor),
            regions=[], references=[], event_provenance={
                "dna_call_id": call_id, "call_set": "PaulinPaz_2024_v2"},
            read_collector=ReadCollector(min_mapping_quality=0, use_soft_clipped_bases=True),
            max_records=5000, max_paths=10, max_queries=100,
            assemble=True, include_regional_candidates=False, dense_support=True)


def _placed_sequence_matches(path, reference):
    for block in path["blocks"]:
        a, b = block["reference_start"] - reference["start"], block["reference_end"] - reference["start"]
        if (block["contig"] != reference["contig"] or block["strand"] != "+"
                or a < 0 or b > len(reference["sequence"])
                or path["sequence"][block["query_start"]:block["query_end"]] != reference["sequence"][a:b]):
            raise ValueError("Assembled placed bases disagree with the pinned GRCh38 reference")


def _junction_ninemers(protein, start, junction):
    """9-mers containing both placed flanks and all unplaced junction bases."""
    left, right = junction["query_interval"]
    return [protein[i:i + 9] for i in range(len(protein) - 8)
            if start + 3 * i < left and start + 3 * (i + 9) > right]


def _protein_matches(genome, peptides):
    matches = defaultdict(list)
    for protein_id in genome.protein_ids():
        sequence = genome.protein_sequence(protein_id)
        for peptide in peptides:
            if peptide in sequence:
                matches[peptide].append(protein_id)
    return matches


def _junction_read_qualities(bam, names, unplaced):
    """Minimum sequenced base quality of the same unplaced junction in each read."""
    reverse = str(Seq(unplaced).reverse_complement())
    found = {}
    with pysam.AlignmentFile(str(bam), "rb") as alignment:
        for read in alignment.fetch(until_eof=True):
            if read.query_name not in names or read.is_supplementary or read.query_qualities is None:
                continue
            sequence = read.query_sequence or ""
            offset = sequence.find(unplaced)
            if offset < 0:
                offset = sequence.find(reverse)
            if offset >= 0:
                found[read.query_name] = min(read.query_qualities[offset:offset + len(unplaced)])
    return dict(sorted(found.items()))


def _locus_annotation(genome, catalog, reference):
    start, end = reference["start"] + 1, reference["end"]
    genes = genome.genes_at_locus("1", start, end)
    exons = [(transcript.id, exon.start, exon.end)
             for gene in genes for transcript in gene.transcripts for exon in transcript.exons
             if exon.start <= end and exon.end >= start]
    orfs = [orf.orf_id for orf in catalog.records if orf.contig in ("1", "chr1")
            and any(a <= end and b >= start for a, b in orf.coding_blocks)]
    return dict(genes=sorted({gene.name for gene in genes}),
                overlapping_ensembl87_exons=sorted(set(exons)),
                overlapping_bundled_orfdata_orfs=sorted(set(orfs)),
                bundled_orfdata_version=catalog.metadata["dataset_version"])


def assess(bam=COLO_BAM, reference_path=REFERENCE_SLICE, genome=None, catalog=None):
    """Assess candidate sequence without promoting RNA reads to protein evidence."""
    reference = _reference_slice(reference_path)
    result = _reconstruction(bam)
    export = export_sv_rna_orfs(result)
    genome = EnsemblRelease(87) if genome is None else genome
    catalog = load_catalog() if catalog is None else catalog
    locus = _locus_annotation(genome, catalog, reference)
    paths = {path["path_id"]: path for path in result["paths"]}
    rows, peptides = [], set()
    for candidate in export["candidates"]:
        for occurrence in candidate["occurrences"]:
            path = paths[occurrence["path_id"]]
            _placed_sequence_matches(path, reference)
            if len(occurrence["junctions"]) != 1:
                raise ValueError("SV3 ORF audit requires one crossed junction per candidate")
            junction = occurrence["junctions"][0]
            start = occurrence["query_interval"][0]
            first = path["blocks"][0]
            if not first["query_start"] <= start < first["query_end"]:
                raise ValueError("SV3 ORF start is outside the placed donor block")
            genomic_start = first["reference_start"] + start - first["query_start"]
            dna = reference["sequence"][genomic_start - reference["start"]:]
            reference_orf = str(Seq(dna[:len(dna) // 3 * 3]).translate(to_stop=True))
            protein = candidate["amino_acids"]
            prefix = next((i for i, pair in enumerate(zip(protein, reference_orf))
                           if pair[0] != pair[1]), min(len(protein), len(reference_orf)))
            ninemers = _junction_ninemers(protein, start, junction)
            peptides.update(ninemers)
            direct = next(j for j in path["junctions"] if j["query_interval"] == junction["query_interval"])
            rows.append(dict(
                candidate_id=candidate["candidate_id"], path_id=path["path_id"],
                amino_acids=protein, amino_acid_length=len(protein),
                orf_query_interval=occurrence["query_interval"],
                genomic_start_zero_based=genomic_start,
                junction_query_interval=junction["query_interval"],
                junction_ninemers=ninemers,
                shared_junction_fragments=direct["direct_support"]["fragments"],
                full_orf_junction_fragments=candidate["rna_support"]["fragments"],
                reference_genomic_orf=reference_orf,
                shared_reference_genomic_prefix_amino_acids=prefix,
                annotated_start=occurrence["annotated_start"],
                initiation_observed=candidate["initiation_observed"],
                translation_observed=candidate["translation_observed"],
                protein_status="exploratory_orf_mature_transcript_unresolved"))
    matches = _protein_matches(genome, peptides)
    for row in rows:
        row["ensembl87_protein_matches_by_junction_ninemer"] = {
            peptide: matches[peptide] for peptide in row["junction_ninemers"]}
    junctions = [junction for path in paths.values() for junction in path["junctions"]
                 if junction["relation"] == "breakpoint_junction"]
    names = sorted({result["observations"][key]["identity"][1]
                    for junction in junctions for key in junction["direct_observations"]})
    unplaced = sorted({junction["unplaced_bases"] for junction in junctions})
    if len(unplaced) != 1:
        raise ValueError("SV3 audit requires one observed unplaced junction sequence")
    qualities = _junction_read_qualities(bam, names, unplaced[0])
    if set(qualities) != set(names):
        raise ValueError("Not every original junction read contains the unplaced sequence with qualities")
    return dict(schema="isovar.colo829_sv3_protein_audit.v1", event_id=result["event_id"],
                original_rna_source=result["source"], reference_slice={
                    key: reference[key] for key in ("contig", "start", "end", "source", "source_fasta_gzip_sha256")},
                locus=locus, path_count=len(paths), assembly_limitations=result["limitations"],
                flanking_voting_fragments=max((path["sequence_evidence"]["voting_support"]["fragments"]
                                               for path in paths.values()), default=0),
                unplaced_junction_bases=unplaced[0], original_junction_read_names=names,
                junction_base_quality_min_by_read=qualities,
                candidates=sorted(rows, key=lambda row: (row["orf_query_interval"], row["candidate_id"])),
                proteoform_call="not_established", orf_export=export)


def _report(assessment, prefix):
    rows = assessment["candidates"]
    matched = sorted({peptide for row in rows for peptide, proteins in
                      row["ensembl87_protein_matches_by_junction_ninemer"].items() if proteins})
    lines = [
        "# COLO829 SV3: RNA ORFs and proteoform status", "",
        "**Proteoform call: not established.** The sequences below are exploratory",
        "ATG-to-stop translations of an assembled RNA path, not measured protein.", "",
        "The SMaHT/Paulin Paz DNA call is `COLO829_SOMATIC_SV3_DEL`. Three distinct",
        "original tumor-RNA fragments support the same breakpoint-compatible join",
        f"with `{assessment['unplaced_junction_bases']}` between its placed flanks.",
        "Their minimum sequenced base qualities over these 19 bases are " + ", ".join(
            f"{name}: {quality}" for name, quality in
            assessment["junction_base_quality_min_by_read"].items()) + ".", "",
        f"Isovar finds {len(rows)} alternative ATG reading frames on a 1,213-base",
        "assembled RNA context.",
        f"The {assessment['flanking_voting_fragments']} flanking voting fragments are regional",
        "coverage; they are not support for either complete isoform. The assembly",
        "reaches `extension_segment_limit`.", "",
        "| Candidate | Amino acids | Shared join fragments | Whole ORF witnesses | "
        "Hypothetical unrearranged genomic ORF | Shared amino-acid prefix | Join-crossing 9-mers |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for index, row in enumerate(rows, 1):
        lines.append(
            f"| {index} | {row['amino_acid_length']} | {row['shared_junction_fragments']} | "
            f"{row['full_orf_junction_fragments']} | {len(row['reference_genomic_orf'])} aa | "
            f"{row['shared_reference_genomic_prefix_amino_acids']} aa | "
            f"{', '.join('`' + peptide + '`' for peptide in row['junction_ninemers'])} |")
    lines.extend(["", "## Candidate amino-acid sequences", ""])
    for index, row in enumerate(rows, 1):
        lines.extend([f"Candidate {index} (`{row['candidate_id']}`):", "", "```text",
                      row["amino_acids"], "```", ""])
    lines += [
        "", "The unrearranged comparisons translate *intronic genomic sequence* from",
        "the same ATGs. They are controls for sequence divergence, not known",
        "reference proteins. " + ("These listed 9-mers also occur in Ensembl 87: " +
        ", ".join(matched) + "." if matched else
        "None of the listed 9-mers occurs in the pinned Ensembl 87 annotated protein sequences."),
        "This is not proteome-wide novelty for unannotated ORFs.", "",
        "Ensembl 87 places the interval within the CNIH3 gene span, but no",
        "annotated exon overlaps the assembled region. The bundled ORFdata",
        f"catalogue ({assessment['locus']['bundled_orfdata_version']}) has no overlapping",
        "coding block; optional ORF catalogues were not assessed. The selected-site",
        "BAM cannot establish a mature splice path from this join to a CNIH3 exon.",
        "Neither candidate has an annotated start or observed initiation. The",
        "lack of a whole-ORF witness limits phasing; no single read is required",
        "to span the assembled sequence. Ribosome and peptide data were not",
        "assessed. A predicted proteoform would require a supported mature RNA",
        "path and frame; protein evidence would validate its expression.", "",
        "The existing `isovar sv-rna --orf-output-prefix` export writes the",
        "[protein FASTA](" + prefix.name + ".protein.fasta),",
        "[nucleotide FASTA](" + prefix.name + ".nucleotide.fasta),",
        "[ORF JSON](" + prefix.name + ".json) and",
        "[ORF TSV](" + prefix.name + ".tsv). This assessment is",
        "[machine readable](" + prefix.name + ".assessment.json).",
        "To resolve one hypothesis, obtain original RNA reads through a",
        "path-discriminating exon join and compare normal/reference splicing;",
        "then seek frame-specific ribosome footprints or a junction-specific",
        "peptide. The RNA and DNA came from separate COLO829 cultures.", "",
    ]
    return "\n".join(lines)


def write_audit(prefix=OUTPUT_PREFIX, bam=COLO_BAM, reference_path=REFERENCE_SLICE):
    """Write the existing portable ORF export and a locus-specific assessment."""
    prefix = Path(prefix)
    assessment = assess(bam=bam, reference_path=reference_path)
    write_sv_rna_orfs(assessment.pop("orf_export"), prefix)
    path = Path(str(prefix) + ".assessment.json")
    path.write_text(json.dumps(assessment, sort_keys=True, indent=2) + "\n")
    Path(str(prefix) + ".assessment.md").write_text(_report(assessment, prefix))
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, default=COLO_BAM)
    parser.add_argument("--reference-slice", type=Path, default=REFERENCE_SLICE)
    parser.add_argument("--output-prefix", type=Path, default=OUTPUT_PREFIX)
    args = parser.parse_args()
    print(write_audit(args.output_prefix, args.bam, args.reference_slice))


if __name__ == "__main__":
    main()
