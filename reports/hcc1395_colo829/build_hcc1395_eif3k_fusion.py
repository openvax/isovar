"""Replay the public HCC1395 EIF3K--CYP39A1 RNA fusion from original reads.

The input BAM contains the 185 unmapped Illumina records in INTEGRATE's public
example. Reference models are pinned separately; no annotation service is used
at report or test time. Only exact original reads spanning the exon join enter
the local assembly, and every reported base must be covered by those reads.
"""

import argparse
from dataclasses import asdict
import gzip
import json
from pathlib import Path
import shutil
import tempfile

import pysam

from isovar.fusion import (
    FusionBlock, FusionBreakpoint, FusionRead, FusionReference, FusionTranscript,
    reconstruct_fusion,
)


ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests/data/hcc1395_eif3k_cyp39a1"
DONOR_ID = "ENST00000248342"
ACCEPTOR_ID = "ENST00000275016"
DONOR_EDGE = 39123318  # GRCh37 interbase, end of EIF3K exon 6
ACCEPTOR_EDGE = 46607405  # GRCh37 interbase, upper edge of CYP39A1 exon 3


def load_references(path):
    """Load sequence-verified, versioned exon and CDS models."""
    with gzip.open(path, "rt") as handle:
        data = json.load(handle)
    return tuple(FusionReference(**model) for model in data["models"])


def _fusion_blocks(donor, acceptor):
    """Lay out the two observed spliced partners in RNA orientation."""
    left = tuple(e for e in donor.exons if e[1] <= DONOR_EDGE)
    right = tuple(e for e in acceptor.exons[::-1] if e[1] <= ACCEPTOR_EDGE)
    if not left or not right or left[-1][1] != DONOR_EDGE or right[0][1] != ACCEPTOR_EDGE:
        raise ValueError("The pinned transcripts do not contain the reported exon join")
    blocks = []
    cursor = 0
    for reference, exons in ((donor, left), (acceptor, right)):
        for start, end in exons:
            blocks.append(FusionBlock(cursor, cursor + end - start, reference.contig,
                                      start, end, reference.strand))
            cursor += end - start
    return tuple(blocks), sum(end - start for start, end in left)


def _subblocks(blocks, start, end):
    """Clip mapped exons to an observed cDNA interval."""
    result = []
    for block in blocks:
        left, right = max(start, block.query_start), min(end, block.query_end)
        if left >= right:
            continue
        if block.strand == "+":
            genomic_start = block.reference_start + left - block.query_start
            genomic_end = genomic_start + right - left
        else:
            genomic_end = block.reference_end - (left - block.query_start)
            genomic_start = genomic_end - (right - left)
        result.append(FusionBlock(left - start, right - start, block.contig,
                                  genomic_start, genomic_end, block.strand))
    return tuple(result)


def path_coordinates(blocks):
    """Map each transcript-oriented query base to one reference position."""
    positions = {}
    for block in blocks:
        for query in range(block.query_start, block.query_end):
            genomic = (block.reference_start + query - block.query_start if block.strand == "+"
                       else block.reference_end - 1 - (query - block.query_start))
            positions[query] = block.contig, genomic
    return positions


def mapped_coordinates_match(record, start, oriented, positions):
    """Check every original aligned base against the nominated exon path."""
    reverse = oriented != record.query_sequence
    pairs = record.get_aligned_pairs(matches_only=True)
    if len(pairs) != len(oriented):
        return False
    return all(positions.get(start + (len(oriented) - 1 - query if reverse else query)) ==
               (record.reference_name, genomic)
               for query, genomic in pairs)


def assemble(original_bam, references, mapped_tiling_bam=None, min_flank=10, min_overlap=20):
    """Return the observed local or complete coding interval of the fusion.

    The reference isoforms nominate a junction and orient each original read;
    they do not contribute any unwitnessed base to the returned sequence.
    Repeated or imperfect placements are excluded from this strict example.
    """
    donor = next(ref for ref in references if ref.transcript_id == DONOR_ID)
    acceptor = next(ref for ref in references if ref.transcript_id == ACCEPTOR_ID)
    blocks, junction = _fusion_blocks(donor, acceptor)
    positions = path_coordinates(blocks)
    acceptor_offset = acceptor.spliced_offset(ACCEPTOR_EDGE - 1)
    full = donor.sequence[:junction] + acceptor.sequence[acceptor_offset:]
    if len(full) != blocks[-1].query_end:
        raise ValueError("Fusion sequence and exon blocks disagree")
    matches = []
    with pysam.AlignmentFile(str(original_bam), check_sq=False) as bam:
        source_count = 0
        for record in bam.fetch(until_eof=True):
            source_count += 1
            seq = record.query_sequence
            candidates = []
            for oriented in (seq, pysam.reverse_complement(seq)):
                start = full.find(oriented)
                while start >= 0:
                    if min(junction - start, start + len(oriented) - junction) >= min_flank:
                        candidates.append((start, oriented))
                    start = full.find(oriented, start + 1)
            if len(candidates) == 1:
                start, oriented = candidates[0]
                matches.append((record.query_name, start, oriented,
                                "INTEGRATE:test-data/unmapped.bam",
                                1 if record.is_read1 else 2 if record.is_read2 else 0))
    if not matches:
        raise ValueError("No uniquely placed original reads span the nominated RNA join")
    mapped_count = 0
    if mapped_tiling_bam is None:
        lo = min(start for _, start, _, _, _ in matches)
        hi = max(start + len(seq) for _, start, seq, _, _ in matches)
    else:
        if donor.cds_start is None or acceptor.cds_end is None:
            raise ValueError("Complete coding reconstruction needs verified CDS models")
        lo = donor.cds_start
        hi = junction + acceptor.cds_end - acceptor_offset
        with pysam.AlignmentFile(str(mapped_tiling_bam)) as bam:
            for record in bam.fetch(until_eof=True):
                mapped_count += 1
                if record.is_unmapped or record.mapping_quality < 20:
                    raise ValueError("The mapped tiling contains an unqualified alignment")
                candidates = []
                for oriented in (record.query_sequence,
                                 pysam.reverse_complement(record.query_sequence)):
                    start = full.find(oriented)
                    while start >= 0:
                        if (start < hi and start + len(oriented) > lo and
                                (start + len(oriented) <= junction and record.reference_name == "19" or
                                 start >= junction and record.reference_name == "6")):
                            candidates.append((start, oriented))
                        start = full.find(oriented, start + 1)
                if len(candidates) != 1:
                    raise ValueError("A mapped tiling read lacks one exact fusion placement")
                start, oriented = candidates[0]
                if not mapped_coordinates_match(record, start, oriented, positions):
                    raise ValueError("A mapped tiling read disagrees with the pinned exon path")
                matches.append((record.query_name, start, oriented,
                                "INTEGRATE:test-data/accepted_hits.bam",
                                1 if record.is_read1 else 2 if record.is_read2 else 0))
    # A coverage union alone can connect two reads that only abut. Require a
    # chain of observed sequence overlaps through the entire returned path.
    intervals = [(max(lo, start), min(hi, start + len(seq)))
                 for _, start, seq, _, _ in matches]
    reach = max((end for start, end in intervals if start <= lo), default=lo)
    while reach < hi:
        next_reach = max((end for start, end in intervals
                          if start <= reach - min_overlap and end > reach), default=reach)
        if next_reach == reach:
            raise ValueError("No read-overlap chain covers the fusion path")
        reach = next_reach
    coverage = [0] * (hi - lo)
    for _, start, seq, _, _ in matches:
        for position in range(max(lo, start) - lo, min(hi, start + len(seq)) - lo):
            coverage[position] += 1
    if not all(coverage):
        raise ValueError("The local fusion path has an unobserved base")
    provenance = dict(
        sample_id="HCC1395", method="reference-guided_exact_original_read_assembly",
        version="1", parameters=dict(minimum_exact_junction_flank_bases=min_flank,
                                     minimum_read_overlap_bases=min_overlap,
                                     original_unmapped_records_screened=source_count,
                                     mapped_tiling_records=mapped_count),
        source="INTEGRATE 2014 test-data.tar.gz / unmapped.bam and accepted_hits.bam",
        contig_id="EIF3K-ENST00000248342--CYP39A1-ENST00000275016",
    )
    fusion = FusionTranscript(
        event_id="EIF3K--CYP39A1", reference_name="GRCh37-lite",
        sequence=full[lo:hi], junction_start=junction - lo, junction_end=junction - lo,
        donor=FusionBreakpoint("19", DONOR_EDGE, "+"),
        acceptor=FusionBreakpoint("6", ACCEPTOR_EDGE, "-"),
        blocks=_subblocks(blocks, lo, hi), provenance=provenance,
    )
    reads = []
    for name, start, seq, source, mate in matches:
        fragment = name.rsplit("/", 1)[0]
        if not mate and name.endswith(("/1", "/2")):
            mate = int(name[-1])
        read_start, read_end = max(lo, start), min(hi, start + len(seq))
        reads.append(FusionRead(
            sample_id="HCC1395", library_id="INTEGRATE-2014-Illumina",
            fragment_id=fragment, read_id=name, source=source,
            source_query_start=read_start - start, cdna_start=read_start - lo,
            sequence=seq[read_start - start:read_end - start],
            blocks=_subblocks(blocks, read_start, read_end), mate_number=mate,
        ))
    return fusion, tuple(reads)


def exact_ont_junctions(fasta_path, fusion, flank=20):
    """Count independent ONT candidates with an exact junction k-mer."""
    if fasta_path is None:
        return None
    motif = fusion.sequence[fusion.junction_start - flank:fusion.junction_end + flank]
    if len(motif) != 2 * flank:
        raise ValueError("The observed fusion interval is too short for the ONT check")
    with gzip.open(fasta_path, "rt") as handle:
        sequences = [line.strip() for line in handle if line and not line.startswith(">")]
    reverse = pysam.reverse_complement(motif)
    return dict(screened_original_reads=len(sequences), exact_junction_reads=sum(
        motif in sequence or reverse in sequence for sequence in sequences),
        flank_bases=flank)


def report(original_bam, models_path, ont_fasta=None, mapped_tiling_bam=None):
    """Generate the Isovar result and a compact evidence summary."""
    references = load_references(models_path)
    fusion, reads = assemble(original_bam, references, mapped_tiling_bam)
    result = reconstruct_fusion(fusion, references, reads, min_fragments=2)
    proteins = [translation["amino_acids"] for translation in result["paths"][0]["translations"]]
    return dict(fusion=fusion, reads=reads, result=result, proteins=proteins,
                ont=exact_ont_junctions(ont_fasta, fusion))


def main(argv=None):
    """Write reviewable JSON, FASTA and Markdown from the pinned public reads."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original-bam", type=Path, default=FIXTURE / "original_unmapped.bam")
    parser.add_argument("--models", type=Path, default=FIXTURE / "models.json.gz")
    parser.add_argument("--ont-fasta", type=Path, default=FIXTURE / "ont_candidate_sequences.fasta.gz")
    parser.add_argument("--mapped-tiling-bam", type=Path, default=FIXTURE / "coding_tiling.bam")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--pdf", action="store_true", help="Also render the five-panel fusion PDF")
    args = parser.parse_args(argv)
    value = report(args.original_bam, args.models, args.ont_fasta, args.mapped_tiling_bam)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "fusion_input.json").write_text(json.dumps(dict(
        fusion=asdict(value["fusion"]), references=[asdict(r) for r in load_references(args.models)],
        reads=[asdict(read) for read in value["reads"]]), indent=2) + "\n")
    (args.output_dir / "fusion_result.json").write_text(json.dumps(value["result"], indent=2) + "\n")
    (args.output_dir / "candidate_proteins.fasta").write_text("".join(
        f">HCC1395_EIF3K_CYP39A1_{index} frame_status={value['result']['status']}\n{protein}\n"
        for index, protein in enumerate(value["proteins"], 1)))
    if args.pdf:
        from isovar.fusion_visualization import save_fusion_figures

        with tempfile.TemporaryDirectory() as scratch:
            directory = save_fusion_figures(value["result"], scratch,
                                            references=load_references(args.models), dpi=180)
            shutil.copyfile(directory / "all-figures.pdf", args.output_dir / "figures.pdf")
    support = value["result"]["paths"][0]["junctions"][0]["direct_support"]["fragments"]
    translation, = value["result"]["paths"][0]["translations"]
    junction_nucleotide = translation["junction_in_translated_cds"][0]
    junction_peptides = [item["sequence"] for item in translation["candidate_peptides"]
                         if len(item["sequence"]) == 9]
    ont = value["ont"]
    (args.output_dir / "report.md").write_text(
        "# HCC1395 EIF3K–CYP39A1 coding fusion candidate\n\n"
        f"- Original fusion-enriched Illumina records screened: "
        f"{value['fusion'].provenance['parameters']['original_unmapped_records_screened']}\n"
        f"- Exact junction-spanning fragments: {support}\n"
        f"- Original mapped Illumina reads in a coding-region overlap chain: "
        f"{value['fusion'].provenance['parameters']['mapped_tiling_records']}\n"
        f"- Fully observed coding RNA path: {len(value['fusion'].sequence)} bases\n"
        f"- Original ONT candidate reads with exact 20+20-base junction sequence: "
        f"{ont['exact_junction_reads']}/{ont['screened_original_reads']}\n"
        f"- Isovar frame status: {value['result']['status']}\n"
        f"- Candidate protein lengths: {', '.join(str(len(p)) for p in value['proteins'])} amino acids\n"
        f"- Fusion junction: coding nucleotide {junction_nucleotide}, within amino-acid codon "
        f"{junction_nucleotide // 3 + 1}\n"
        f"- Example junction-spanning 9-mers: {', '.join(junction_peptides[:5])}\n\n"
        "Every base of the coding path is witnessed by an exact original read, and "
        "the selected reads form a chain with at least 20 bases of overlap. "
        "The exon path is reference guided; six independent ONT candidate reads "
        "also contain the exact junction sequence. For this supplied coding path, "
        "Isovar finds one annotated frame "
        "with an observed start and stop. Alternative RNA isoforms remain possible. "
        "RNA evidence does not by itself establish DNA causation or protein "
        "expression. See candidate_proteins.fasta for the sequence.\n\n"
        + ("[Five-panel fusion figure](figures.pdf) shows the junction, annotated "
           "exon context, protein and example peptides.\n\n" if args.pdf else "") +
        "The Illumina input is a fusion-enriched subset, so these counts are not "
        "whole-sample prevalence estimates. The ONT candidates came from a local "
        "two-chromosome alignment and are checked here by exact original query "
        "sequence at the join. The [INTEGRATE public example]"
        "(https://sourceforge.net/p/integrate-fusion/wiki/example/) supplies the "
        "Illumina RNA and tumor/normal DNA sources; its [primary publication]"
        "(https://pmc.ncbi.nlm.nih.gov/articles/PMC4691743/) describes the "
        "multiple EIF3K–CYP39A1 isoforms. The independent ONT run is "
        "[SRR23456709](https://www.ncbi.nlm.nih.gov/sra/SRR23456709).\n")


if __name__ == "__main__":
    main()
