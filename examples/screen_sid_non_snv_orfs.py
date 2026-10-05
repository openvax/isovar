"""Screen acquired Sid indels for RNA-derived shifted coding sequence.

No new source downloads. Protein windows require observed RNA and an annotated
frame; full-window witnesses are counted separately from contributing reads.
"""

import argparse
from hashlib import sha256
import logging
from pathlib import Path

from Bio.Seq import Seq
import pysam
from pyensembl import EnsemblRelease
from varcode import Variant

from isovar import ProteinSequenceCreator, ReadCollector, run_isovar
from isovar.read_identity import fragment_ids, count_reads
from examples.sid_sv_audit.inventory import digest, identity, read_json, write_json
from examples.sid_sv_audit.reconstruct import reconstruction_deadline, implementation_identity


def full_window_support(protein):
    """Check the complete translated interval at the focal allele in each read."""
    witnesses, quality_witnesses, unknown_quality, checked = {}, {}, set(), 0
    for translation in protein.translations:
        orf = translation.variant_orf
        offset = orf.offset_to_first_complete_codon
        nt = orf.cdna_sequence[offset:offset + 3 * (len(protein.amino_acids) + int(protein.ends_with_stop_codon))]
        aa = str(Seq(nt).translate(table=2 if translation.reference_context.mitochondrial else 1))
        expected = protein.amino_acids + ("*" if protein.ends_with_stop_codon else "")
        if aa != expected:
            raise ValueError("Independent translation disagrees")
        checked += 1
        for read in translation.reads:
            reverse = translation.reference_context.strand == "-"
            sequence = str(Seq(read.sequence).reverse_complement()) if reverse else read.sequence
            start = (len(read.suffix) if reverse else len(read.prefix)) - (orf.variant_cdna_interval_start - offset)
            if start < 0 or sequence[start:start + len(nt)] != nt:
                continue
            for key in fragment_ids([read]):
                fid = sha256(repr(key).encode()).hexdigest()
                witnesses[fid] = dict(fragment_id=fid, nucleotide_sha256=sha256(nt.encode()).hexdigest())
                scores = read.quality_scores
                if scores is None:
                    unknown_quality.add(fid)
                else:
                    scores = tuple(reversed(scores)) if reverse else scores
                    interval = scores[start:start + len(nt)]
                    if any(q is None for q in interval):
                        unknown_quality.add(fid)
                    elif min(interval) >= 20:
                        quality_witnesses[fid] = witnesses[fid]
    return dict(full_window_fragments=len(witnesses), q20_full_window_fragments=len(quality_witnesses),
                fragments_with_unavailable_qualities=len(unknown_quality),
                independent_translation_checks=checked, witnesses=[witnesses[k] for k in sorted(witnesses)])


def screen(directory, output):
    selection = read_json(directory / "priority-selection.json")
    genome = EnsemblRelease(115)
    installed = genome.inspect_data()
    if not installed["installed"]:
        raise ValueError("Ensembl 115 annotation must already be installed and indexed")
    # Pin both source files and the indexes actually consumed by PyEnsembl.
    reference_pins = [dict(path=info.path, sha256=digest(info.path))
                      for _, info in sorted(installed["files"].items())]
    engine = implementation_identity()
    creator = ProteinSequenceCreator(protein_sequence_length=120, variant_sequence_assembly=False,
                                     max_protein_sequences_per_variant=0)
    collector = ReadCollector(min_mapping_quality=20, merge_overlapping_fragments=False,
                              use_secondary_alignments=False)
    for sid in selection["source_ids"]:
        for receipt_path in sorted((directory / "small-variants" / sid).glob("*.json")):
            receipt = read_json(receipt_path)
            if "bam" not in receipt:
                continue
            bam = receipt["bam"]
            if digest(bam["path"]) != bam["sha256"]:
                raise ValueError("Small-variant input changed")
            for name in receipt["request"]["targets"]:
                row = selection["small_variants"][name]
                variant = row["variant"]
                if not (any(len(ref) != len(alt) for _, _, ref, alt in variant["alleles"]) or
                        "splice_" in variant["annotations"].get("consequence", "")):
                    continue
                for number, allele in enumerate(variant["alleles"]):
                    request = dict(allele=allele, target_id=name, source_id=sid, bam_sha256=bam["sha256"],
                                   selection_sha256=digest(directory / "priority-selection.json"),
                                   engine_id=engine, parameters=creator.settings(), min_mapping_quality=20,
                                   merge_overlapping_fragments=False, use_secondary_alignments=False,
                                   convert_ucsc_contig_names=True, code_sha256=digest(__file__), max_seconds=90)
                    request["reference_source_pins"] = reference_pins
                    path = output / sid / (identity(request) + ".json.gz")
                    if path.exists():
                        if read_json(path)["request"] != request:
                            raise ValueError("Non-SNV screen request changed")
                        continue
                    record = dict(request=request, acquisition_status=receipt["status"],
                                  acquisition_receipt_sha256=digest(receipt_path),
                                  input_limitations=receipt.get("receipt", {}).get("limits", []))
                    try:
                        with reconstruction_deadline(90), pysam.AlignmentFile(bam["path"]) as handle:
                            annotation_variant = Variant(*allele, genome=genome, convert_ucsc_contig_names=True)
                            # Validate the annotation contig before absence of a frame is interpreted.
                            annotation_variant.transcripts
                            result, = run_isovar([annotation_variant], handle,
                                                read_collector=collector, protein_sequence_creator=creator)
                            proteins = []
                            for p in result.sorted_protein_sequences:
                                proteins.append(dict(amino_acids=p.amino_acids, frameshift=p.frameshift,
                                    mutation_interval=[p.mutation_start_idx, p.mutation_end_idx],
                                    mutant_amino_acids=p.mutant_amino_acids,
                                    ends_with_stop_codon=p.ends_with_stop_codon,
                                    contributing_fragments=p.num_supporting_fragments,
                                    transcript_ids=sorted(p.transcript_ids),
                                    **full_window_support(p)))
                            record.update(status="screened", allele_support={a: dict(
                                fragments=len(fragment_ids(getattr(result, a + "_reads"))),
                                reads=count_reads(getattr(result, a + "_reads"))) for a in ("ref", "alt", "other")},
                                proteins=proteins)
                    except (ValueError, TimeoutError) as error:
                        record.update(status="screen_error", error=dict(type=type(error).__name__, message=str(error)))
                    write_json(path, record)
                    print(sid[:12], name, number, record["status"], len(record.get("proteins", [])), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    logging.disable(logging.WARNING)
    screen(args.directory, args.output)


if __name__ == "__main__":
    main()
