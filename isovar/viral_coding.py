"""Apply ORFdata coding references to observed RNA without donating sequence."""

from hashlib import sha256
from pathlib import Path

import orfdata
from orfdata import SourceRNA, TranscriptORF, viral_orfs, viral_references

from .dna import reverse_complement_dna


def viral_coding_references(genome, contigs):
    """Load exact-sequence ORFdata models and report unavailable coverage."""
    available = {ref.sequence_sha256: ref for ref in viral_references()}
    models, coverage = [], []
    for contig in contigs:
        sequence = str(genome.fasta[contig]).upper()
        fingerprint = sha256(sequence.encode("ascii")).hexdigest()
        source = available.get(fingerprint)
        row = dict(contig=contig, sequence_sha256=fingerprint,
                   status="available" if source else "unavailable",
                   reason=None if source else "no_sequence_qualified_ORFdata_reference")
        if source:
            loaded = viral_orfs(sequence, contig=contig)
            models.extend(loaded)
            row.update(reference=source.to_dict(), coding_models=len(loaded))
        else:
            row["coding_models"] = 0
        coverage.append(row)
    return models, dict(provider="orfdata", version=orfdata.__version__, contigs=coverage,
                        scope="source_CDS_subset_not_complete_transcript_annotation")


def _reference_path_status(observed, reference_positions, offsets):
    matches = [(q, offsets[p]) for q, p in enumerate(observed) if p in offsets]
    if not matches:
        return "unplaced"
    first, low = matches[0]
    last, high = matches[-1]
    if low > high or list(observed[first:last + 1]) != reference_positions[low:high + 1]:
        return "altered"
    if high < len(reference_positions) - 1 and last < len(observed) - 1:
        return "altered"
    return "compatible_at_observed_positions"


def viral_coding_predictions(sequence_id, sequence, positions, models):
    """Translate observed RNA from a source-placed frame, keeping partial ends.

    Source association anchors a frame; it does not establish a complete
    source transcript, the actual RNA strand, initiation or protein expression.
    Alternative splice/indel paths remain marked rather than assigned the
    source's original protein. Both orientations are considered explicitly.
    """
    readings = {}
    size = len(sequence)
    views = [("+", sequence, list(positions)),
             ("-", reverse_complement_dna(sequence),
              [(p[0], p[1], "-" if p[2] == "+" else "+") if p is not None else None
               for p in reversed(positions)])]
    for orientation, observed_sequence, observed_positions in views:
        for model in models:
            reference_positions = [(model.transcript.contig, p - 1, model.transcript.strand)
                                   for p in model.genomic_positions[model.start:model.end]]
            offsets = {p: i for i, p in enumerate(reference_positions)}
            for query in range(size - 2):
                offset = offsets.get(observed_positions[query])
                if (offset is None or offset % 3 or offset >= len(reference_positions) - 3
                        or observed_positions[query:query + 3] != reference_positions[offset:offset + 3]):
                    continue
                complete_start = offset == 0 and observed_sequence[query:query + 3] == model.start_codon
                rna = SourceRNA(id=sequence_id + ":" + orientation, sequence=observed_sequence,
                                reference=model.context, source="Isovar observed RNA path",
                                biotype="observed_RNA_fragment")
                predicted = TranscriptORF.at(
                    rna, query, orf_id=sequence_id + ":" + model.orf_id + ":" + orientation,
                    source="Observed RNA; frame anchored by " + model.source,
                    translation_table=model.genetic_code.table_id, start_complete=complete_start)
                end = predicted.end
                low, high = (query, end) if orientation == "+" else (size - end, size - query)
                protein = predicted.protein_sequence
                key = (low, high, orientation, protein, complete_start, predicted.stop_complete)
                identifier = "isovar_viral_orf_" + sha256(repr((sequence_id, key)).encode()).hexdigest()[:24]
                row = readings.setdefault(key, dict(
                    id=identifier, query_start=low, query_end=high, orientation=orientation,
                    amino_acids=protein, coding_sequence=predicted.coding_sequence,
                    complete_start=complete_start, observed_stop=predicted.stop_complete,
                    status="complete_CDS" if complete_start and predicted.stop_complete else "partial_CDS",
                    translation_table=model.genetic_code.table_id,
                    evidence="reference_anchored_prediction_from_observed_RNA", reference_orfs=[]))
                coding_positions = observed_positions[query:end]
                association = dict(
                    protein_id=model.orf_id, gene_name=model.gene_name, source=model.source,
                    reference_context=model.context.to_dict(), reference_key=model.context.key,
                    coding_blocks=[list(block) for block in model.coding_blocks],
                    reference_amino_acids=model.protein_sequence,
                    protein_matches_reference=(protein == model.protein_sequence)
                    if complete_start and predicted.stop_complete else None,
                    anchor=dict(query_offset=query, reference_coding_offset=offset,
                                coordinates="0-based-oriented-RNA"),
                    coding_path=_reference_path_status(coding_positions, reference_positions, offsets),
                    start_status="observed_source_start" if complete_start else
                    "source_start_altered" if offset == 0 else "source_start_uncovered",
                    stop_status=("source_stop" if coding_positions[-3:] == reference_positions[-3:] else
                                 "alternative_stop") if predicted.stop_complete else "stop_uncovered_or_lost")
                row["reference_orfs"].append(association)
                break
    return sorted(readings.values(), key=lambda row: (row["query_start"], row["query_end"],
                                                      row["orientation"], row["amino_acids"]))


def write_viral_proteins(result, path):
    """Export complete and partial predicted coding products with explicit status.

    Parameters
    ----------
    result : dict
        Result from ``reconstruct_viral_sequences``.
    path : str or pathlib.Path
        Output protein FASTA, empty when no coding products are available.
    """
    with Path(path).open("w") as handle:
        for sequence in result["sequences"]:
            for coding in sequence["coding"]:
                if not coding["amino_acids"]:
                    continue
                sources = ",".join(sorted({r["protein_id"] for r in coding["reference_orfs"]}))
                handle.write(">%s sequence_id=%s status=%s orientation=%s reference_ORF_anchors=%s\n%s\n" % (
                    coding["id"], sequence["id"], coding["status"], coding["orientation"], sources,
                    coding["amino_acids"]))
