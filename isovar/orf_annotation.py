"""Overlay shared reference ORF evidence on reconstructed SV RNA candidates."""

from copy import deepcopy

from varcode.effect_hypotheses import RealizedTranscriptProduct
from varcode.orf import annotate_orf_product

from orfdata import load_catalog


def annotate_sv_orf_evidence(result, transcripts, *, catalog=None):
    """Return a copy with per-candidate reference-ORF annotations.

    ``transcripts`` supplies full reference transcript objects for validating
    each catalogue model. No annotation download or strongest-isoform selection
    is implicit. Isovar path blocks are zero-based half-open; the Varcode
    interchange uses one-based base origins. Each candidate keeps its own
    complete-interval RNA support. Reference MS never proves mutant translation.
    """
    catalog = load_catalog() if catalog is None else catalog
    if result["reference_name"] != catalog.genome_build:
        raise ValueError("RNA result and ORF catalogue assemblies differ")
    by_id = {str(t.id).split(".")[0]: t for t in transcripts}
    annotated = deepcopy(result)
    for path in annotated.get("paths", ()):
        sequence = path["sequence"]
        origins = [None] * len(sequence)
        for block in path["blocks"]:
            start, end = block["query_start"], block["query_end"]
            a, b = block["reference_start"], block["reference_end"]
            if not 0 <= start <= end <= len(sequence) or end - start != b - a:
                raise ValueError("Invalid reconstructed RNA block")
            positions = range(a + 1, b + 1) if block["strand"] == "+" else range(b, a, -1)
            for i, position in zip(range(start, end), positions):
                if origins[i] is not None:
                    raise ValueError("Overlapping reconstructed RNA blocks")
                origins[i] = (block["contig"], position, block["strand"])
        for candidate in path.get("exploratory_orfs", {}).get("candidates", ()):
            start, end = candidate["query_interval"]
            origin = origins[start]
            annotations = []
            if origin is not None:
                contig, position, strand = origin
                for record in catalog.query(contig=contig, start=position, end=position,
                                            genome_build=result["reference_name"]):
                    mapped = record.map_genomic_position(position, genome_build=result["reference_name"])
                    if record.strand != strand or mapped["coding_nucleotide_offset"] != 0:
                        continue
                    transcript = by_id.get(record.transcript_id)
                    if transcript is None:
                        annotations.append(dict(orf_id=record.orf_id, status="unresolved",
                                                reason="reference_transcript_not_supplied"))
                        continue
                    # Candidate-specific evidence, not the path/junction total.
                    product = RealizedTranscriptProduct(
                        transcript, sequence, None,
                        evidence=dict(base_origins=origins, path_id=path["path_id"],
                                      event_id=result["event_id"], source=result.get("source"),
                                      sample_id=result.get("sample_id"),
                                      candidate_query_interval=[start, end],
                                      rna_support=candidate["full_interval_support"],
                                      rna_polarity="unresolved",
                                      initiation_observed=candidate["initiation_observed"],
                                      translation_observed=candidate["translation_observed"]))
                    try:
                        annotation = annotate_orf_product(record.bind(transcript), product)
                    except ValueError as error:
                        annotation = dict(orf_id=record.orf_id, status="unresolved", reason=str(error))
                    if (annotation.get("status") == "conditional_translation"
                            and annotation["product"]["protein_sequence"] != candidate["amino_acids"]):
                        raise ValueError("Mapped ORF translation disagrees with RNA candidate")
                    annotations.append(annotation)
            candidate["reference_orf_annotations"] = annotations
            candidate["reference_orf_lookup"] = "matched" if annotations else "no_matching_reference_start"
    annotated["orf_catalogue_metadata"] = catalog.metadata
    return annotated
