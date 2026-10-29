"""Variant-free, targeted RNA reconstruction with conservative attribution.

Uses the same observed placements, exact overlaps and original segment
identities as SV reconstruction. Reference models are competitors, never
sequence donors. Transcript abundance is estimated separately by Salmon.
"""

from collections import defaultdict
from dataclasses import asdict
from hashlib import sha256
from tempfile import TemporaryDirectory


from .antigen_reference import IsoformTarget
from .genetic_code import standard_genetic_code
from .read_collector import ReadCollector
from .sv_rna import _Model, _runs
from .rna_reconstruction import _contains, _discover, _support
from .rna_library import _original_library_context


def _coding(sequence, positions, models):
    readings = {}
    for model in models:
        reference = model.reference
        if reference.cds_start is None:
            continue
        # Transfer a frame through an observed, collinear complete codon.
        # Never choose a longest ORF or fill an uncovered start with reference.
        for q in range(len(sequence) - 2):
            offset = model.offsets.get(positions[q])
            if (offset is None or not reference.cds_start <= offset < reference.cds_end - 3
                    or (offset - reference.cds_start) % 3
                    or list(positions[q:q + 3]) != model.positions[offset:offset + 3]):
                continue
            complete_start = offset == reference.cds_start and sequence[q:q + 3] in standard_genetic_code.start_codons
            protein, stop = standard_genetic_code.translate(sequence[q:], first_codon_is_start=complete_start)
            if not protein:
                break
            key = (q, protein, complete_start, stop)
            row = readings.setdefault(key, dict(
                query_start=q, query_end=q + 3 * (len(protein) + stop), amino_acids=protein,
                complete_start=complete_start, observed_stop=stop,
                status="complete_CDS" if complete_start and stop else "partial_CDS", transcript_ids=[]))
            row["transcript_ids"].append(reference.transcript_id)
            break
    return sorted(readings.values(), key=lambda row: (row["query_start"], row["amino_acids"]))


def _candidate(sequence, positions, targets, models, annotated=()):
    positions = tuple(positions)
    genes = sorted({t.gene_id for t in targets if any(
        p is not None and p[0] == t.reference.contig and p[2] == t.reference.strand
        and t.reference.exons[0][0] <= p[1] < t.reference.exons[-1][1] for p in positions)})
    if annotated:
        genes = sorted({t.gene_id for t in targets if t.reference.transcript_id in annotated})
    key = repr((sequence, positions)).encode()
    reaches_ends = any(positions[0] == model.positions[0] and positions[-1] == model.positions[-1]
                       for model in models)
    return dict(id="isovar_isoform_" + sha256(key).hexdigest()[:24], sequence=sequence,
                placements=positions, gene_ids=genes, annotated_transcript_ids=sorted(annotated),
                origin="reference_model" if annotated else "observed_assembly",
                reaches_annotated_transcript_ends=reaches_ends,
                coding=_coding(sequence, positions, models))


def _protein_groups(candidates, support):
    groups = defaultdict(lambda: dict(isoform_ids=[], fully_covered_isoform_ids=[], whole_cds_segments=set()))
    for row in candidates:
        products = {c["amino_acids"] for c in row["coding"] if c["status"] == "complete_CDS"}
        for protein in sorted(products):
            group = groups[protein]
            group["isoform_ids"].append(row["id"])
            matching = [c for c in row["coding"] if c["amino_acids"] == protein and c["status"] == "complete_CDS"]
            if any(c["all_bases_observed"] for c in matching):
                group["fully_covered_isoform_ids"].append(row["id"])
            group["whole_cds_segments"].update(tuple(key) for c in matching for key in c["whole_cds_segments"])
    output = []
    for protein, group in sorted(groups.items()):
        ids = set(group["isoform_ids"])
        group["whole_cds_segments"] = [list(key) for key in sorted(group["whole_cds_segments"])]
        group["compatible_fragments"] = sum(c["fragments"] for c in support["equivalence_classes"] if ids.intersection(c["isoform_ids"]))
        output.append(dict(amino_acids=protein, **group))
    return output


def reconstruct_antigen_isoforms(bam, *, targets, sample_id, source, reference_identity=None,
                                max_discovery_reads=10000, max_paths=256, min_overlap=30,
                                scratch_dir=None):
    """Reconstruct targeted RNA without a nominated mutation.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Complete declared input; BAM, SAM and CRAM use the shared original-template index.
    targets : iterable of IsoformTarget
        All annotated competitors of the explicitly selected genes. Include
        related loci when shared reads could otherwise be misattributed.
    sample_id, source : str
        Explicit namespaces retained alongside original RG/QNAME identities.
    reference_identity : dict, optional
        Alignment/annotation identity and file checksums supplied by the caller.
    max_discovery_reads, max_paths, min_overlap : int
        Discovery resource limits and minimum exact overlap. All eligible
        records are subsequently checked for support, including omitted reads.

    Returns
    -------
    dict
        Versioned sequence/evidence result. Compatible fragments are not TPM or
        independent molecules. Noisy reads may support structure without an
        exact sequence match; this initial workflow leaves them unassigned.
    """
    targets = tuple(targets)
    if not targets or not all(isinstance(t, IsoformTarget) for t in targets):
        raise ValueError("Provide at least one IsoformTarget")
    targets = tuple(sorted(targets, key=lambda t: t.reference.transcript_id))
    if not sample_id or not source:
        raise ValueError("Explicit sample_id and source are required")
    samples = {row["SM"] for row in bam.header.to_dict().get("RG", []) if row.get("SM")}
    if len(samples) > 1:
        raise ValueError("The BAM declares multiple samples; reconstruct each sample separately")
    if len({t.reference.reference_name for t in targets}) != 1:
        raise ValueError("Targets must share one alignment assembly")
    for name, value in (("max_discovery_reads", max_discovery_reads), ("max_paths", max_paths), ("min_overlap", min_overlap)):
        if type(value) is not int or value < 1:
            raise ValueError("%s must be a positive integer" % name)
    names = [t.reference.transcript_id for t in targets]
    if len(set(names)) != len(names):
        raise ValueError("Reference transcript IDs must be unique")
    lengths = dict(zip(bam.references, bam.lengths))
    for target in targets:
        r = target.reference
        if r.contig not in lengths or r.exons[-1][1] > lengths[r.contig]:
            raise ValueError("Target transcript lies outside the BAM reference: %s" % r.transcript_id)
    collector = ReadCollector(use_secondary_alignments=True, min_mapping_quality=0,
                              use_soft_clipped_bases=False, use_reads_without_base_qualities=True)
    regions = defaultdict(set)
    for t in targets:
        r = t.reference
        regions[r.contig].add((r.exons[0][0], r.exons[-1][1], r.strand))
    models = [_Model(t.reference) for t in targets]
    with TemporaryDirectory(dir=scratch_dir, prefix="isovar-isoforms-") as directory, _original_library_context(
            bam, collector, directory, library_index=None, reference_name=targets[0].reference.reference_name, read_filter_identity=None) as groups:
        statistics = groups.statistics
        paths, discovery = _discover(groups, regions, collector, max_discovery_reads, max_paths, min_overlap)
        candidates = {}
        for target, model in zip(targets, models):
            r = target.reference
            key = r.sequence, tuple(model.positions)
            if key in candidates:
                candidates[key]["annotated_transcript_ids"].append(r.transcript_id)
                candidates[key]["gene_ids"] = sorted(set(candidates[key]["gene_ids"]) | {target.gene_id})
            else:
                candidates[key] = _candidate(*key, targets, models, [r.transcript_id])
        partial = []
        for sequence, positions in paths:
            key = sequence, positions
            if key in candidates:
                continue
            row = _candidate(*key, targets, models)
            containers = [c["id"] for k, c in candidates.items() if _contains(k, key)]
            if containers or not row["reaches_annotated_transcript_ends"]:
                row["compatible_models"] = sorted(containers)
                partial.append(row)
            else:
                candidates[key] = row
        candidates = sorted(candidates.values(), key=lambda c: c["id"])
        support = _support(groups, candidates + partial, regions, collector, bam.header.to_dict(), sample_id, source)
        full_ids = {c["id"] for c in candidates}
        support["fragments_without_full_transcript_assignment"] = sum(
            row["fragments"] for row in support["equivalence_classes"] if not full_ids.intersection(row["isoform_ids"]))
    for row in candidates + partial:
        # Compact placed runs leave insertions explicitly unplaced; query
        # sequence length remains available to reconstruct the placement map.
        row["mapping"] = [list(run) for run in _runs(row.pop("placements"))]
    return dict(schema="isovar.antigen_isoforms.v1", sample_id=sample_id, source=source,
                reference=reference_identity or {}, targets=[asdict(t) for t in targets],
                input=dict(records=statistics['records_scanned'], unique_records=statistics['unique_records'],
                           eligible_segments=statistics['eligible_segments'], excluded_records=statistics['excluded_records'], header=bam.header.to_dict()),
                discovery=discovery, isoforms=candidates, partial_sequences=partial, support=support,
                protein_groups=_protein_groups(candidates + partial, support), quantification=None,
                limitations=["exact_overlap_reconstruction", "RNA_does_not_establish_translation",
                             "fragment_counts_are_not_molecule_counts",
                             "unreported_alignments_not_recovered",
                             "overlap_assembly_does_not_establish_end_to_end_phase"])
