"""Observed viral-associated RNA paths without a mutation or transcript GTF."""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Iterable
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import pysam

from pyensembl import Genome

from .antigen_reference import GDC_REFERENCE, gdc_viral_panel, validate_alignment_reference
from .audit_bundle import digest
from .read_collector import ReadCollector
from .rna_reconstruction import _assemble, _observations, _support
from .sv_rna import _GeometryCache, _collection_parameters, _runs
from .rna_library import RnaLibraryIndex, _original_library_context
from .viral_coding import viral_coding_predictions, viral_coding_references


def _discover_viruses(groups, regions, collector, max_discovery_reads, max_paths, min_overlap):
    # The shared placement index selects complete templates per contig without
    # rescanning the original library or depending on its physical record order.
    candidates, discovery, reasons = {}, [], Counter()
    for contig in sorted(regions):
        chosen, considered, selected, cache = {}, 0, 0, _GeometryCache()
        for template in groups.templates([(contig, start, end) for start, end, _ in regions[contig]]):
            for _, observations in _observations(template, collector, reasons, regions, cache):
                observed = [o for o in observations if any(p is not None and p[0] == contig for p in o.positions)]
                if not observed:
                    continue
                considered += 1
                if selected >= max_discovery_reads:
                    continue
                selected += 1
                orientations = {}
                for o in observed:
                    key = o.key[:-1]
                    previous = orientations.get(key)
                    if previous is None or (tuple(_runs(o.positions)), o.sequence) < (
                            tuple(_runs(previous.positions)), previous.sequence):
                        orientations[key] = o
                for o in orientations.values():
                    chosen.setdefault((o.sequence, o.positions), o)
        paths, notes = _assemble(chosen, max_paths, min_overlap)
        discovery.append(dict(contig=contig, eligible_observed_segments=considered,
                              selected_segments=selected, distinct_observations=len(chosen),
                              paths=len(paths), truncated=considered > selected or bool(notes),
                              limitations=sorted(notes)))
        for sequence, positions in paths:
            key = sequence, positions
            names = sorted({p[0] for p in positions if p is not None})
            candidates.setdefault(key, dict(
                id="isovar_viral_" + sha256(repr(key).encode()).hexdigest()[:24],
                sequence=sequence, placements=positions, coding=[], origin="observed_assembly",
                viral_contigs=sorted(set(names).intersection(regions)),
                other_contigs=sorted(set(names).difference(regions)),
                rna_strand="unknown", transcript_completeness="unknown",
                coding_status="unannotated", orientation="canonical_reference_orientation"))
    return sorted(candidates.values(), key=lambda row: row["id"]), dict(
        scope="per_viral_contig", template_order="read_group_query_name", max_discovery_reads=max_discovery_reads, max_paths=max_paths,
        peak_discovery_observations=max((r["distinct_observations"] for r in discovery), default=0),
        min_overlap=min_overlap, contigs=discovery, truncated=any(r["truncated"] for r in discovery),
        observation_rejections=dict(reasons))


def _contig_support(contigs, candidates, support, discovery):
    counts = {name: dict(aligned_fragments=0, aligned_segments=0, unambiguous_contig_fragments=0,
                         sequence_assigned_fragments=0) for name in contigs}
    by_id = {row["id"]: row for row in candidates}
    for witness in support["witnesses"]:
        segments = defaultdict(set)
        for record in witness["records"]:
            if record["contig"] in counts and not record["sam_flag"] & 4:
                segments[record["contig"]].add(record["sam_flag"] & 192 if record["sam_flag"] & 1 else 0)
        assigned = {name for cid in witness["compatible_sequences"] for name in by_id[cid]["viral_contigs"]}
        for name, identities in segments.items():
            row = counts[name]
            row["aligned_fragments"] += 1
            row["aligned_segments"] += len(identities)
            row["unambiguous_contig_fragments"] += (
                len(segments) == 1 and not witness["competing_placement_outside_targets"])
            row["sequence_assigned_fragments"] += name in assigned
    searches = {row["contig"]: row for row in discovery["contigs"]}
    return [dict(contig=name, **counts[name],
                 unassigned_fragments=counts[name]["aligned_fragments"] - counts[name]["sequence_assigned_fragments"],
                 sequence_ids=[c["id"] for c in candidates if name in c["viral_contigs"]],
                 status=("no_eligible_alignments" if not counts[name]["aligned_fragments"] else
                         "discovery_truncated" if searches[name]["truncated"] else
                         "sequences_reconstructed" if any(name in c["viral_contigs"] for c in candidates) else
                         "aligned_but_not_reconstructed")) for name in sorted(contigs)]


def reconstruct_viral_sequences(
        bam: pysam.AlignmentFile, *, genome_fasta: str | Path, sample_id: str, source: str,
        contigs: Iterable[str] | None = None, reference_profile: str | None = GDC_REFERENCE,
        max_discovery_reads: int = 10000, max_paths: int = 256, min_overlap: int = 30,
        scratch_dir: str | Path | None = None, reference_cache: str | Path | None = None,
        read_collector: ReadCollector | None = None, library_index: RnaLibraryIndex | None = None,
        read_filter_identity: str | None = None) -> dict[str, Any]:
    """Recover viral-associated nucleotide paths from the complete declared input.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        RNA aligned competitively to the supplied human/viral reference. All
        available mates and alternative placements should remain in the input.
    genome_fasta : str or Path
        Local alignment reference, used for identity checks, never gap filling.
    sample_id, source : str
        Sample and input-library identities retained with original witnesses.
    contigs : iterable of str, optional
        Viral contig names. Defaults to GDC's complete published viral panel.
        Must be supplied explicitly for a custom reference.
    reference_profile : str or None
        GDC_REFERENCE checks the exact combined dictionary and viral DNA MD5s.
        None selects a custom reference; missing BAM M5 remains unverified.
    max_discovery_reads, max_paths, min_overlap : int
        Per-contig discovery limits and minimum exact overlap. The subsequent
        support pass always considers the full declared input.
    scratch_dir, reference_cache : str or Path, optional
        Temporary disk-index location and PyEnsembl genomic-FASTA cache root.
    read_collector : ReadCollector, optional
        Shared record/sequence policy for discovery and support; defaults to
        ReadCollector.for_rna_paths(). Original mate identities remain separate;
        the small-variant merge_overlapping_fragments override is unsupported.
    library_index : RnaLibraryIndex, optional
        Caller-owned index of the same BAM/reference/eligibility policy. Reuse
        avoids rescanning the full library for each reconstruction task.
    read_filter_identity : str, optional
        Stable custom-filter identity required when reusing a library index.

    Returns
    -------
    dict
        Observed sequences, mappings, fragment support and original witnesses.
        The canonical sequence orientation does not establish RNA strand or
        transcript ends. No unannotated ORFs, protein abundance or viral absence
        are inferred. Exact overlaps retain both alleles and sequencing errors.
        Available sequence-qualified ORFdata CDSs anchor coding predictions;
        these do not establish initiation, actual RNA strand or protein expression.
    """
    if not sample_id or not source:
        raise ValueError("Explicit sample_id and source are required")
    samples = {r["SM"] for r in bam.header.to_dict().get("RG", []) if r.get("SM")}
    if len(samples) > 1:
        raise ValueError("The BAM declares multiple samples; reconstruct each sample separately")
    for name, value in (("max_discovery_reads", max_discovery_reads), ("max_paths", max_paths), ("min_overlap", min_overlap)):
        if type(value) is not int or value < 1:
            raise ValueError("%s must be a positive integer" % name)
    panel = gdc_viral_panel() if reference_profile == GDC_REFERENCE else None
    if reference_profile not in (None, GDC_REFERENCE):
        raise ValueError("Unknown alignment reference profile")
    if contigs is None:
        if panel is None:
            raise ValueError("Declare viral contigs explicitly for a custom reference")
        contigs = [r["name"] for r in panel["contigs"]]
    if isinstance(contigs, str):
        raise ValueError("contigs must be a collection of names, not one string")
    contigs = list(contigs)
    if not contigs or any(not isinstance(c, str) or not c for c in contigs):
        raise ValueError("Select at least one viral contig")
    contigs = sorted(set(contigs))
    if panel and set(contigs).difference(r["name"] for r in panel["contigs"]):
        raise ValueError("Selected contig is not in GDC's published viral panel")
    lengths = dict(zip(bam.references, bam.lengths))
    if set(contigs).difference(lengths):
        raise ValueError("Selected viral contigs are missing from the BAM; competitive realignment is required")
    fasta_path = Path(genome_fasta).expanduser().resolve()
    if not fasta_path.is_file():
        raise ValueError("A local alignment FASTA is required")
    with Genome(reference_name="GRCh38" if panel else "custom", annotation_name="unannotated",
                genome_fasta_path_or_url=str(fasta_path), cache_directory_path=reference_cache) as genome:
        identity = validate_alignment_reference(bam, genome, profile=reference_profile, contigs=contigs)
        coding_models, coding_reference = viral_coding_references(genome, contigs)
    identity["files"] = dict(genome_fasta=dict(path=str(fasta_path), sha256=digest(fasta_path)))
    if panel:
        identity["viral_panel"] = panel
    regions = {name: [(0, lengths[name], "+")] for name in contigs}
    collector = read_collector if read_collector is not None else ReadCollector.for_rna_paths()
    if not collector.merge_overlapping_fragments:
        raise ValueError('merge_overlapping_fragments is not a control for original RNA paths')
    with TemporaryDirectory(dir=scratch_dir, prefix="isovar-viral-") as directory, _original_library_context(
            bam, collector, directory, library_index=library_index,
            reference_name="GRCh38" if panel else None, read_filter_identity=read_filter_identity) as groups:
        statistics = groups.statistics
        candidates, discovery = _discover_viruses(
            groups, regions, collector, max_discovery_reads, max_paths, min_overlap)
        covered_contigs = {r.transcript.contig for r in coding_models}
        for candidate in candidates:
            candidate["coding"] = viral_coding_predictions(
                candidate["id"], candidate["sequence"], candidate["placements"], coding_models)
            candidate["coding_status"] = (
                "reference_anchored_prediction" if candidate["coding"] else
                "no_coding_overlap" if set(candidate["viral_contigs"]).intersection(covered_contigs)
                else "unavailable_reference")
        support = _support(groups, candidates, regions, collector, bam.header.to_dict(),
                           sample_id, source, label="sequence")
        counts = _contig_support(contigs, candidates, support, discovery)
    for row in candidates:
        row["mapping"] = [list(run) for run in _runs(row.pop("placements"))]
    return dict(schema="isovar.viral_sequences.v1", sample_id=sample_id, source=source,
                collection=dict(_collection_parameters(collector), read_filter_identity=read_filter_identity),
                reference=identity, coding_reference=coding_reference, selected_contigs=contigs,
                input=dict(records=statistics['records_scanned'], unique_records=statistics['unique_records'],
                           eligible_segments=statistics['eligible_segments'], excluded_records=statistics['excluded_records'], header=bam.header.to_dict()),
                discovery=discovery, sequences=candidates, contigs=counts, support=support,
                limitations=["exact_overlap_reconstruction", "RNA_strand_and_transcript_ends_unknown",
                             "coding_reference_is_a_source_CDS_subset", "protein_predictions_are_not_expression_evidence",
                             "no_viral_absence_or_tumor_origin_claim",
                             "fragment_counts_are_not_molecule_counts", "unreported_alignments_not_recovered",
                             "overlap_assembly_does_not_establish_end_to_end_phase"])


def write_viral_fasta(result, path):
    """Export only observed sequence paths, including explicitly partial paths."""
    with Path(path).open("w") as out:
        for row in result["sequences"]:
            out.write(">%s viral_contigs=%s rna_strand=unknown transcript_completeness=unknown\n%s\n" % (
                row["id"], ",".join(row["viral_contigs"]), row["sequence"]))
