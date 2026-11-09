"""Measure local RNA features through the existing original-path machinery."""

from collections import Counter, defaultdict
from tempfile import TemporaryDirectory

from .read_collector import ReadCollector
from .rna_feature_support import RNAFeature, measure_rna_feature
from .rna_library import original_library_context
from .rna_path_linkage import OriginalPathMatch
from .rna_reconstruction import original_template_observations, placement_uncertainty
from .rna_observation import ObservationGeometryCache
from .rna_path_geometry import PathPositions, position_runs, overlap_offsets
from .rna_evidence import content_identifier



def _local_matches(sequence, positions, observation, identity):
    """Exact runs of one original placement, preserving alternative path offsets."""
    offsets = sorted({offset for _, offset in overlap_offsets(
        position_runs(positions), 0, observation.key, position_runs(observation.positions))})
    for offset in offsets:
        low, high = max(0, offset), min(len(sequence), offset + len(observation.sequence))
        key = content_identifier("feature_placement", [observation.key, offset])
        start = None
        for q in range(low, high):
            i = q - offset
            exact = (sequence[q] in "ACGT" and sequence[q] == observation.sequence[i]
                     and positions[q] == observation.positions[i])
            if exact and start is None:
                start = q
            if not exact and start is not None:
                yield OriginalPathMatch(identity, key, start, q)
                start = None
        if start is not None:
            yield OriginalPathMatch(identity, key, start, high)

def measure_rna_path_features(
        bam, *, sequence, positions, features, sample_id, source,
        profile=None, read_collector=None, library_index=None,
        read_filter_identity=None, scratch_dir=None):
    """Count a batch of features on one supplied nucleotide path from original RNA.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Declared input alignment. A regional input remains regional; its record
        count never supplies a whole-library denominator automatically.
    sequence : str
        Nucleotide path to assess, supplied by reconstruction or an explicitly
        identified reference/alternate model. Supplying a model does not assert
        that all its bases were observed.
    positions : sequence of tuple or None
        Per-base ``(contig, zero_based_position, strand)`` placements in path
        orientation. Inserted bases use None. Existing coordinate-aware exact
        matching is reused; no new BAM/CIGAR parser is introduced.
    features : iterable of RNAFeature
        Local allele, junction, coding or compound features on this same path.
        Complete and partial support are separate feature definitions.
    sample_id, source : str
        Source-qualified evidence identity.
    profile : RNASamplingProfile, optional
        Independently supplied library geometry and denominators.
    read_collector : ReadCollector, optional
        Existing shared original-record/sequence eligibility policy.
    library_index : RnaLibraryIndex, optional
        Reuse a prepared index for multiple targets in the same input.
    read_filter_identity : str, optional
        Required identity for a custom filter with a reused index.
    scratch_dir : path, optional
        Parent of the temporary index directory when no index is supplied.

    Returns
    -------
    list of dict
        Validated, JSON-compatible common feature measurements. Attribution is
        compatible to this supplied coordinate path; other transcript/locus
        assignments are not silently declared absent. Known original-placement
        uncertainty is retained with the eligibility and source identity.
    """
    features = tuple(features)
    if not features or any(not isinstance(f, RNAFeature) for f in features):
        raise ValueError("At least one RNAFeature is required")
    if len({f.identifier for f in features}) != len(features) or len({f.path_id for f in features}) != 1:
        raise ValueError("Features require distinct identities on the same RNA path")
    sequence = sequence.upper()
    positions = PathPositions(None if p is None else tuple(p) for p in positions)
    if any(p is not None and (len(p) != 3 or not isinstance(p[0], str) or not p[0]
            or type(p[1]) is not int or p[1] < 0 or p[2] not in ("+", "-")) for p in positions):
        raise ValueError("RNA path positions require contig, nonnegative position and strand")
    if len(sequence) != len(positions) or not sequence or set(sequence) - set("ACGTN"):
        raise ValueError("RNA path requires matching nucleotide sequence and placements")
    if any(max(b for _, b in f.intervals) > len(sequence) or
           f.path_length not in (None, len(sequence)) for f in features):
        raise ValueError("Feature coordinates disagree with the supplied RNA path")
    if any(set(sequence[a:b]) - set("ACGT") for f in features for a, b in f.intervals):
        raise ValueError("Required feature bases must be resolved; unknown bases are not measured absence")
    regions = defaultdict(list)
    for _, _, contig, start, end, strand in position_runs(positions):
        regions[contig].append((start, end, strand))
    if not regions:
        raise ValueError("RNA path needs coordinate anchors for original-read matching")
    lengths = dict(zip(bam.references, bam.lengths))
    if any(contig not in lengths or end > lengths[contig]
           for contig, spans in regions.items() for _, end, _ in spans):
        raise ValueError("RNA path placements are outside the supplied alignment reference")
    collector = read_collector or ReadCollector.for_rna_paths()
    matches, rejections, uncertainty = [], Counter(), Counter()
    cache = ObservationGeometryCache()
    with TemporaryDirectory(prefix="isovar-feature-", dir=scratch_dir) as directory:
        with original_library_context(bam, collector, directory, library_index=library_index,
                reference_name=None, read_filter_identity=read_filter_identity) as groups:
            for template in groups.templates([(contig, a, b) for contig, spans in regions.items()
                                               for a, b, _ in spans]):
                uncertainty.update(placement_uncertainty(template))
                for identity, observations in original_template_observations(template, collector, rejections, regions, cache):
                    for observation in observations:
                        matches.extend(_local_matches(sequence, positions, observation, identity))
            provenance = dict(input_identity=groups.prepared.identity,
                              input_extent="declared_input",
                              attribution="compatible_coordinate_path",
                              matching_scope="declared_feature_intervals",
                              observation_rejections=dict(rejections),
                              placement_uncertainty=dict(uncertainty),
                              read_policy=groups.prepared.identity["eligibility"])
    return [measure_rna_feature(f, matches, evidence_scope=(sample_id, source),
                                profile=profile, provenance=provenance) for f in features]
