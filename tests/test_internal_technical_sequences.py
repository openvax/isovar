"""Simulated internal linkers: preserve originals and never synthesize a join.

These are matching/eligibility controls, not observed PacBio adapter positives.
The unchanged real Sid terminal fixtures remain in test_read_end_inference.
"""

from collections import Counter
from dataclasses import replace
from itertools import product

import pysam
import pytest

from isovar import Adapter, ReadCollector, ReadEndProfile, infer_read_ends, read_sequence_view_from_alignment
from isovar.dna import reverse_complement_dna
from isovar.read_end_inference import ReadEndInferenceCache, adapter_matches
from isovar.read_sequence_origin import GenomicQueryEvidence, genomic_query_evidence_from_alignment
from isovar.rna_library import prepare_rna_library
from isovar.sv_path_routing import PathEvidenceQuery, route_path_evidence
from isovar.sv_rna import _build_segment, _record_gaps
from isovar.sv_rna_assessment import aligned_exon_clip_contains, assess_sv_rna_paths, evidence_from_records
from isovar.sv_rna_path_support import exact_path_matches, interval_spots, overlap_chain_search
from tests.test_sv_rna import Scenario, link, record, write_bam


LEFT = 'GCTACGTCGCTG'
RIGHT = 'CTGACGTGACCG'
MOTIF = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
ADAPTER = Adapter('simulated-linker', MOTIF, role='linker', placement='internal',
                  source='synthetic regression sequence')
PROFILE = ReadEndProfile('simulated-concatemer', adapters=(ADAPTER,))
SEQUENCE = LEFT + MOTIF + RIGHT


def collector():
    return ReadCollector(read_end_profile=PROFILE, trim_adapters=True, use_soft_clipped_bases=True)


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('observed', [MOTIF, MOTIF[:8] + 'T' + MOTIF[9:],
                                      MOTIF[:8] + 'C' + MOTIF[8:], MOTIF[:8] + MOTIF[9:]])
def test_internal_full_motif_exposes_separate_original_intervals(observed, reverse):
    sequence = LEFT + observed + RIGHT
    if reverse:
        sequence = reverse_complement_dna(sequence)
    view = infer_read_ends(sequence, range(len(sequence)), PROFILE, reverse=reverse, trim_adapters=True)
    assert view.sequence == sequence  # Existing contiguous sequence semantics are unchanged.
    assert view.original_qualities == tuple(range(len(sequence)))
    assert len(view.usable_intervals) == 2
    assert all(annotation.kind == 'linker' and annotation.side == 'internal' for annotation in view.annotations)
    assert all(annotation.origin == 'reference_unavailable' for annotation in view.annotations)
    assert not view.contains_usable_interval(0, len(sequence))
    assert [sequence[a:b] for a, b in view.usable_intervals] == (
        [reverse_complement_dna(RIGHT), reverse_complement_dna(LEFT)] if reverse else [LEFT, RIGHT])
    unchanged = infer_read_ends(sequence, profile=PROFILE, reverse=reverse)
    assert unchanged.usable_intervals == ((0, len(sequence)),)


def test_distant_noisy_copy_is_not_hidden_by_a_perfect_copy():
    noisy = MOTIF[:8] + 'T' + MOTIF[9:]
    sequence = LEFT + MOTIF + RIGHT + noisy + LEFT
    view = infer_read_ends(sequence, profile=PROFILE, trim_adapters=True)
    assert len(view.annotations) == 2
    assert {annotation.errors for annotation in view.annotations} == {0, 1}
    assert [sequence[a:b] for a, b in view.usable_intervals] == [LEFT, RIGHT, LEFT]


def _distance(first, second):
    row = list(range(len(second) + 1))
    for index, a in enumerate(first, 1):
        next_row = [index]
        for j, b in enumerate(second, 1):
            next_row.append(min(row[j] + 1, next_row[-1] + 1, row[j - 1] + (a != b)))
        row = next_row
    return row[-1]


def test_internal_boundaries_match_an_independent_edit_distance_oracle():
    adapter = Adapter('small', 'ACGT', min_overlap=4, max_error_rate=0.3,
                      placement='internal', min_internal_flank=1)
    for letters in product('ACGT', repeat=5):
        sequence = 'G' + ''.join(letters) + 'C'
        scores = {(a, b): _distance(adapter.sequence, sequence[a:b])
                  for a in range(1, 6) for b in range(a + 1, 7)}
        best = min(scores.values())
        expected = {span for span, distance in scores.items() if distance == best and distance <= 1}
        observed = {(hit.start, hit.end) for hit in adapter_matches(sequence, adapter, 4)}
        assert observed == expected, sequence


def test_internal_partials_and_short_flanks_do_not_establish_a_join():
    for sequence in (LEFT + MOTIF[:15] + RIGHT, MOTIF + RIGHT, LEFT + MOTIF):
        assert not adapter_matches(sequence, ADAPTER, 20)


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('equality', ['MD', '='])
def test_adapter_like_genomic_sequence_remains_usable(reverse, equality):
    sequence = reverse_complement_dna(SEQUENCE) if reverse else SEQUENCE
    read = record('genomic', '1', 1000, f'{len(sequence)}' + ('=' if equality == '=' else 'M'),
                  sequence, flag=16 if reverse else 0)
    if equality == 'MD':
        read.set_tag('MD', str(len(sequence)))
    original = read.to_string()
    view = read_sequence_view_from_alignment(read, PROFILE, trim_adapters=True)
    assert view.annotations[0].origin == 'reference_supported'
    assert view.barrier_intervals == ()
    assert view.usable_intervals == ((0, len(sequence)),)
    assert read.to_string() == original


def test_genomic_competitor_is_fresh_for_each_placement_even_with_cached_sequence():
    read = record('same', '1', 1000, f'{len(SEQUENCE)}M', SEQUENCE)
    cache = ReadEndInferenceCache()
    unresolved = read_sequence_view_from_alignment(read, PROFILE, trim_adapters=True, cache=cache)
    read.set_tag('MD', str(len(SEQUENCE)))
    resolved = read_sequence_view_from_alignment(read, PROFILE, trim_adapters=True, cache=cache)
    assert unresolved.barrier_intervals and not resolved.barrier_intervals
    assert unresolved.annotations[0].origin == 'alignment_competitor'
    assert resolved.annotations[0].origin == 'reference_supported'
    assert (cache.hits, cache.misses) == (1, 1)


@pytest.mark.parametrize('gap', ['I', 'D', 'N'])
def test_cigar_insertions_and_reference_gaps_do_not_claim_contiguous_genomic_identity(gap):
    if gap == 'I':
        cigar, md = f'12M{len(MOTIF)}I12M', '24'
    else:
        cigar = f'20M30{gap}{len(SEQUENCE) - 20}M'
        md = '20^' + 'A' * 30 + str(len(SEQUENCE) - 20) if gap == 'D' else str(len(SEQUENCE))
    read = record('join', '1', 1000, cigar, SEQUENCE)
    read.set_tag('MD', md)
    before = read.to_string()
    view = read_sequence_view_from_alignment(read, PROFILE, trim_adapters=True)
    assert view.annotations[0].origin == 'alignment_competitor'
    assert view.barrier_intervals
    assert read.to_string() == before


def test_explicit_reference_evidence_is_length_checked_and_source_identified():
    evidence = GenomicQueryEvidence(exact_intervals=((12, 12 + len(MOTIF)),), sources=('reference:checksum:interval',))
    view = infer_read_ends(SEQUENCE, profile=PROFILE, trim_adapters=True, genomic_evidence=evidence)
    assert not view.barrier_intervals
    with pytest.raises(ValueError, match='outside'):
        infer_read_ends(SEQUENCE, profile=PROFILE, genomic_evidence=replace(evidence, exact_intervals=((0, 1000),)))
    with pytest.raises(ValueError, match='source'):
        GenomicQueryEvidence(exact_intervals=((0, 1),))


def test_reference_context_is_lazy_when_no_technical_motif_matches():
    def not_needed():
        raise AssertionError('Unmatched records must not expand CIGAR/MD bases')
    assert not infer_read_ends(LEFT + RIGHT, profile=PROFILE, genomic_evidence=not_needed).annotations


def test_local_allele_assembly_uses_one_side_and_keeps_original_quality_footprint():
    read = record('local', '1', 1000, f'12M{len(MOTIF)}I100N12M', SEQUENCE)
    before = read.to_string()
    left = collector().locus_read_from_pysam_aligned_segment(read, 1003, 1004)
    right = collector().locus_read_from_pysam_aligned_segment(read, 1115, 1116)
    assert left.sequence == LEFT and right.sequence == RIGHT
    assert left.read_base0_start_inclusive == right.read_base0_start_inclusive == 3
    assert len(left.quality_scores) == len(right.quality_scores) == 12
    assert left.source_read_views[0][1].original_sequence == SEQUENCE
    assert read.to_string() == before


@pytest.mark.parametrize('reverse', [False, True])
def test_sv_observations_preserve_parent_identity_without_cross_linker_gaps(reverse):
    sequence = reverse_complement_dna(SEQUENCE) if reverse else SEQUENCE
    read = record('parent', '1', 1000, f'12M{len(MOTIF)}I100N12M', sequence, flag=16 if reverse else 0)
    before = read.to_string()
    reasons = Counter()
    observations = _build_segment(('a', 'parent', 0), [read], collector(), reasons, {})
    assert len(observations) == 4  # Two intervals, each in two orientations; one parent.
    assert {observation.identity for observation in observations} == {('a', 'parent', 0)}
    assert all(not observation.breaks and len(observation.sequence) == 12 for observation in observations)
    assert len({observation.key for observation in observations}) == 4
    assert _record_gaps(read) and not _record_gaps(read, collector())
    assert reasons['internal_technical_join_not_spanning_evidence'] == 1
    assert read.to_string() == before


def test_linked_supplementary_records_cannot_bridge_the_same_internal_barrier():
    records = link([
        record('parent', '1', 1000, f'12M{len(SEQUENCE) - 12}S', SEQUENCE),
        record('parent', '2', 2000, f'{len(SEQUENCE) - 12}S12M', SEQUENCE, flag=2048),
    ])
    originals = [read.to_string() for read in records]
    observations = _build_segment(('a', 'parent', 0), records, collector(), Counter(), {})
    assert len(observations) == 4
    assert all(not observation.breaks and len(observation.sequence) == 12 for observation in observations)
    assert {observation.fragment for observation in observations} == {('a', 'parent')}
    assert all(len(observation.records) == 2 for observation in observations)
    assert [read.to_string() for read in records] == originals


def test_original_bam_and_raw_fastq_cannot_restore_a_spanning_witness():
    read = record('parent', '1', 1000, f'12M{len(MOTIF)}I100N12M', SEQUENCE)
    raw = [('a', 'parent', 0, SEQUENCE)]
    evidence = evidence_from_records([read], (), raw, read_collector=collector())
    assert {row[0] for row in evidence['reads']} == {('a', 'parent')}
    matches = exact_path_matches(SEQUENCE, evidence['reads'], min_match=10, seed_size=8)
    assert interval_spots(matches, 0, len(SEQUENCE)) == []
    assert not overlap_chain_search(matches, 0, len(SEQUENCE), min_overlap=1)['chain']
    assert {row[2] for row in evidence['reads']} == {LEFT, RIGHT}
    row = dict(query_sequence=SEQUENCE, query_aligned_interval=[0, 12],
               usable_query_intervals=[[0, 12], [12 + len(MOTIF), len(SEQUENCE)]])
    assert not aligned_exon_clip_contains(row, SEQUENCE[5:25])


def test_withheld_genomic_competitor_remains_an_explicit_alternative():
    from isovar.sv_rna_assessment import _alternative_placements

    read = record('competitor', '1', 1000, f'{len(SEQUENCE)}M', SEQUENCE)
    evidence = evidence_from_records([read], (), read_collector=collector())
    alternatives = list(_alternative_placements(SEQUENCE, [None] * len(SEQUENCE),
                                                [('a', 'competitor')], evidence['placements']))
    assert len(alternatives) == 1
    assert alternatives[0]['sequence_policy_eligible'] is False
    assert alternatives[0]['original_alignment']['query_sequence'] == SEQUENCE


def test_routed_assessment_matches_full_library_with_explicit_sequence_policy(tmp_path):
    scenario = Scenario()
    reads = scenario.reads()
    # A real reference-like exact control and a simulated motif-containing raw end.
    reads.append(record('technical', '1', 1000, f'12M{len(SEQUENCE) - 12}S', SEQUENCE))
    bam_path = write_bam(tmp_path / 'original.bam', reads)
    result = scenario.run(bam_path)
    references = (scenario.donor_ref, scenario.acceptor_ref)
    policy = collector()
    raw = [('a', 'technical', 0, SEQUENCE)]
    with pysam.AlignmentFile(bam_path) as bam:
        expected = assess_sv_rna_paths(result, bam, references, raw_reads=raw, read_collector=policy)
    query = PathEvidenceQuery('view', tuple(path['sequence'] for path in result['paths']), references)
    with prepare_rna_library(bam_path, tmp_path / 'cache').open() as library:
        with route_path_evidence(library, (query,), tmp_path / 'routes.sqlite', raw_reads=raw,
                                 read_collector=policy) as routed:
            evidence = routed.evidence('view')
    actual = assess_sv_rna_paths(result, None, references, evidence=evidence, read_collector=policy)
    assert actual == expected
    assert assess_sv_rna_paths(result, None, references, evidence=evidence,
                               read_collector=policy, **actual['parameters']) == actual
    with pytest.raises(ValueError, match='sequence policy'):
        assess_sv_rna_paths(result, None, references, evidence=evidence)


def test_dense_and_bounded_reconstruction_apply_the_same_internal_policy(tmp_path):
    scenario = Scenario()
    profile = ReadEndProfile('simulated-join-motif', adapters=(Adapter(
        'synthetic', scenario.sequence[235:265], placement='internal', role='linker'),))
    policy = ReadCollector(read_end_profile=profile, trim_adapters=True, use_soft_clipped_bases=True)
    bam_path = write_bam(tmp_path / 'original.bam', scenario.reads())
    bounded = scenario.run(bam_path, read_collector=policy)
    dense = scenario.run(bam_path, read_collector=policy, dense_support=True)
    assert {path['path_id'] for path in bounded['paths']} == {path['path_id'] for path in dense['paths']}
    for before, after in zip(bounded['paths'], dense['paths']):
        assert [join['direct_support']['fragments'] for join in before['junctions']] == [
            join['direct_support']['fragments'] for join in after['junctions']]
        assert [value['amino_acids'] for value in before['translations']] == [
            value['amino_acids'] for value in after['translations']]


def test_missing_md_is_unknown_not_a_negative_reference_measurement():
    read = record('unknown', '1', 1000, f'{len(SEQUENCE)}M', SEQUENCE)
    evidence = genomic_query_evidence_from_alignment(read)
    assert not evidence.exact_intervals
    assert evidence.limitations == ('MD_reference_claim_unavailable',)


@pytest.mark.parametrize('md', ['99999999999999', '57broken', 57, '0A0'])
def test_malformed_md_does_not_discard_the_original_read(md):
    read = record('malformed', '1', 1000, f'{len(SEQUENCE)}M', SEQUENCE)
    read.set_tag('MD', md)
    before = read.to_string()
    evidence = genomic_query_evidence_from_alignment(read)
    assert not evidence.exact_intervals
    assert evidence.limitations == ('MD_reference_claim_unavailable_or_malformed',)
    assert read.to_string() == before


def test_conflicting_cigar_and_md_claims_do_not_rescue_a_technical_candidate():
    read = record('conflicting', '1', 1000, f'{len(SEQUENCE)}=', SEQUENCE)
    read.set_tag('MD', '20A' + str(len(SEQUENCE) - 21))
    view = read_sequence_view_from_alignment(read, PROFILE, trim_adapters=True)
    assert view.barrier_intervals
    assert view.genomic_evidence.limitations == ('CIGAR_MD_equality_claims_conflict',)


@pytest.mark.parametrize('pt,status', [(20, 'estimate'), (-1, 'estimation_failed'),
                                      (0, 'estimation_failed'), ('bad', 'invalid_tag')])
def test_upstream_tail_estimates_and_failures_stay_distinct(pt, status):
    read = record('tail', '1', 1000, f'{len(SEQUENCE)}M', SEQUENCE)
    read.set_tag('pt', pt)
    before = read.to_string()
    view = read_sequence_view_from_alignment(read)
    assert view.upstream_poly_a_status == status
    assert view.upstream_poly_a_length == (None if isinstance(pt, str) else pt)
    assert read.to_string() == before


@pytest.mark.parametrize('reverse', [False, True])
def test_explicit_adapter_order_preserves_conflicts_without_reordering_reads(reverse):
    other = Adapter('second', 'CTACAGGTCCGATGCTATCCAGGTCAT', role='linker', placement='internal')
    profile = replace(PROFILE, adapters=(ADAPTER, other), expected_adapter_order=(ADAPTER.name, other.name))
    sequence = LEFT + MOTIF + RIGHT + other.sequence + LEFT
    if reverse:
        sequence = reverse_complement_dna(sequence)
    view = infer_read_ends(sequence, profile=profile, reverse=reverse, trim_adapters=True)
    assert view.adapter_order_status == 'consistent'
    wrong = replace(profile, expected_adapter_order=profile.expected_adapter_order[::-1])
    assert infer_read_ends(sequence, profile=wrong, reverse=reverse).adapter_order_status == 'inconsistent'
    assert view.original_sequence == sequence
    assert infer_read_ends(SEQUENCE, profile=profile).adapter_order_status == 'insufficient_evidence'


@pytest.mark.parametrize('same_parent', [True, False])
def test_two_direction_extension_cannot_switch_parent_segment(same_parent):
    from isovar.sv_rna import RnaObservation, _ObservationIndex, _extend_both

    observations = []
    for key, sequence, start, interval in [('first', 'CGTT', 100, (0, 4)), ('second', 'AACG', 98, (10, 14))]:
        positions = tuple(('1', offset, '+') for offset in range(start, start + 4))
        identity = ('a', 'parent' if same_parent else key, 0)
        for reverse in (False, True):
            observations.append(RnaObservation(key + ('-' if reverse else '+'), identity, (key,),
                reverse_complement_dna(sequence) if reverse else sequence,
                tuple((c, p, '-') for c, p, _ in reversed(positions)) if reverse else positions,
                (), reverse, interval, False, False))
    paths = list(_extend_both('CG', (('1', 100, '+'), ('1', 101, '+')),
                            _ObservationIndex(observations), 2, (1, 0, 0), 8, [], set()))
    assert sorted(path[0] for path in paths) == (['AACG', 'CGTT'] if same_parent else ['AACGTT'])
