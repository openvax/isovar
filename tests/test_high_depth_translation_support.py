"""Equivalence and high-depth runtime regression tests for #544."""

import pickle
import random
import time

import pytest

from isovar.allele_read import AlleleRead
from isovar.assembly import merge_identical_sequences
from isovar.protein_sequence_creator import ProteinSequenceCreator
from isovar.variant_sequence import VariantSequence
from isovar.variant_sequence_creator import VariantSequenceCreator


def _per_group_reference_sequences(creator, reads):
    """Independent per-group, per-length baseline matching pre-#544 semantics."""
    groups = creator._transcript_read_groups(reads)
    if groups is None:
        group_reads_list = [list(reads)]
    else:
        group_reads_list = [compatible_reads for _, compatible_reads in groups]

    sequences = []
    for compatible_reads in group_reads_list:
        if creator.protein_sequence_preference == "support":
            sequences.extend(
                creator._variant_sequence_creator.reads_to_variant_sequences(
                    None, compatible_reads
                )
            )
        else:
            per_group = []
            for length in creator.candidate_context_lengths():
                sub = (
                    creator._variant_sequence_creator
                    if length == creator.protein_sequence_length
                    else creator._make_variant_sequence_creator(3 * length + 2)
                )
                per_group.extend(sub.reads_to_variant_sequences(None, compatible_reads))
            sequences.extend(merge_identical_sequences(per_group))
    return merge_identical_sequences(sequences)


def _snapshot(sequences):
    """Capture every observable property of candidate VariantSequences."""
    return [
        (
            seq.prefix,
            seq.alt,
            seq.suffix,
            seq.compatible_transcript_ids,
            frozenset(
                (r, r.compatible_transcript_ids)
                for r in seq.reads
            ),
            seq.coverage().tolist(),
            dict(seq._source_alignments),
        )
        for seq in sequences
    ]


@pytest.mark.parametrize("alt", ["C", "", "CTG"])
@pytest.mark.parametrize("preference", ["support", "balanced"])
@pytest.mark.parametrize("assembly", [False, True])
def test_shared_index_matches_per_group_baseline_across_modes(alt, preference, assembly):
    rng = random.Random(5440 + len(alt) + (10 if preference == "balanced" else 0) + int(assembly))
    ref_prefix = "".join(rng.choices("ACGT", k=90))
    ref_suffix = "".join(rng.choices("ACGT", k=90))
    tx_pool = [
        frozenset({"T1", "T2", "T3", "T4"}),
        frozenset({"T1", "T2"}),
        frozenset({"T2", "T3"}),
        frozenset({"T3", "T4"}),
        frozenset({"T1"}),
        None,
    ]
    reads = []
    for i in range(60):
        p_len = rng.randint(0 if alt else 2, 40)
        s_len = rng.randint(0 if alt else 2, 40)
        pref = ref_prefix[-p_len:] if p_len else ""
        suff = ref_suffix[:s_len]
        if i % 7 == 0 and pref:
            pos = rng.randrange(len(pref))
            pref = pref[:pos] + ("T" if pref[pos] != "T" else "A") + pref[pos + 1 :]
        if not (pref or alt or suff):
            continue
        seg = ("", f"read-{i}", 64)
        reads.append(
            AlleleRead(
                prefix=pref,
                allele=alt,
                suffix=suff,
                name=f"read-{i}",
                source_alignments=((seg, (0, 1000 - p_len, f"{len(pref) + len(alt) + len(suff)}M", False)),),
                compatible_transcript_ids=rng.choice(tx_pool),
            )
        )

    # Also add a merged mate pair + retained single mate (multi-occurrence segment)
    # and a competing alternative alignment placement (conflicting segment).
    shared_seg_r1 = ("", "merged-frag", 64)
    shared_seg_r2 = ("", "merged-frag", 128)
    reads.append(
        AlleleRead(
            prefix=ref_prefix[-25:],
            allele=alt,
            suffix=ref_suffix[:25],
            name="merged-frag",
            source_read_count=2,
            source_alignments=(
                (shared_seg_r1, (0, 975, "51M", False)),
                (shared_seg_r2, (0, 975, "51M", True)),
            ),
            compatible_transcript_ids=frozenset({"T1", "T2", "T3", "T4"}),
        )
    )
    reads.append(
        AlleleRead(
            prefix=ref_prefix[-15:],
            allele=alt,
            suffix=ref_suffix[:15],
            name="merged-frag",
            source_read_count=1,
            source_alignments=((shared_seg_r1, (0, 975, "51M", False)),),
            compatible_transcript_ids=frozenset({"T1", "T2"}),
        )
    )

    creator = ProteinSequenceCreator(
        protein_sequence_length=25,
        protein_context_peptide_length=15,
        protein_sequence_preference=preference,
        variant_sequence_assembly=assembly,
        min_variant_sequence_coverage=2,
    )
    optimized = creator.variant_sequences_from_reads(None, reads)
    baseline = _per_group_reference_sequences(creator, reads)
    assert _snapshot(optimized) == _snapshot(baseline)

    # Now add a competing alternative placement for one segment and verify exact
    # equivalence on the conflict-aware path too.
    conflict_seg = ("", "multi-mapped", 64)
    reads_with_conflict = list(reads) + [
        AlleleRead(
            prefix=ref_prefix[-20:],
            allele=alt,
            suffix=ref_suffix[:20],
            name="multi-mapped",
            source_alignments=((conflict_seg, (0, 980, "41M", False)),),
            compatible_transcript_ids=frozenset({"T1", "T2"}),
        ),
        AlleleRead(
            prefix="G" + ref_prefix[-19:],
            allele=alt,
            suffix=ref_suffix[:20],
            name="multi-mapped",
            source_alignments=((conflict_seg, (0, 990, "41M", False)),),
            compatible_transcript_ids=frozenset({"T1", "T2"}),
        ),
    ]
    optimized_conflict = creator.variant_sequences_from_reads(None, reads_with_conflict)
    baseline_conflict = _per_group_reference_sequences(creator, reads_with_conflict)
    assert _snapshot(optimized_conflict) == _snapshot(baseline_conflict)


def test_allele_read_and_variant_sequence_hash_caching_and_pickle_compat():
    read_a = AlleleRead(
        prefix="ACGT",
        allele="C",
        suffix="TGCA",
        name="r1",
        source_alignments=((("", "r1", 64), (0, 100, "9M", False)),),
        compatible_transcript_ids=frozenset({"T1", "T2"}),
    )
    read_b = read_a.with_compatible_transcript_ids(frozenset({"T1"}))
    assert read_a == read_b
    assert hash(read_a) == hash(read_b)
    assert hash(read_a) == hash(read_a._values)

    restored = pickle.loads(pickle.dumps(read_a))
    assert restored == read_a
    assert hash(restored) == hash(read_a)

    vs = VariantSequence("ACGT", "C", "TGCA", [read_a])
    assert hash(vs) == hash(vs._values)
    restored_vs = pickle.loads(pickle.dumps(vs))
    assert restored_vs == vs
    assert hash(restored_vs) == hash(vs)


def test_high_depth_multi_transcript_locus_runs_quickly_and_preserves_all_paths():
    """Synthetic CDKN2A-like high-depth frameshift locus across 4 transcript groups."""
    rng = random.Random(544)
    ref_prefix = "".join(rng.choices("ACGT", k=120))
    ref_suffix = "".join(rng.choices("ACGT", k=120))
    alt_branch_suffix = "T" + ref_suffix[1:120] if ref_suffix[0] != "T" else "A" + ref_suffix[1:120]
    tx_sets = [
        frozenset({"T1", "T2", "T3", "T4"}),
        frozenset({"T1", "T2"}),
        frozenset({"T2", "T3"}),
        frozenset({"T3", "T4"}),
        frozenset({"T1"}),
        frozenset({"T4"}),
    ]
    reads = []
    for i in range(3000):
        p_len = rng.randint(12, 74)
        s_len = rng.randint(12, 74)
        pref = ref_prefix[-p_len:]
        suff = alt_branch_suffix[:s_len] if i < 200 else ref_suffix[:s_len]
        if i >= 200 and i % 15 == 0:
            idx = rng.randrange(len(pref))
            pref = pref[:idx] + ("A" if pref[idx] != "A" else "C") + pref[idx + 1 :]
        txs = tx_sets[0] if i < 2600 else rng.choice(tx_sets[1:])
        reads.append(
            AlleleRead(
                prefix=pref,
                allele="C",
                suffix=suff,
                name=f"r{i}",
                source_alignments=(
                    (("", f"r{i}", 64), (0, 21971154 - p_len, f"{p_len + 1 + s_len}M", False)),
                ),
                compatible_transcript_ids=txs,
            )
        )

    creator = ProteinSequenceCreator(
        max_protein_sequences_per_variant=1,
        variant_sequence_assembly=False,
        protein_sequence_preference="support",
    )
    t0 = time.perf_counter()
    sequences = creator.variant_sequences_from_reads(None, reads)
    elapsed = time.perf_counter() - t0

    assert elapsed < 15.0, f"High-depth locus took too long: {elapsed:.2f}s"
    assert sequences
    # All 4 transcript groups must remain represented.
    observed_tx_groups = {seq.compatible_transcript_ids for seq in sequences}
    assert observed_tx_groups == {
        frozenset({"T1"}),
        frozenset({"T2"}),
        frozenset({"T3"}),
        frozenset({"T4"}),
    }
    # Both main and competing branch suffixes must survive with expected support.
    assert any(seq.suffix.startswith(ref_suffix[:12]) and len(seq.reads) > 2000 for seq in sequences)
    assert any(seq.suffix.startswith(alt_branch_suffix[:12]) and len(seq.reads) >= 150 for seq in sequences)
    # Check that single-group VariantSequenceCreator matches for one group slice.
    t1_reads = [
        r.with_compatible_transcript_ids(frozenset({"T1"}))
        for r in reads
        if "T1" in r.compatible_transcript_ids
    ]
    t1_expected = merge_identical_sequences(
        VariantSequenceCreator(
            min_variant_sequence_coverage=creator.min_variant_sequence_coverage,
            preferred_sequence_length=creator._cdna_sequence_length,
            variant_sequence_assembly=False,
            min_flanking_sequence_length=creator.min_transcript_prefix_length,
        ).reads_to_variant_sequences(None, t1_reads)
    )
    t1_actual = [seq for seq in sequences if seq.compatible_transcript_ids == frozenset({"T1"})]
    assert _snapshot(t1_actual) == _snapshot(t1_expected)
