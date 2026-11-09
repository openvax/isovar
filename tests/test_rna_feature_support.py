"""Observation accounting and geometry, independent of antigen category."""

from dataclasses import replace
import json
import random

import pytest

from isovar.rna_feature_support import (
    RNAFeature, RNASamplingGeometry, RNASamplingProfile, effective_rna_feature_lengths,
    measure_rna_feature, validate_rna_feature_support,
)
from isovar.rna_path_linkage import OriginalPathMatch


def match(name, start, end, mate=0, placement="one"):
    return OriginalPathMatch(("RG", name, mate), placement, start, end)


def profile(read_length=100, *, depth=1000000):
    return RNASamplingProfile("whole-input-test", (RNASamplingGeometry(((0, read_length),), read_length),),
                              depth, depth, "whole_library")


def feature(start=500, end=501, kind="allele", **kwargs):
    return RNAFeature("feature", "path", ((start, end),), kind, **kwargs)


@pytest.mark.parametrize("length", [50, 100, 150])
@pytest.mark.parametrize("size,anchor", [(1, 0), (1, 10), (0, 10), (20, 10), (75, 0), (300, 10)])
def test_se_analytic_opportunities(length, size, anchor):
    f = feature(500 - anchor, 500 + size + anchor)
    expected = max(0, length - size - 2 * anchor + 1)
    assert effective_rna_feature_lengths(f, profile(length)) == dict(reads=expected, fragments=expected)


def brute_opportunity(f, geometry):
    """Independent base-set oracle for the optimized interval calculation."""
    requirements = [set(range(a, b)) for a, b in f.intervals]
    reads = fragments = 0
    upper = f.path_length - geometry.template_length + 1 if f.path_length is not None else 1200
    lower = 0 if f.path_length is not None else -500
    for start in range(lower, upper):
        ends = [set(range(start + a, start + b)) for a, b in geometry.segments]
        def accepts(bases):
            covered = [required <= bases for required in requirements]
            return all(covered) if f.mode == "all" else any(covered)
        reads += sum(accepts(end) for end in ends)
        fragments += (accepts(set.union(*ends)) if f.fragment_mode == "union"
                      else any(accepts(end) for end in ends))
    return reads / len(geometry.segments), fragments


@pytest.mark.parametrize("fragment_mode", ["union", "any_read"])
def test_geometry_matches_independent_oracle(fragment_mode):
    rng = random.Random(11)
    for _ in range(60):
        insert = rng.randint(60, 250)
        length = rng.randint(20, min(insert, 100))
        geometry = RNASamplingGeometry(((0, length), (insert - length, insert)), insert)
        start = rng.randint(0, 150)
        f = feature(start, start + rng.randint(1, 80), path_length=300, fragment_mode=fragment_mode)
        if rng.choice([True, False]):
            f = replace(f, intervals=(*f.intervals, (240, 250)), mode=rng.choice(["all", "any"]))
        expected = brute_opportunity(f, geometry)
        actual = effective_rna_feature_lengths(f, RNASamplingProfile("test", (geometry,)))
        assert (actual["reads"], actual["fragments"]) == expected


def test_single_end_counts_and_all_normalizations_agree():
    record = measure_rna_feature(feature(), [match("a", 450, 550), match("b", 490, 590)],
                                 evidence_scope=("sample", "input"), profile=profile())
    values = record["values"]
    for suffix in ("", "_per_million", "_per_effective_kb", "_per_effective_kb_per_million"):
        assert values["reads" + suffix] == values["fragments" + suffix]
    assert values["reads"] == 2
    assert values["reads_per_million"] == 2
    assert values["reads_per_effective_kb_per_million"] == 20
    assert validate_rna_feature_support(json.loads(json.dumps(record))) == json.loads(json.dumps(record))


def test_pair_union_and_read_opportunities_are_not_doubled():
    p = RNASamplingProfile("PE", (RNASamplingGeometry(((0, 100), (50, 150)), 150),),
                           2000000, 1000000, "whole_library")
    f = feature()
    assert effective_rna_feature_lengths(f, p) == dict(reads=100, fragments=150)
    record = measure_rna_feature(f, [match("a", 450, 550, 64), match("a", 480, 580, 128)],
                                 evidence_scope=("S", "L"), profile=p)
    assert record["values"]["reads"] == 2
    assert record["values"]["fragments"] == 1
    assert record["values"]["reads_per_million"] == record["values"]["fragments_per_million"] == 1


def test_joint_template_support_does_not_create_supporting_reads():
    f = feature(500, 650, kind="coding_interval")
    p = RNASamplingProfile("PE", (RNASamplingGeometry(((0, 100), (75, 175)), 175),))
    record = measure_rna_feature(f, [match("a", 490, 590, 64), match("a", 565, 665, 128)],
                                 evidence_scope=("S", "L"), profile=p)
    assert record["values"]["reads"] == 0
    assert record["values"]["fragments"] == 1
    assert record["values"]["effective_length_reads"] == 0
    assert record["values"]["effective_length_fragments"] == 26
    assert record["values"]["reads_per_effective_kb"] is None
    assert len(record["joint_only_fragment_ids"]) == 1


def test_individually_spanning_fragment_rule_matches_denominator():
    f = feature(500, 650, kind="coding_interval", fragment_mode="any_read")
    p = RNASamplingProfile("PE", (RNASamplingGeometry(((0, 100), (75, 175)), 175),))
    rows = [match("a", 490, 590, 64), match("a", 565, 665, 128)]
    record = measure_rna_feature(f, rows, evidence_scope=("S", "L"), profile=p)
    assert record["values"]["fragments"] == record["values"]["reads"] == 0
    assert record["values"]["effective_length_fragments"] == 0
    assert not record["joint_only_fragment_ids"]
    shorter = replace(f, intervals=((500, 550),))
    supported = measure_rna_feature(shorter, rows, evidence_scope=("S", "L"), profile=p)
    assert supported["values"]["fragments"] == supported["values"]["reads"] == 1
    assert supported["values"]["effective_length_fragments"] == 102
    # For SE both rules coincide, including normalized values.
    for mode in ("union", "any_read"):
        se = measure_rna_feature(replace(shorter, fragment_mode=mode), [match("s",490,590)],
            evidence_scope=("S", "L"), profile=profile())
        assert se["values"]["reads_per_effective_kb_per_million"] == se["values"]["fragments_per_effective_kb_per_million"]


def test_disconnected_bases_are_not_tiled_by_unsequenced_insert():
    rows = [match("a", 450, 550, 64), match("a", 600, 700, 128)]
    record = measure_rna_feature(feature(500, 650), rows, evidence_scope=("S", "L"))
    assert record["values"]["fragments"] == 0
    # A phase-only compound feature need not claim the intervening bases.
    f = replace(feature(), intervals=((500, 501), (640, 641)), kind="compound")
    linked = measure_rna_feature(f, rows, evidence_scope=("S", "L"))
    assert linked["values"]["fragments"] == 1
    assert linked["values"]["reads"] == 0


def test_alternative_placements_and_separate_templates_cannot_create_phase():
    f = replace(feature(), intervals=((500, 501), (640, 641)), kind="compound")
    for rows in ([match("a", 450, 550, 64, "left"), match("a", 600, 700, 64, "right")],
                 [match("a", 450, 550, 64), match("b", 600, 700, 128)]):
        assert measure_rna_feature(f, rows, evidence_scope=("S", "L"))["values"]["fragments"] == 0


def test_partial_union_counts_a_template_once():
    f = replace(feature(), intervals=((500, 501), (640, 641)), kind="compound", mode="any")
    rows = [match("a", 450, 550, 64), match("a", 600, 700, 128)]
    record = measure_rna_feature(f, rows + rows, evidence_scope=("S", "L"))
    assert record["values"]["reads"] == 2
    assert record["values"]["fragments"] == 1
    assert record["feature"]["mode"] == "any"


@pytest.mark.parametrize("kind", ["allele", "junction", "coding_interval", "compound"])
def test_feature_label_does_not_change_measurement(kind):
    record = measure_rna_feature(feature(kind=kind), [match("a", 450, 550)],
                                 evidence_scope=("S", "L"), profile=profile())
    assert record["values"]["reads_per_effective_kb_per_million"] == 10


def test_depth_and_geometry_normalization_recover_same_sampling_density():
    # Uniform one read at every start: a point gets 100 witnesses, a 75-nt
    # interval gets 26. Their length-normalized support must agree.
    rows = [match(str(i), i, i + 100) for i in range(300, 800)]
    values = [measure_rna_feature(feature(500, end), rows, evidence_scope=("S", "L"),
                                  profile=profile())["values"] for end in (501, 575)]
    assert [v["reads"] for v in values] == [100, 26]
    assert values[0]["reads_per_effective_kb_per_million"] == values[1]["reads_per_effective_kb_per_million"]
    doubled = rows + [replace(r, identity=("RG2", *r.identity[1:])) for r in rows]
    v = measure_rna_feature(feature(), doubled, evidence_scope=("S", "L"),
                            profile=profile(depth=2000000))["values"]
    assert v["reads"] == 200
    assert v["reads_per_effective_kb_per_million"] == values[0]["reads_per_effective_kb_per_million"]


def test_missing_and_unobservable_are_distinct():
    r = measure_rna_feature(feature(), [], evidence_scope=("S", "L"))
    assert r["values"]["reads"] == 0
    assert r["values"]["reads_per_million"] is None
    assert r["unavailable"]["reads_per_effective_kb"] == "sampling_geometry_unavailable"
    r = measure_rna_feature(feature(500, 800), [], evidence_scope=("S", "L"), profile=profile())
    assert r["values"]["reads_per_million"] == 0
    assert r["unavailable"]["reads_per_effective_kb"] == "feature_unobservable"
    with pytest.raises(ValueError, match="zero effective"):
        measure_rna_feature(feature(500, 800), [match("a", 450, 850)], evidence_scope=("S", "L"), profile=profile())


def test_original_allele_call_restriction_and_source_identity():
    f = replace(feature(), observation_keys=("passed",))
    rows = [match("a", 450, 550, placement="passed"), match("b", 450, 550, placement="failed")]
    a = measure_rna_feature(f, rows, evidence_scope=("S", "L"))
    b = measure_rna_feature(f, rows, evidence_scope=("S", "another-library"))
    assert a["values"]["reads"] == 1
    assert a["read_ids"] != b["read_ids"]


def test_denominator_and_json_integrity():
    with pytest.raises(ValueError, match="whole_library"):
        RNASamplingProfile("regional", library_reads=100)
    with pytest.raises(ValueError, match="Single-end"):
        replace(profile(), library_reads=2000000)
    record = measure_rna_feature(feature(), [], evidence_scope=("S", "L"), profile=profile())
    record["values"]["reads"] = 100
    with pytest.raises(ValueError, match="identity mismatch"):
        validate_rna_feature_support(record)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_original_bam_coding_junction_and_allele_features(tmp_path, strand):
    import pysam
    from isovar import measure_rna_path_features, prepare_rna_library, ReadCollector
    from tests.test_antigen_isoforms import target_models
    from tests.test_sv_rna import aligned, write_bam

    targets, paths = target_models(strand)
    sequence, positions = targets[0].reference.sequence, paths[0]
    reads = aligned("shared", sequence[150:250], positions[150:250])
    reads += aligned("pair", sequence[120:220], positions[120:220], flag=65)
    reads += aligned("pair", sequence[190:290], positions[190:290], flag=129)
    # A mismatching original must not support the nominated nucleotide allele.
    wrong = sequence[150:200] + ("A" if sequence[200] != "A" else "T") + sequence[201:250]
    reads += aligned("other_allele", wrong, positions[150:250])
    bam_path = write_bam(tmp_path / "features.bam", reads + reads[:1])
    features = [RNAFeature("allele", "T1", ((200, 201),), "allele"),
                RNAFeature("junction", "T1", ((170, 190),), "junction"),
                RNAFeature("coding", "T1", ((160, 260),), "coding_interval")]
    policy = ReadCollector.for_rna_paths()
    prepared = prepare_rna_library(bam_path, tmp_path / "index", read_collector=policy)
    with pysam.AlignmentFile(bam_path) as bam, prepared.open() as index:
        records = measure_rna_path_features(bam, sequence=sequence, positions=positions,
            features=features, sample_id="S", source="L", library_index=index)
    assert [(r["values"]["reads"], r["values"]["fragments"]) for r in records] == [(3, 2), (3, 3), (0, 1)]
    assert all(validate_rna_feature_support(r) == r for r in records)
    assert records[2]["joint_only_fragment_ids"]
    assert all(r["values"]["library_fragments"] is None for r in records)
    assert records[0]["provenance"]["read_policy"] == prepared.identity["eligibility"]
    assert records[0]["provenance"]["attribution"] == "compatible_coordinate_path"


def test_rna_features_cli_keeps_library_identity_and_unknown_denominator(tmp_path):
    from isovar.cli.commands import run
    from tests.test_sv_rna import aligned, write_bam
    sequence = "ACGT" * 25
    positions = [("1", p, "+") for p in range(100, 200)]
    bam = write_bam(tmp_path / "input.bam", aligned("one", sequence, positions))
    request = tmp_path / "request.json"
    from dataclasses import asdict
    request.write_text(json.dumps(dict(paths=[dict(sequence=sequence, positions=positions,
        features=[asdict(RNAFeature("local", "path", ((40, 60),), "coding_interval"))])])))
    output = tmp_path / "support.json"
    args = ["rna-features", "--bam", str(bam), "--features", str(request), "--sample-id", "S",
            "--source-id", "L", "--output", str(output)]
    run(args)
    record, = json.loads(output.read_text())["measurements"]
    assert record["evidence_scope"] == ["S", "L"]
    assert record["values"]["fragments"] == record["values"]["reads"] == 1
    assert record["values"]["fragments_per_million"] is None
    validate_rna_feature_support(record)
    with pytest.raises(SystemExit):
        run(args)


@pytest.mark.parametrize("field,value", [("normalization_model", "unknown"),
    ("geometry_scope", "complete_path"), ("evidence_scope", ["S"]),
    ("unavailable", {"reads": "missing"})])
def test_saved_measurements_reject_inconsistent_rehashed_metadata(field, value):
    from isovar.rna_evidence import content_identifier
    record = measure_rna_feature(feature(), [match("a", 450, 550)], evidence_scope=("S", "L"))
    record[field] = value
    del record["measurement_id"]
    record["measurement_id"] = content_identifier("rna_feature", record)
    with pytest.raises(ValueError):
        validate_rna_feature_support(record)


def test_mixed_length_weighted_interval_geometry_and_unbounded_edges():
    geometries = (RNASamplingGeometry(((0, 60),), 60, 2),
                  RNASamplingGeometry(((40, 120), (0, 80)), 120, 3),
                  RNASamplingGeometry(((0, 30), (150, 200)), 200, 0.5))
    for path_length in (None, 300):
        for start, end in ((0, 1), (5, 75), (100, 175), (290, 300)):
            f = feature(start, end, path_length=path_length)
            expected = [brute_opportunity(f, g) for g in geometries]
            got = effective_rna_feature_lengths(f, RNASamplingProfile('mixture', geometries))
            assert got['reads'] == pytest.approx(sum(g.weight * len(g.segments) * v[0]
                for g, v in zip(geometries, expected)) / sum(g.weight * len(g.segments) for g in geometries))
            assert got['fragments'] == pytest.approx(sum(g.weight * v[1]
                for g, v in zip(geometries, expected)) / sum(g.weight for g in geometries))


def test_later_library_normalization_retains_raw_identity_and_feature():
    from isovar import normalize_rna_feature_support
    original = measure_rna_feature(feature(), [match("a", 450, 550)], evidence_scope=("S", "L"))
    normalized = normalize_rna_feature_support(original, profile())
    assert original['values']['reads_per_million'] is None
    assert normalized['values']['reads_per_million'] == 1
    assert normalized['feature'] == original['feature']
    assert normalized['read_ids'] == original['read_ids']
    assert normalized['fragment_ids'] == original['fragment_ids']
    assert normalized['measurement_id'] != original['measurement_id']
    assert normalized['provenance']['normalizes_measurement_id'] == original['measurement_id']


def test_local_feature_ignores_external_mismatch_but_not_internal_mismatch(tmp_path):
    import pysam
    from isovar import measure_rna_path_features
    from tests.test_sv_rna import aligned, write_bam
    sequence = "ACGT" * 25
    positions = [("1", i, "+") for i in range(100, 200)]
    read_sequence = sequence[:50] + "A" + sequence[51:]
    assert sequence[50] != "A"
    path = write_bam(tmp_path / "local.bam", aligned("read", read_sequence, positions))
    features = [RNAFeature("left", "path", ((20, 30),), "coding_interval"),
                RNAFeature("whole", "path", ((20, 80),), "coding_interval"),
                RNAFeature("compound", "path", ((20, 30), (70, 80)), "compound")]
    with pysam.AlignmentFile(path) as bam:
        records = measure_rna_path_features(bam, sequence=sequence, positions=positions,
            features=features, sample_id="S", source="L")
    assert [r['values']['reads'] for r in records] == [1, 0, 1]
    assert [r['values']['fragments'] for r in records] == [1, 0, 1]


def test_repeated_path_offsets_do_not_phase_distinct_copies(tmp_path):
    import pysam
    from isovar import measure_rna_path_features
    from tests.test_sv_rna import aligned, write_bam
    sequence = "ACGT" * 25
    positions = [("1", i, "+") for i in range(100, 200)]
    path = write_bam(tmp_path / "copies.bam", aligned("read", sequence, positions))
    features = [RNAFeature("copies", "duplication", ((40, 45), (140, 145)), "compound"),
                RNAFeature("partial", "duplication", ((40, 45), (140, 145)), "compound", mode="any")]
    with pysam.AlignmentFile(path) as bam:
        records = measure_rna_path_features(bam, sequence=sequence * 2, positions=positions * 2,
            features=features, sample_id="S", source="L")
    assert [r['values']['reads'] for r in records] == [0, 1]
    assert [r['values']['fragments'] for r in records] == [0, 1]


@pytest.mark.parametrize("event", ["insertion", "deletion", "splice", "fusion"])
def test_local_bam_features_use_alternate_rna_coordinates(tmp_path, event):
    import pysam
    from isovar import measure_rna_path_features
    from tests.test_sv_rna import aligned, record, write_bam
    rng = random.Random(22)
    sequence = "".join(rng.choice("ACGT") for _ in range(100))
    positions = [("1", i, "+") for i in range(100, 200)]
    if event == "insertion":
        sequence = sequence[:40] + "AACGT" + sequence[40:]
        positions = positions[:40] + [None] * 5 + positions[40:]
        reads = [record("read", "1", 100, "40M5I60M", sequence)]
        interval = (35, 50)
    elif event in ("deletion", "splice"):
        gap = 10 if event == "deletion" else 100
        positions = positions[:40] + [("1", i + gap, "+") for i in range(140, 200)]
        reads = [record("read", "1", 100, "40M%d%s60M" % (gap, "D" if event == "deletion" else "N"), sequence)]
        interval = (35, 45)
    else:
        positions = positions[:40] + [("2", i, "+") for i in range(500, 560)]
        reads = aligned("read", sequence, positions)
        interval = (35, 45)
    path = write_bam(tmp_path / "event.bam", reads)
    f = RNAFeature(event, "alternate-path", (interval,), "allele" if event == "insertion" else "junction")
    with pysam.AlignmentFile(path) as bam:
        measured, = measure_rna_path_features(bam, sequence=sequence, positions=positions,
            features=[f], sample_id="S", source="L", profile=profile(len(sequence)))
    assert measured['values']['reads'] == measured['values']['fragments'] == 1
    assert measured['values']['effective_length_reads'] == len(sequence) - (interval[1] - interval[0]) + 1
