"""Control-derived, directional phasing calibration."""

import json

import pytest
from varcode import Variant

from isovar.allele_read import AlleleRead
from isovar.cli.isovar_calibrate_phasing import run as run_calibration
from isovar.phasing import _four_gamete
from isovar.phasing_calibration import (
    PhasingControl, calibrate_phasing_error_rates, one_sided_error_bound,
)
from isovar.read_evidence import ReadEvidence
from isovar.variant_helpers import trim_variant


META = dict(platform="ONT", basecaller="dorado-1", library="PCR-free-RNA",
            biological_context="osteosarcoma", reference="GRCh38", context="AAAAAA")


def _control(variant, truth, group, source, correct, wrong, other=0):
    start, ref, alt = trim_variant(variant)

    def reads(allele, names):
        return [AlleleRead("A", allele, "T", name,
                           source_alignments=(((group, name, 0), (0, start, "10M", False)),))
                for name in names]

    ref_names = ["r%d" % i for i in range(correct if truth == "ref" else wrong)]
    alt_names = ["a%d" % i for i in range(wrong if truth == "ref" else correct)]
    evidence = ReadEvidence(start, ref, alt, reads(ref, ref_names), reads(alt, alt_names),
                            reads("X", ["o%d" % i for i in range(other)]))
    return PhasingControl(variant, truth, group, source, **META, evidence=evidence)


def _target(variant, groups=("rg1",)):
    return {variant: {group: META for group in groups}}


def test_exact_bound_and_directional_control_rates_change_wrong_definite_call():
    variant = Variant("1", 100, "A", "G", normalize_contig_names=False)
    controls = [
        _control(variant, "ref", "rg1", "ref-source", 90, 10, other=12),
        _control(variant, "alt", "rg1", "alt-source", 99, 1),
    ]
    result = calibrate_phasing_error_rates(controls, _target(variant))
    ref_to_alt, alt_to_ref = result.rates[variant]
    assert ref_to_alt > 0.1 > alt_to_ref
    assert result.report[0]["wrong"] == 10
    assert result.report[0]["other"] == 12
    assert result.profile_rows()[0]["ref_to_alt"] == ref_to_alt
    assert _four_gamete(3, 47, 2, 48, phasing_error_rate=0.01) is True
    assert _four_gamete(3, 47, 2, 48, first_error_rates=result.rates[variant],
                        second_error_rates=result.rates[variant]) is None
    assert one_sided_error_bound(0, 100) == pytest.approx(1 - 0.05 ** .01)


def test_mixed_groups_take_worst_bound_and_deduplicate_ambiguous_fragments():
    variant = Variant("1", 100, "AAAA", "AAA", normalize_contig_names=False)
    controls = [
        _control(variant, "ref", "rg1", "normal-1", 100, 1),
        _control(variant, "alt", "rg1", "mutant-1", 100, 1),
        _control(variant, "ref", "rg2", "normal-2", 90, 10),
        _control(variant, "alt", "rg2", "mutant-2", 80, 20),
    ]
    evidence = controls[0].evidence
    evidence.ref_reads.append(evidence.ref_reads[0])  # alternative placement
    evidence.alt_reads.append(evidence.ref_reads[1])  # conflicting mate/call
    result = calibrate_phasing_error_rates(controls, _target(variant, ("rg1", "rg2")))
    assert result.rates[variant][0] > 0.1
    assert result.rates[variant][1] > 0.2
    assert result.report[0]["ambiguous"] == 1
    assert result.report[0]["callable"] == 100
    assert result.report[0]["confidence"] == pytest.approx(0.9875)
    assert {row["read_group"] for row in result.report} == {"rg1", "rg2"}
    assert _four_gamete(4, 48, 1, 47, phasing_error_rate=0.01) is True
    assert _four_gamete(4, 48, 1, 47, first_error_rates=result.rates[variant],
                        second_error_rates=result.rates[variant]) is None
    assert _four_gamete(1, 3, 33, 63, phasing_error_rate=0.01) is False
    assert _four_gamete(1, 3, 33, 63, first_error_rates=result.rates[variant],
                        second_error_rates=result.rates[variant]) is None


def test_multiple_control_sources_take_worst_source_bound():
    variant = Variant("1", 100, "A", "G", normalize_contig_names=False)
    controls = [
        _control(variant, "ref", "rg1", "normal-low", 100, 0),
        _control(variant, "ref", "rg1", "normal-high", 90, 10),
        _control(variant, "alt", "rg1", "mutant", 100, 0),
    ]
    result = calibrate_phasing_error_rates(controls, _target(variant))
    assert result.rates[variant][0] == max(row["upper_bound"] for row in result.report
                                          if row["truth"] == "ref")
    assert result.rates[variant][0] > 0.1
    assert len(result.report) == 3
    assert all(row["confidence"] == pytest.approx(1 - .05 / 3) for row in result.report)


def test_missing_truth_and_low_coverage():
    variant = Variant("1", 100, "A", "G", normalize_contig_names=False)
    controls = [_control(variant, "ref", "rg1", "normal", 10, 1),
                _control(variant, "alt", "rg1", "mutant", 10, 1)]
    with pytest.raises(ValueError, match="Missing alt control"):
        calibrate_phasing_error_rates(controls[:1], _target(variant))
    sparse = [_control(variant, "ref", "rg1", "normal", 1, 0),
              _control(variant, "alt", "rg1", "mutant", 1, 0)]
    assert min(calibrate_phasing_error_rates(sparse, _target(variant)).rates[variant]) > 0.97


def test_missing_group_context_and_truth_fail_closed():
    variant = Variant("1", 100, "A", "G", normalize_contig_names=False)
    controls = [_control(variant, "ref", "rg1", "normal", 10, 1),
                _control(variant, "alt", "rg1", "mutant", 10, 1)]
    with pytest.raises(ValueError, match="Missing ref control"):
        calibrate_phasing_error_rates(controls, _target(variant, ("rg1", "rg2")))
    with pytest.raises(ValueError, match="differs from target"):
        calibrate_phasing_error_rates(controls, {variant: {"rg1": {**META, "context": "TAC"}}})
    with pytest.raises(ValueError, match="cannot have both"):
        calibrate_phasing_error_rates([controls[0], _control(variant, "alt", "rg1", "normal", 10, 1)],
                                      _target(variant))
    with pytest.raises(ValueError, match="No callable"):
        calibrate_phasing_error_rates([_control(variant, "ref", "rg1", "normal", 0, 0), controls[1]],
                                      _target(variant))
    with pytest.raises(ValueError, match="bound is 1"):
        calibrate_phasing_error_rates([_control(variant, "ref", "rg1", "normal", 0, 10), controls[1]],
                                      _target(variant))


def test_read_groups_and_legacy_reads_are_rejected():
    variant = Variant("1", 100, "A", "G", normalize_contig_names=False)
    control = _control(variant, "ref", "rg1", "normal", 10, 1)
    alt_control = _control(variant, "alt", "rg1", "mutant", 10, 1)
    control.evidence.ref_reads[0].source_alignments = ()
    with pytest.raises(ValueError, match="source identities"):
        calibrate_phasing_error_rates([control, alt_control], _target(variant))


def test_cli_calibrates_from_truth_labeled_bams(tmp_path):
    import pysam

    def write_bam(name, ref_count, alt_count, untagged=False):
        path = tmp_path / (name + ".bam")
        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "1", "LN": 1000}],
                  "RG": [{"ID": group, "PL": "ONT"} for group in ("rg1", "rg2")]}
        with pysam.AlignmentFile(str(path), "wb", header=header) as out:
            for group in ("rg1", "rg2"):
                for allele, count in (("T", ref_count), ("G", alt_count)):
                    for i in range(count):
                        read = pysam.AlignedSegment(out.header)
                        read.query_name = "%s-%s-%s-%d" % (name, group, allele, i)
                        read.query_sequence = "AAA" + allele + "AAAA"
                        read.flag = 0
                        read.reference_id = 0
                        read.reference_start = 0
                        read.mapping_quality = 60
                        read.cigarstring = "8M"
                        read.query_qualities = pysam.qualitystring_to_array("I" * 8)
                        read.set_tag("RG", group)
                        out.write(read)
            if untagged:
                read = pysam.AlignedSegment(out.header)
                read.query_name = name + "-untagged"
                read.query_sequence = "AAATAAAA"
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 0
                read.mapping_quality = 60
                read.cigarstring = "8M"
                read.query_qualities = pysam.qualitystring_to_array("I" * 8)
                out.write(read)
        pysam.index(str(path))
        return str(path)

    ref_bam = write_bam("normal", 30, 5)
    alt_bam = write_bam("mutant", 2, 30)
    target_bam = write_bam("target", 25, 5)
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">1\n" + "AAATAAA" + "A" * 993 + "\n")
    pysam.faidx(str(fasta))
    metadata = {**META, "context": "AAATAAA"}
    manifest = {"reference_fasta": str(fasta), "targets": [{
        "variant": {"contig": "1", "start": 4, "ref": "T", "alt": "G"},
        "context_start": 1,
        "target_bam": target_bam,
        "read_groups": {"rg1": metadata, "rg2": metadata},
        "controls": [
            {"truth": truth, "source": source, "read_group": group, "bam": bam, **metadata}
            for group in ("rg1", "rg2")
            for truth, source, bam in (("ref", "normal", ref_bam),
                                       ("alt", "mutant", alt_bam))
        ],
    }]}
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    profile, report = tmp_path / "profile.json", tmp_path / "report.json"
    run_calibration(["--manifest", str(manifest_path), "--output", str(profile),
                     "--report", str(report)])
    rows = json.loads(profile.read_text())
    audit = json.loads(report.read_text())
    assert rows[0]["ref_to_alt"] > rows[0]["alt_to_ref"]
    assert audit["controls"][0]["wrong"] == 5
    assert audit["controls"][1]["wrong"] == 2
    assert [(row["read_group"], row["wrong"]) for row in audit["controls"]] == [
        ("rg1", 5), ("rg1", 2), ("rg2", 5), ("rg2", 2)]
    assert {row["bam"] for row in audit["control_files"]} == {ref_bam, alt_bam}
    from isovar.cli.isovar_calibrate_phasing import _calibrate
    manifest["targets"][0]["controls"][0]["bam"] = target_bam
    with pytest.raises(ValueError, match="independent"):
        _calibrate(manifest)
    manifest["targets"][0]["controls"][0]["bam"] = ref_bam
    manifest["targets"][0]["read_groups"] = {"other": metadata}
    with pytest.raises(ValueError, match="read groups differ"):
        _calibrate(manifest)
    manifest["targets"][0]["read_groups"] = {"rg1": metadata, "rg2": metadata}
    manifest["targets"][0]["target_bam"] = write_bam("untagged", 25, 5, untagged=True)
    with pytest.raises(ValueError, match="uncalibrated read group"):
        _calibrate(manifest)
    manifest["targets"][0]["read_groups"] = {
        "rg1": {**metadata, "context": "AAAAAAA"}, "rg2": metadata}
    with pytest.raises(ValueError, match="sequence context differs"):
        _calibrate(manifest)
