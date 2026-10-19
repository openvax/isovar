"""Derive per-locus phasing error bounds from truth-labeled control BAMs."""

import argparse
from contextlib import ExitStack
from hashlib import sha256
import json
from pathlib import Path

from pysam import AlignmentFile, FastaFile
from varcode import Variant

from .. import __version__
from ..phasing_calibration import PhasingControl, calibrate_phasing_error_rates
from ..read_collector import ReadCollector
from ..read_evidence import ReadEvidence
from ..read_identity import source_read_ids

_COLLECTOR_OPTIONS = frozenset((
    "use_secondary_alignments", "use_duplicate_reads", "min_mapping_quality",
    "use_soft_clipped_bases", "merge_overlapping_fragments",
    "use_reads_without_base_qualities",
))


def _read_groups(bam):
    return {entry["ID"]: entry for entry in bam.header.to_dict().get("RG", ())}


def _check_group(bam, group, metadata):
    header_group = _read_groups(bam).get(group)
    if header_group is None:
        raise ValueError("Read group %s is absent from %s" % (group, bam.filename))
    if header_group.get("PL", "").upper() != metadata["platform"].upper():
        raise ValueError("Read group %s platform differs from BAM header" % group)


def _group_evidence(evidence, group=None):
    buckets, observed = [], set()
    for allele in ("ref", "alt", "other"):
        selected = []
        for read in getattr(evidence, allele + "_reads"):
            groups = {key[0] for key in source_read_ids(read)}
            if len(groups) != 1:
                raise ValueError("Collected reads need one unambiguous read group")
            observed.update(groups)
            if group is None or group in groups:
                selected.append(read)
        buckets.append(selected)
    if group is None:
        return observed
    return ReadEvidence(evidence.trimmed_base1_start, evidence.trimmed_ref,
                        evidence.trimmed_alt, *buckets)


def _calibrate(manifest, base_path=Path(".")):
    if not isinstance(manifest, dict) or not isinstance(manifest.get("targets"), list):
        raise ValueError("Manifest needs a targets list")
    collector_options = manifest.get("read_collector_options", {})
    if not isinstance(collector_options, dict) or set(collector_options) - _COLLECTOR_OPTIONS:
        raise ValueError("Unsupported read_collector_options")
    for key, value in collector_options.items():
        if key == "min_mapping_quality":
            if type(value) is not int or value < 0:
                raise ValueError("min_mapping_quality must be a nonnegative integer")
        elif type(value) is not bool:
            raise ValueError("%s must be boolean" % key)
    collector = ReadCollector(**collector_options)
    target_metadata, controls = {}, []
    with ExitStack() as stack:
        reference_path = str((base_path / manifest["reference_fasta"]).resolve())
        reference_fasta = stack.enter_context(FastaFile(reference_path))
        bams, evidence_cache = {}, {}
        for target in manifest["targets"]:
            locus = target["variant"]
            variant = Variant(locus["contig"], locus["start"], locus["ref"],
                              locus["alt"], normalize_contig_names=False)
            if variant in target_metadata:
                raise ValueError("Duplicate target variant")
            target_metadata[variant] = target["read_groups"]
            context_start = target["context_start"]
            if type(context_start) is not int or context_start < 1:
                raise ValueError("context_start must be a positive 1-based position")
            for metadata in target["read_groups"].values():
                context = metadata["context"]
                context_end = context_start + len(context) - 1
                variant_end = variant.start + max(len(variant.ref), 1) - 1
                if variant.start < context_start or variant_end > context_end:
                    raise ValueError("Sequence context must span the target variant")
                observed = reference_fasta.fetch(variant.contig, context_start - 1,
                                                 context_end).upper()
                if observed != context.upper():
                    raise ValueError("Declared sequence context differs from reference FASTA")
                reference_allele = reference_fasta.fetch(
                    variant.contig, variant.start - 1, variant.start - 1 + len(variant.ref)).upper()
                if reference_allele != variant.ref.upper():
                    raise ValueError("Target reference allele differs from reference FASTA")
            target_path = str((base_path / target["target_bam"]).resolve())
            if target_path not in bams:
                bams[target_path] = stack.enter_context(AlignmentFile(target_path))
            target_bam = bams[target_path]
            if set(_read_groups(target_bam)) != set(target["read_groups"]):
                raise ValueError("Target BAM read groups differ from the manifest")
            for group, metadata in target["read_groups"].items():
                _check_group(target_bam, group, metadata)
            target_sq = {sq["SN"]: sq for sq in target_bam.header.to_dict()["SQ"]}
            if variant.contig not in target_sq:
                raise ValueError("Target BAM has no sequence for %s" % variant.contig)
            if target_sq[variant.contig]["LN"] != reference_fasta.get_reference_length(variant.contig):
                raise ValueError("Target BAM and reference FASTA contig lengths differ")
            target_evidence = collector.read_evidence_for_variant(variant, target_bam)
            if _group_evidence(target_evidence) - set(target["read_groups"]):
                raise ValueError("Target reads contain an uncalibrated read group")
            for entry in target["controls"]:
                path = str((base_path / entry["bam"]).resolve())
                if path == target_path:
                    raise ValueError("A control BAM must be independent of the target BAM")
                if path not in bams:
                    bams[path] = stack.enter_context(AlignmentFile(path))
                control_bam = bams[path]
                _check_group(control_bam, entry["read_group"], entry)
                control_sq = {sq["SN"]: sq for sq in control_bam.header.to_dict()["SQ"]}
                if variant.contig not in control_sq:
                    raise ValueError("Control BAM has no sequence for %s" % variant.contig)
                target_contig, control_contig = target_sq[variant.contig], control_sq[variant.contig]
                if target_contig["LN"] != control_contig["LN"] or (
                        target_contig.get("M5") != control_contig.get("M5")):
                    raise ValueError("Control and target BAM reference sequences differ")
                evidence_key = (path, variant)
                if evidence_key not in evidence_cache:
                    evidence_cache[evidence_key] = collector.read_evidence_for_variant(
                        variant, control_bam)
                    if _group_evidence(evidence_cache[evidence_key]) - set(_read_groups(control_bam)):
                        raise ValueError("Control reads contain an undeclared read group")
                evidence = _group_evidence(evidence_cache[evidence_key], entry["read_group"])
                controls.append(PhasingControl(
                    variant=variant, truth=entry["truth"],
                    read_group=entry["read_group"], source=entry["source"],
                    platform=entry["platform"], basecaller=entry["basecaller"],
                    library=entry["library"],
                    biological_context=entry["biological_context"],
                    reference=entry["reference"], context=entry["context"],
                    evidence=evidence))
    return calibrate_phasing_error_rates(
        controls, target_metadata, confidence=manifest.get("confidence", 0.95))


def run(args=None, *, prog=None):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path,
                        help="JSON target loci and independent true-ref/true-alt control BAMs")
    parser.add_argument("--output", required=True, type=Path,
                        help="JSON profile accepted by isovar --phasing-error-rates")
    parser.add_argument("--report", required=True, type=Path,
                        help="JSON audit report with per-group counts and confidence bounds")
    options = parser.parse_args(args)
    if options.output == options.report:
        parser.error("--output and --report must be different files")
    try:
        manifest_bytes = options.manifest.read_bytes()
        manifest = json.loads(manifest_bytes)
        calibration = _calibrate(manifest, options.manifest.resolve().parent)
        with options.output.open("w") as handle:
            json.dump(calibration.profile_rows(), handle, indent=2)
            handle.write("\n")
        files = [dict(contig=target["variant"]["contig"], start=target["variant"]["start"],
                      ref=target["variant"]["ref"], alt=target["variant"]["alt"],
                      source=entry["source"], truth=entry["truth"],
                      read_group=entry["read_group"],
                      bam=str((options.manifest.resolve().parent / entry["bam"]).resolve()))
                 for target in manifest["targets"] for entry in target["controls"]]
        with options.report.open("w") as handle:
            json.dump({"isovar_version": __version__,
                       "manifest_sha256": sha256(manifest_bytes).hexdigest(),
                       "confidence": calibration.confidence,
                       "read_collector_options": manifest.get("read_collector_options", {}),
                       "reference_fasta": str((options.manifest.resolve().parent /
                                                manifest["reference_fasta"]).resolve()),
                       "target_files": [dict(**target["variant"],
                                             bam=str((options.manifest.resolve().parent /
                                                      target["target_bam"]).resolve()))
                                        for target in manifest["targets"]],
                       "control_files": files, "controls": calibration.report}, handle, indent=2)
            handle.write("\n")
    except (OSError, KeyError, TypeError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    run()
