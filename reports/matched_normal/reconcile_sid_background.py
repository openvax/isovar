"""Compare a new normal receipt with Sid's prior supplied background and flags.

This report consumes library evidence; it never changes completed RNA results.
All removals are described under the new scope/policy, not automatically called
errors. A separate original-VCF audit identifies the earlier multiallelic-index
bug within the acquired territory.
"""

import argparse
from collections import Counter
from collections.abc import Mapping, Sequence
import csv
from functools import lru_cache
from pathlib import Path
from typing import Any

import pysam
from varcode import normalize_contig_name

from isovar.genomic_intervals import GenomicIntervalIndex
from isovar.normal_evidence import literal_allele_key, reconcile_normal_states
from isovar.normal_sources import AlleleKey, MatchedNormalEvidence
from isovar.sv_rna_audit import file_sha256, read_json, write_json


def allele_key(row: Mapping[str, Any]) -> AlleleKey:
    """Use the library's trimmed, assembly-aware literal identity for both inputs."""
    return literal_allele_key(row['contig'], row['start'], row['ref'], row['alt'], reference_name='GRCh38')


@lru_cache(maxsize=2048)
def canonical_contig(name: str) -> str:
    """Normalize each source vocabulary once while unioning large region lists."""
    return normalize_contig_name(name, reference_name='GRCh38')


def legacy_unindexed_alternates(sources: Sequence[Mapping[str, Any]], prior_keys: set[AlleleKey]
                               ) -> list[dict[str, Any]]:
    """Identify listed ALT alleles without a corresponding called GT index."""
    called, uncalled = set(), {}
    for source in sources:
        with pysam.VariantFile(source['path'], index_filename=source['index_path']) as reader:
            for record in reader:
                if set(record.filter) - {'PASS'}:
                    continue
                genotype = record.samples[source['vcf_sample']].get('GT') or ()
                indexes = {i for i in genotype if i is not None and i > 0}
                if not indexes:
                    continue
                for index, alt in enumerate(record.alts or (), 1):
                    if set(record.ref + alt) - set('ACGT'):
                        continue
                    key = literal_allele_key(record.contig, record.pos, record.ref, alt, reference_name='GRCh38')
                    if key not in prior_keys:
                        continue
                    if index in indexes:
                        called.add(key)
                    elif key in prior_keys:
                        uncalled.setdefault(key, []).append(dict(source_id=source['source_id'], contig=record.contig,
                            position=record.pos, ref=record.ref, alts=list(record.alts), genotype=list(genotype),
                            uncalled_alt_index=index))
    return [dict(allele=list(key), records=uncalled[key]) for key in sorted(set(uncalled) - called)]


def compare_background(assessment: MatchedNormalEvidence, old_background: Sequence[Mapping[str, Any]],
                       selected: Sequence[Mapping[str, Any]]) -> tuple[dict[str, Any], dict[str, Any], set[str]]:
    """Compare genotype-indexed backgrounds within declared requested/acquired scope."""
    old_keys = {allele_key(row) for row in old_background}
    baseline_ids = {source_id for row in old_background for source_id in row['origins']}
    baseline_sources = [row for row in selected if row['source_id'] in baseline_ids]
    if {row['source_id'] for row in baseline_sources} != baseline_ids:
        raise ValueError('The new assessment must include every previously supplied normal source')
    acquired = GenomicIntervalIndex({(canonical_contig(contig), start, end)
        for source in selected for contig, start, end in source['acquired_regions']})
    old_in_scope = {allele_key(row) for row in old_background if acquired.covers(
        canonical_contig(row['contig']), row['start'] - 1,
        row['start'] - 1 + max(1, len(row['ref'])))}
    current_keys = {allele_key(row) for row in assessment['background_variants']}
    baseline_current = {allele_key(row) for row in assessment['background_variants'] if any(
        origin['source_id'] in baseline_ids for origin in row['origins'])}
    requested = GenomicIntervalIndex((canonical_contig(contig), start, end)
        for contig, start, end in assessment.get('background_regions', ()))
    previous_in_requested = sum(requested.covers(canonical_contig(row['contig']),
        row['start'] - 1, row['start'] - 1 + max(1, len(row['ref']))) for row in old_background)
    unindexed = legacy_unindexed_alternates(baseline_sources, old_keys)
    summary = dict(previous_total=len(old_keys), previous_in_acquired_scope=len(old_in_scope),
        previous_in_requested_scope=previous_in_requested if assessment.get('background_regions') else None,
        current_total=len(current_keys), current_from_baseline_sources=len(baseline_current),
        current_only_from_additional_sources=len(current_keys - baseline_current),
        newly_supplied=len(current_keys - old_keys), previous_in_scope_not_retained=len(old_in_scope - current_keys),
        legacy_listed_but_uncalled=len(unindexed))
    artifacts = {
        'new_background_alleles.json.gz': [list(key) for key in sorted(current_keys - old_keys)],
        'previous_background_not_retained_in_scope.json.gz': [list(key) for key in sorted(old_in_scope - current_keys)],
        'legacy_unindexed_alternates.json.gz': unindexed,
    }
    return summary, artifacts, baseline_ids


def compare_nominations(assessment: MatchedNormalEvidence, old_variants: Mapping[AlleleKey, Mapping[str, Any]],
                        selected: Sequence[Mapping[str, Any]], baseline_ids: set[str]
                        ) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    """Separate changed origin-review flags from new RNA or protein reconstruction."""
    transitions, states = Counter(), Counter()
    source_states = {row['source_id']: Counter() for row in selected}
    rejections = {row['source_id']: Counter() for row in selected}
    changes = []
    old_conflicts, new_conflicts = 0, 0
    mitochondrial = []
    for row in assessment['alleles']:
        key = tuple(row['normalized_allele'])
        old = old_variants[key]
        old_flag, new_flag = bool(old.get('origin_conflict')), bool(row['origin_conflict'])
        old_conflicts += old_flag
        new_conflicts += new_flag
        states[row['state']] += 1
        baseline = reconcile_normal_states(outcome['state'] for outcome in row['source_assessments']
                                           if outcome['source_id'] in baseline_ids)
        transitions[(baseline, row['state'])] += 1
        if old_flag != new_flag:
            changes.append(dict(contig=row['contig'], start=row['start'], ref=row['ref'], alt=row['alt'],
                previous_origin_conflict=old_flag, new_origin_conflict=new_flag, new_state=row['state'],
                baseline_same_policy_state=baseline))
        for outcome in row['source_assessments']:
            source_states[outcome['source_id']][outcome['state']] += 1
            for record in outcome.get('records', []):
                rejections[outcome['source_id']].update(record['rejection_reasons'])
        if key[0] == 'MT':
            mitochondrial.append({field: value for field, value in row.items() if field != 'nomination'})
    if len(old_variants) != len(assessment['alleles']):
        raise ValueError('Previous and new nomination inventories differ')
    summary = dict(nominated_alleles=len(assessment['alleles']), normal_states=dict(states),
        previous_origin_conflicts=old_conflicts, new_origin_conflicts=new_conflicts,
        newly_flagged=sum(row['new_origin_conflict'] for row in changes),
        no_longer_flagged=sum(row['previous_origin_conflict'] for row in changes),
        same_policy_state_transitions=[dict(baseline=old, new=new, alleles=count) for (old, new), count in sorted(transitions.items())])
    per_source = [dict(source_id=source['source_id'], specimen_id=source['specimen_id'], library_id=source['library_id'],
        caller=source['caller'], evidence_kind=source['evidence_kind'], policy=source['policy'],
        states=dict(source_states[source['source_id']]),
        rejected_record_query_reasons=dict(rejections[source['source_id']])) for source in selected]
    return summary, {'origin_changes.json.gz': changes, 'mitochondrial_evidence.json.gz': mitochondrial}, per_source


def reconcile(assessment_path: Path, previous_variants: Path, previous_background: Path,
              output: Path) -> dict[str, Any]:
    """Write source summaries and every changed nomination/background identity."""
    assessment = read_json(assessment_path)
    if assessment['schema'] != 'isovar.matched_normal.v1' or assessment['sample_id'] != 'Sid':
        raise ValueError('Expected the versioned Sid matched-normal assessment')
    old_variants = {allele_key(row): row for row in read_json(previous_variants)}
    selected = [row for row in assessment['source_inventory'] if row['selected']]
    background, artifacts, baseline_ids = compare_background(assessment, read_json(previous_background), selected)
    nominations, nomination_artifacts, per_source = compare_nominations(assessment, old_variants, selected, baseline_ids)
    artifacts.update(nomination_artifacts)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if any((output / name).exists() for name in (*artifacts, 'reconciliation.json')):
        raise ValueError('Reconciliation already exists; choose a new output directory')
    for name, data in artifacts.items():
        write_json(output / name, data)
    summary = dict(schema='isovar.normal_reconciliation.v1', sample_id='Sid',
        inputs={name: dict(path=str(path), sha256=file_sha256(path)) for name, path in (
            ('assessment', assessment_path), ('previous_variants', previous_variants), ('previous_background', previous_background))},
        software=assessment['software'], background_flank=assessment['background_flank'],
        background_region_count=len(assessment.get('background_regions', ())),
        background_source_scopes=assessment.get('background_source_scopes', []),
        selected_sources=len(selected), inventory_sources=len(assessment['source_inventory']),
        specimens=sorted({s['specimen_id'] for s in selected}), per_source=per_source,
        baseline_source_ids=sorted(baseline_ids), background=background, **nominations,
        artifacts={name: dict(sha256=file_sha256(output / name), entries=len(data)) for name, data in artifacts.items()},
        interpretation='Versioned DNA-normal/background reassessment. Earlier RNA/protein receipts retain their original supplied scope. '
            'A changed origin flag is not a new RNA reconstruction. Rejection counts are record-query observations and may overlap. '
            'Background comparison is limited to declared acquired intervals; absence is never tumor specificity.')
    write_json(output / 'reconciliation.json', summary)
    with (output / 'source_inventory.tsv').open('w', newline='') as stream:
        columns = ['source_id', 'selected', 'specimen_id', 'library_id', 'vcf_sample', 'caller', 'evidence_kind', 'selection_reason', 'url']
        writer = csv.DictWriter(stream, columns, delimiter='\t')
        writer.writeheader()
        for source in assessment['source_inventory']:
            writer.writerow({**{key: source.get(key) for key in columns[:-1]}, 'url': source['provenance']['url']})
    print({key: summary[key] for key in ('normal_states', 'previous_origin_conflicts', 'new_origin_conflicts', 'newly_flagged', 'no_longer_flagged', 'background')})
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assessment', type=Path, required=True)
    parser.add_argument('--previous-variants', type=Path, required=True)
    parser.add_argument('--previous-background', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    options = parser.parse_args()
    reconcile(options.assessment, options.previous_variants, options.previous_background, options.output_dir)
