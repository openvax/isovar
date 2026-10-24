"""Report original-read local SV paths and pinned unrearranged GRCh38 arms."""

import argparse
import csv
import json
from pathlib import Path

from reports.hcc1395_colo829.build_sv_junction_audit import COLO_BAM, EVENTS, HCC_BAM, _reconstruct


HERE = Path(__file__).parent
FLANKS = HERE.parents[1] / 'tests/data/hcc1395_sv_junction/sv_reference_flanks.json'
PREFIX = HERE / 'sv_transcript_evidence'
FIELDS = ('sample', 'call_id', 'variant_class', 'path_id', 'path_bases',
          'direct_junction_fragments', 'direct_junction_fragments_mapq20',
          'junction_12_flank_exact_fragments', 'junction_12_flank_specificity',
          'local_exact_fragments', 'local_interval', 'local_flank_bases',
          'local_specificity', 'local_read_ids', 'local_piece_mapqs',
          'local_missing_quality_reads', 'full_path_exact_fragments',
          'junction_left', 'junction_inserted', 'junction_right',
          'reference_donor_before', 'reference_donor_after',
          'reference_acceptor_before', 'reference_acceptor_after',
          'limitations')


def build(hcc_bam=HCC_BAM, colo_bam=COLO_BAM, flanks=FLANKS):
    """Replay the four original-read events with a direct RNA junction."""
    controls = json.loads(Path(flanks).read_text())
    if controls['schema'] != 'isovar.sv_reference_flanks.v1':
        raise ValueError('Unknown reference flank schema')
    rows, sequences = [], {}
    for event in EVENTS[:4]:
        sample, call_id, variant_class, *_ = event
        bam = hcc_bam if sample == 'HCC1395' else colo_bam
        result = _reconstruct(bam, event, 0)
        subset = _reconstruct(bam, event, 20)
        high_mapq = max((j['direct_support']['fragments'] for p in subset['paths']
                         for j in p['junctions'] if j['relation'] == 'breakpoint_junction'), default=0)
        control = controls['events'][call_id]
        for path in result['paths']:
            junction = next((j for j in path['junctions']
                             if j['relation'] == 'breakpoint_junction'), None)
            if junction is None:
                continue
            if junction['left'] != control['left'] or junction['right'] != control['right']:
                raise ValueError('Observed junction differs from pinned GRCh38 flank')
            start, end = junction['query_interval']
            local = path['event_local_support']
            observation_keys = {witness['observation'] for witness in local['witnesses']}
            observed = [result['observations'][key] for key in observation_keys]
            record_ids = {rid for observation in observed for rid in observation['records']}
            row = dict(sample=sample, call_id=call_id, variant_class=variant_class,
                       path_id=path['path_id'], path_bases=len(path['sequence']),
                       direct_junction_fragments=junction['direct_support']['fragments'],
                       direct_junction_fragments_mapq20=high_mapq,
                       junction_12_flank_exact_fragments=path['event_junction_window_support']['fragments'],
                       junction_12_flank_specificity=path['event_junction_window_support']['specificity'],
                       local_exact_fragments=local['fragments'],
                       local_interval='%d:%d' % tuple(local['query_interval']),
                       local_flank_bases=local['flank_bases'],
                       local_specificity=local['specificity'],
                       local_read_ids=';'.join(sorted({observation['identity'][1] for observation in observed})),
                       local_piece_mapqs=';'.join(map(str, sorted({
                           result['record_evidence'][rid]['mapping_quality_raw'] for rid in record_ids}))),
                       local_missing_quality_reads=local['missing_quality_reads'],
                       full_path_exact_fragments=path['full_interval_support']['fragments'],
                       junction_left=path['sequence'][max(0, start - 20):start],
                       junction_inserted=junction['unplaced_bases'],
                       junction_right=path['sequence'][end:end + 20],
                       reference_donor_before=control['donor_before'],
                       reference_donor_after=control['donor_after'],
                       reference_acceptor_before=control['acceptor_before'],
                       reference_acceptor_after=control['acceptor_after'],
                       limitations=';'.join(result['limitations']))
            rows.append(row)
            sequences[path['path_id']] = path['sequence']
    return rows, sequences, controls


def write(prefix=PREFIX, **kwargs):
    """Write machine-readable rows, full path sequences, and an assessment."""
    rows, sequences, controls = build(**kwargs)
    prefix = Path(prefix)
    with prefix.with_suffix('.tsv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    with prefix.with_suffix('.fasta').open('w') as stream:
        for row in rows:
            sequence = sequences[row['path_id']]
            stream.write('>%s %s %s %s\n' % (
                row['path_id'], row['sample'], row['call_id'], row['variant_class']))
            for start in range(0, len(sequence), 80):
                stream.write(sequence[start:start + 80] + '\n')
    counts = {sample: sum(row['sample'] == sample for row in rows)
              for sample in ('HCC1395', 'COLO829')}
    lines = [
        '# Original RNA evidence for selected SV transcript paths', '',
        'These are **bounded, breakpoint-compatible RNA path hypotheses**, not a',
        'complete atlas of SV transcripts. The selected original-read fixtures cover',
        'three HCC1395 and one COLO829 DNA events with direct junction observations.',
        f'They yield **{len(rows)} retained paths** ({counts["HCC1395"]} HCC1395;',
        f'{counts["COLO829"]} COLO829). Only **{sum(row["local_exact_fragments"] >= 2 for row in rows)}**',
        'have at least two exact local path witnesses. HCC1395 chr9 reaches',
        '`path_limit`, so this is a retained-hypothesis count rather than a',
        'count of biological isoforms.', '',
        '**Evidence units:** direct fragments are distinct original RNA fragments',
        'making the join, shared by paths. Local exact fragments each match one',
        'reported window across the join with compatible genomic placements.',
        'The 12-base-flank window is reported separately when it cannot',
        'distinguish competing paths. A `unique_among_retained_paths` window',
        'differs in sequence from other retained event paths; uniqueness is',
        'limited by the path limit and does',
        'not establish a full isoform. Exact matching can undercount noisy ONT',
        'reads; a path distinguished by a single ONT base can also be a read',
        'error rather than a biological RNA branch.',
        'The [TSV](sv_transcript_evidence.tsv) retains original local witness',
        'read IDs, raw alignment MAPQs, quality availability, intervals, and',
        'assembly limitations. The collector allows missing base qualities;',
        'no base-quality minimum is imposed on these counts.',
        'Full-path exact support requires one read to match the entire candidate,',
        'which is a much stricter observation than assembly from overlapping reads.',
        'MAPQ 20 is a sensitivity subset, not a gate for the primary SV search.', '',
        '| Sample / event | Path | RNA bases | Direct join (all / MAPQ20) | 12-flank exact | Selected local exact / window | Full exact | Limits |',
        '|---|---|---:|---:|---:|---:|---:|---|',
    ]
    for row in rows:
        label = row['call_id'].replace('COLO829_SOMATIC_', '')
        lines.append('| %s %s | `%s` | %d | %d / %d | %d (%s) | %d / %d flank (%s) | %d | %s |' % (
            row['sample'], label, row['path_id'][:8], row['path_bases'],
            row['direct_junction_fragments'], row['direct_junction_fragments_mapq20'],
            row['junction_12_flank_exact_fragments'],
            {'shared_or_unresolved': 'shared',
             'unique_among_retained_paths': 'unique',
             'no_competing_retained_path': 'one retained'}[row['junction_12_flank_specificity']],
            row['local_exact_fragments'], row['local_flank_bases'],
            {'unique_among_retained_paths': 'unique',
             'no_competing_retained_path': 'one retained',
             'shared_or_unresolved': 'shared'}[row['local_specificity']],
            row['full_path_exact_fragments'], row['limitations'].replace(';', ', ') or 'none'))
    lines += ['', '## Observed join versus unrearranged reference', '',
              'The RNA column is `20 observed upstream bases | unplaced join bases |',
              '20 observed downstream bases`. Reference columns are oriented genomic',
              'controls at the observed left and right anchors: donor continuation',
              'and acceptor predecessor, respectively. These are **genomic arms**,',
              'not asserted mature reference transcript sequences. The full assembled',
              'candidate sequences are in [`sv_transcript_evidence.fasta`](sv_transcript_evidence.fasta).', '',
              '| Event / path | RNA join | Reference donor arm | Reference acceptor arm |',
              '|---|---|---|---|']
    for row in rows:
        lines.append('| %s `%s` | `%s|%s|%s` | `%s|%s` | `%s|%s` |' % (
            row['call_id'], row['path_id'][:8], row['junction_left'],
            row['junction_inserted'] or '∅', row['junction_right'],
            row['reference_donor_before'], row['reference_donor_after'],
            row['reference_acceptor_before'], row['reference_acceptor_after']))
    lines += ['', 'The GRCh38 controls are pinned in',
              '[`sv_reference_flanks.json`](../../tests/data/hcc1395_sv_junction/sv_reference_flanks.json)',
              'from the [NCBI no-alt analysis set](%s), uncompressed FASTA SHA-256 `%s`.' % (
                  controls['source_url'], controls['source_fasta_sha256']), '',
              '## Remaining SV screen signals', '',
              'The HCC1395 whole-run screen has six nominated SV records with',
              'split-alignment candidates. The three outside this table are',
              'chr1 records associated with POU2F1–NME7 or nearby complex',
              'breaks; all three cite the **same single ONT read**',
              '`SRR23456709.320596`. Its selected fixture contains two actual',
              'alignment records and SA tags pointing to another placement;',
              'an SA hint alone cannot reconstruct that missing piece. These',
              'are one-read hypotheses, not three independent supported RNA',
              'transcripts. The remaining COLO829 SV screen records have mate',
              'context or breakend coverage but no independently observed',
              'junction from which this method can assemble a path. The',
              '[HCC1395](HCC1395_ATLAS.md) and [COLO829](COLO829_ATLAS.md)',
              'source atlases enumerate the larger screened',
              'call sets; this table covers all directly replayed joins with',
              'at least two original fragments in the pinned selected-site',
              'fixtures.', '',
              '## Stronger fusion path from a separate HCC1395 source', '',
              'The EIF3K–CYP39A1 fusion-enriched INTEGRATE RNA source supports one',
              '531-aa predicted proteoform with 149 distinct Illumina junction',
              'fragments, 19 tiling Illumina reads covering the 1,596-nt coding path,',
              'and four independent ONT reads with exact junction sequence on',
              'collinear start-to-stop alignments. Matched tumor DNA has 18 reciprocal',
              'pairs between flanking introns versus zero normal pairs. These are',
              'separate data and cannot be added to the GRCh38 selected-site counts.',
              'See the [fusion report](eif3k_cyp39a1/report.md) and',
              '[protein FASTA](eif3k_cyp39a1/candidate_proteins.fasta).', '',
              '## What SV3 still lacks', '',
              'SV3 has three original split-junction fragments at all MAPQs and',
              'three exact 12-base-flank local witnesses, but zero MAPQ20 join',
              'fragments and zero whole-path witnesses. Six opposite-mate pairs',
              'are context, not split-junction support. The 1,213-base candidate',
              'reaches `extension_segment_limit`, with no observed mature splice',
              'connection to a CNIH3 exon. Its locus is intronic under Ensembl 87',
              'and lacks a bundled ORFdata coding block. Two ATG-to-stop translations',
              '(37 and 162 aa) share the same three nucleotide join witnesses;',
              'neither has an established initiation site, reading frame, or full',
              'ORF witness. A separate 68,553-representative polished isoform',
              'check found no exact SV3 27-mer. Thus the RNA join is observed,',
              'while a mature transcript, somatic causation and proteoform remain',
              'unresolved. See the [SV3 ORF assessment](colo829_sv3_candidate_orfs.assessment.md).', '',
              'The original BAMs are selected-site fixtures, not all HCC1395 or',
              'COLO829 RNA reads. The HCC and COLO DNA/RNA sources also come from',
              'different cell cultures or passages. The [HCC1395](HCC1395_ATLAS.md)',
              'and [COLO829](COLO829_ATLAS.md) atlases give the screened source',
              'universe and events without direct reconstructed joins.', '',
              'Regenerate the table from pinned inputs with',
              '`python -m reports.hcc1395_colo829.build_sv_transcript_evidence`.',
              'The exact-window policy and its known ONT error sensitivity are',
              'specified in `SV_TRANSCRIPT_EVIDENCE_SPEC.md`.', '']
    lines += ['The [versioned exon path audit](sv_exon_path_evidence.md) is generated',
              'by the same command and keeps exon connection and frame status',
              'separate from local join evidence.', '']
    prefix.with_suffix('.md').write_text('\n'.join(lines))
    from reports.hcc1395_colo829.build_sv_exon_path_evidence import write as write_exon

    write_exon(prefix.with_name("sv_exon_path_evidence"),
               hcc_bam=kwargs.get("hcc_bam", HCC_BAM),
               colo_bam=kwargs.get("colo_bam", COLO_BAM))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hcc-bam', type=Path, default=HCC_BAM)
    parser.add_argument('--colo-bam', type=Path, default=COLO_BAM)
    parser.add_argument('--reference-flanks', type=Path, default=FLANKS)
    parser.add_argument('--output-prefix', type=Path, default=PREFIX)
    args = parser.parse_args()
    write(args.output_prefix, hcc_bam=args.hcc_bam, colo_bam=args.colo_bam,
          flanks=args.reference_flanks)


if __name__ == '__main__':
    main()
