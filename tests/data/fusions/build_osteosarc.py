"""Extract observed, directly spanning windows; never synthesize RNA bases."""
import argparse
import collections
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pysam
from pyensembl import EnsemblRelease
from isovar.chimeric_alignment import compatible_phasing_alignments, source_alignment_paths_from_pysam
from isovar.read_identity import source_alignments_from_pysam
from isovar.default_parameters import USE_READS_WITHOUT_BASE_QUALITIES
from tests.osteosarc_protein_helpers import reverse_complement

EVENTS = [('TPST1--CRCP', 'chr7', 66205522, 'chr7', 66127704, ['TPST1'], ['CRCP']),
          ('FOXO3--STRADA-CCDC47', 'chr6', 108611245, 'chr17', 63745980, ['FOXO3'], ['STRADA', 'CCDC47']),
          ('PARD3B--CDKN2B-AS1-CDKN2B', 'chr2', 205310103, 'chr9', 22007648, ['PARD3B'], ['CDKN2B-AS1', 'CDKN2B'])]


def blocks(mapping, lo, hi):
    result = []
    for q in range(lo, hi):
        if q not in mapping:
            continue
        chrom, pos, strand = mapping[q]
        contiguous = result and (result[-1]['reference_end'] == pos if strand == '+'
                                 else result[-1]['reference_start'] == pos + 1)
        if result and result[-1]['contig'] == chrom and result[-1]['strand'] == strand and result[-1]['query_end'] == q - lo and contiguous:
            result[-1]['query_end'] += 1
            result[-1]['reference_end'] = max(result[-1]['reference_end'], pos + 1)
            result[-1]['reference_start'] = min(result[-1]['reference_start'], pos)
        else:
            result.append(dict(query_start=q-lo, query_end=q-lo+1, contig=chrom,
                               reference_start=pos, reference_end=pos+1, strand=strand))
    return result


def references(g, genes):
    models = []
    for gene in genes:
        for tid in g.transcript_ids_of_gene_name(gene):
            t = g.transcript_by_id(tid)
            if t.sequence is None or set(t.sequence) - set('ACGT'):
                continue
            start = min(t.start_codon_spliced_offsets) if t.complete else None
            end = max(t.stop_codon_spliced_offsets)+1 if t.complete else None
            models.append(dict(transcript_id=tid, reference_name='GRCh38', annotation='Ensembl 87',
                               contig='chr'+t.contig, strand=t.strand, exons=sorted((a-1,b) for a,b in t.exon_intervals),
                               sequence=t.sequence, cds_start=start, cds_end=end))
    return models


def query_mapping(read):
    """Actual aligned bases in original-query coordinates, including hard clips.

    SA CIGARs are approximate in minimap2 and MUST NOT supply base mappings.
    """
    hard = read.cigartuples[0][1] if read.cigartuples[0][0] == 5 else 0
    length = read.infer_read_length()
    return {(length - 1 - q - hard if read.is_reverse else q + hard):
            (read.reference_name, r, '-' if read.is_reverse else '+')
            for q, r in read.get_aligned_pairs(matches_only=True)}


def segment_key(read):
    return (read.get_tag('RG') if read.has_tag('RG') else '', read.query_name,
            read.flag & 0xc0 if read.is_paired else 0)


def linked_records(first, second, records):
    """Validate observed supplementary pieces, never paired mates/alternatives."""
    if segment_key(first) != segment_key(second) or first is second:
        return False
    placements = collections.defaultdict(set)
    observations = []
    for read in records:
        for segment, placement in source_alignments_from_pysam(read, read.query_name):
            placements[segment].add(placement)
    for read in (first, second):
        source = source_alignments_from_pysam(read, read.query_name)
        pairs = read.get_aligned_pairs(matches_only=True)
        q, _ = pairs[len(pairs) // 2]
        paths = source_alignment_paths_from_pysam(read, source, (q, q + 1))
        observations.append(SimpleNamespace(source_alignments=source, source_alignment_paths=paths))
    return compatible_phasing_alignments(dict(observations[0].source_alignments),
                                        *observations, placements)


def extract(read, c1, p1, c2, p2, flank, records=(), strands=('+', '+'),
            use_reads_without_base_qualities=USE_READS_WITHOUT_BASE_QUALITIES):
    if read.flag & (4 | 256 | 512 | 1024 | 2048) or not read.has_tag('SA') or read.mapping_quality < 20:
        return None
    sequence = read.get_forward_sequence()
    if sequence is None or any(op == 5 for op, n in read.cigartuples):
        return None
    candidates = []
    for alignment in records:
        if (segment_key(alignment) != segment_key(read) or alignment.flag & (4 | 256 | 512 | 1024)
                or alignment.mapping_quality < 20 or alignment.infer_read_length() != len(sequence)):
            continue
        observed = alignment.get_forward_sequence()
        if observed is not None:
            clips = alignment.cigartuples[::-1] if alignment.is_reverse else alignment.cigartuples
            offset = clips[0][1] if clips[0][0] == 5 else 0
            if sequence[offset:offset + len(observed)] != observed:
                return None
        candidates.append((alignment, query_mapping(alignment)))
    donors = [(r, m, q) for r, m in candidates for q,p in m.items() if p[:2] == (c1, p1-1)]
    acceptors = [(r, m, q) for r, m in candidates for q,p in m.items() if p[:2] == (c2, p2-1)]
    if len(donors) != 1 or len(acceptors) != 1:
        return None
    first, donor, q0 = donors[0]
    second, acceptor, q1 = acceptors[0]
    if not linked_records(first, second, records):
        return None
    reverse = donor[q0][2] != strands[0]
    quality = read.get_forward_qualities() if read.query_qualities is not None else None
    if reverse:
        sequence = reverse_complement(sequence)
        quality = quality[::-1] if quality is not None else None
        def reverse_mapping(mapping):
            return {len(sequence)-q-1:(c,p,'-' if s=='+' else '+') for q,(c,p,s) in mapping.items()}
        donor, acceptor = reverse_mapping(donor), reverse_mapping(acceptor)
        q0, q1 = len(sequence)-q0-1, len(sequence)-q1-1
    if acceptor[q1][2] != strands[1]:
        return None
    q0 += 1
    if not q0 <= q1 or q1-q0 > 100:
        return None
    lo, hi = q0-flank, q1+flank
    if lo < 0 or hi > len(sequence):
        return None
    if ((quality is None and not use_reads_without_base_qualities) or
            (quality is not None and min(quality[lo:hi]) < 10)):
        return None
    mapping = {q:p for q,p in donor.items() if lo <= q < q0}
    mapping.update({q:p for q,p in acceptor.items() if q1 <= q < hi})
    local = sequence[lo:hi]
    if set(local)-set('ACGT'):
        return None
    return dict(sequence=local, junction_start=flank, junction_end=q1-lo,
                blocks=blocks(mapping,lo,hi), source_query_start=lo,
                minimum_base_quality=min(quality[lo:hi]) if quality is not None else None,
                base_quality_status='available' if quality is not None else 'missing', read_id=read.query_name,
                source_reverse_complement=reverse != read.is_reverse,
                fragment_id=(read.get_tag('CB') + '/' + read.get_tag('UB')
                             if read.has_tag('CB') and read.has_tag('UB') else read.query_name),
                original_sam=read.to_string(), partner_sam=second.to_string() if first is read else first.to_string())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--alignments', type=Path, required=True,
                        help='Acquired regional BAM directory, with source-id/regions.bam and receipt.json')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    g=EnsemblRelease(87)
    summaries=[]
    for name,c1,p1,c2,p2,left_genes,right_genes in EVENTS:
        refs=references(g,left_genes+right_genes)
        for sid,sample in [('T1-ONT-tagged','T1'),('T2-ONT-tagged','T2')]:
            source=args.alignments/sid/'regions.bam'
            receipt=json.loads((source.parent/'receipt.json').read_text())
            selected=None
            by_segment=collections.defaultdict(list)
            with pysam.AlignmentFile(source) as bam:
                for r in bam:
                    by_segment[segment_key(r)].append(r)
            for flank in (60,45,30,20):
                groups=collections.defaultdict(list)
                for records in by_segment.values():
                    for r in records:
                        row=extract(r,c1,p1,c2,p2,flank,records)
                        if row:
                            key=json.dumps({k:row[k] for k in ('sequence','junction_start','junction_end','blocks')},sort_keys=True)
                            if row['read_id'] not in {x['read_id'] for x in groups[key]}:
                                groups[key].append(row)
                ranking=sorted(groups.values(), key=lambda rows:(-len({r['fragment_id'] for r in rows}), rows[0]['sequence']))
                counts=[len({r['fragment_id'] for r in rows}) for rows in ranking]
                print(name,sample,flank,counts[:15],flush=True)
                if ranking and counts[0] >= 2:
                    selected=ranking[0]
                    break
            if not selected:
                summaries.append(dict(event=name,sample=sample,status='no_repeated_exact_window_Q10'))
                continue
            row=selected[0]
            entry=dict(fusion=dict(event_id=name,reference_name='GRCh38',
                sequence=row['sequence'],junction_start=row['junction_start'],junction_end=row['junction_end'],
                donor=dict(contig=c1,position=p1,strand='+'),acceptor=dict(contig=c2,position=p2-1,strand='+'),
                blocks=row['blocks'],provenance=dict(sample_id=sample, method='direct identical original RNA windows',
                version='2',parameters=dict(min_mapq=20,min_base_quality=10,flank=flank,observed_partner_CIGARs=True),
                source=receipt['url'],contig_id=name+'-'+sample+'-local-window',
                acquisition=receipt,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())),
                references=refs,reads=[dict(sample_id=sample,library_id=sample+'-ONT',fragment_id=r['fragment_id'],
                  read_id=r['read_id'],source=receipt['url'],source_query_start=r['source_query_start'],
                  cdna_start=0,sequence=r['sequence'],blocks=r['blocks']) for r in selected],
                original_records=[dict(sam=r['original_sam'],partner_sam=r['partner_sam'],
                    source_reverse_complement=r['source_reverse_complement'], minimum_base_quality=r['minimum_base_quality']) for r in selected],
                note='Partial RNA window; selected repeated sequence group, not all event support or an assembled full transcript.')
            entry['reference_names'] = {r['transcript_id']:g.transcript_by_id(r['transcript_id']).name for r in refs}
            path=args.output/(name+'-'+sample+'.input.json.gz')
            path.write_bytes(gzip.compress((json.dumps(entry,sort_keys=True)+'\n').encode(), mtime=0))
            summaries.append(dict(event=name,sample=sample,input=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                  fragments=len({r['fragment_id'] for r in selected}),
                                  flank=flank,sequence=row['sequence']))
    (args.output/'manifest.json').write_text(json.dumps(summaries,indent=2)+'\n')


if __name__ == '__main__':
    main()
