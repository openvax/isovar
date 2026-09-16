"""Extract observed, directly spanning windows; never synthesize RNA bases."""
import argparse
import collections
import gzip
import hashlib
import json
from pathlib import Path

import pysam
from pyensembl import EnsemblRelease

EVENTS = [('TPST1--CRCP', 'chr7', 66205522, 'chr7', 66127704, ['TPST1'], ['CRCP']),
          ('FOXO3--STRADA-CCDC47', 'chr6', 108611245, 'chr17', 63745980, ['FOXO3'], ['STRADA', 'CCDC47']),
          ('PARD3B--CDKN2B-AS1-CDKN2B', 'chr2', 205310103, 'chr9', 22007648, ['PARD3B'], ['CDKN2B-AS1', 'CDKN2B'])]


def blocks(mapping, lo, hi):
    result = []
    for q in range(lo, hi):
        if q not in mapping:
            continue
        chrom, pos = mapping[q]
        if result and result[-1]['contig'] == chrom and result[-1]['query_end'] == q - lo and result[-1]['reference_end'] == pos:
            result[-1]['query_end'] += 1
            result[-1]['reference_end'] += 1
        else:
            result.append(dict(query_start=q-lo, query_end=q-lo+1, contig=chrom,
                               reference_start=pos, reference_end=pos+1, strand='+'))
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


def extract(read, c1, p1, c2, p2, flank):
    if read.flag & (4 | 256 | 512 | 1024 | 2048) or not read.has_tag('SA') or read.mapping_quality < 20:
        return None
    sequence = read.query_sequence
    if sequence is None or any(op == 5 for op, n in read.cigartuples):
        return None
    alignments = [(read.reference_name, read.reference_start, '-' if read.is_reverse else '+', read.cigarstring, read.mapping_quality)]
    for sa in read.get_tag('SA').strip(';').split(';'):
        chrom, position, strand, cigar, mapq, nm = sa.split(',')
        alignments.append((chrom, int(position)-1, strand, cigar, int(mapq)))
    candidates = []
    for chrom, position, strand, cigar, mapq in alignments:
        if mapq < 20 or chrom not in (c1, c2):
            continue
        # These three explicit events are +/+. In BAM SEQ is already oriented
        # with the primary reference, irrespective of the sequencer direction.
        if strand != ('-' if read.is_reverse else '+'):
            continue
        alignment = pysam.AlignedSegment()
        alignment.reference_start = position
        alignment.cigarstring = cigar.replace('H', 'S')
        if alignment.infer_query_length() != len(sequence):
            continue
        mapping = {q:(chrom, r) for q,r in alignment.get_aligned_pairs(matches_only=True)}
        candidates.append(mapping)
    donors = [(m, q) for m in candidates for q,p in m.items() if p == (c1, p1-1)]
    acceptors = [(m, q) for m in candidates for q,p in m.items() if p == (c2, p2-1)]
    if len(donors) != 1 or len(acceptors) != 1:
        return None
    donor, q0 = donors[0]
    acceptor, q1 = acceptors[0]
    q0 += 1
    if not q0 <= q1 or q1-q0 > 100:
        return None
    lo, hi = q0-flank, q1+flank
    if lo < 0 or hi > len(sequence):
        return None
    quality = read.query_qualities
    if quality is None or min(quality[lo:hi]) < 10:
        return None
    mapping = {q:p for q,p in donor.items() if lo <= q < q0}
    mapping.update({q:p for q,p in acceptor.items() if q1 <= q < hi})
    local = sequence[lo:hi]
    if set(local)-set('ACGT'):
        return None
    return dict(sequence=local, junction_start=flank, junction_end=q1-lo,
                blocks=blocks(mapping,lo,hi), source_query_start=lo,
                minimum_base_quality=min(quality[lo:hi]), read_id=read.query_name,
                fragment_id=(read.get_tag('CB') + '/' + read.get_tag('UB')
                             if read.has_tag('CB') and read.has_tag('UB') else read.query_name),
                original_sam=read.to_string())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--alignments', type=Path, required=True,
                        help='Acquired regional BAM directory, with source-id/regions.bam and receipt.json')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    g=EnsemblRelease(87)
    summaries=[]
    for name,c1,p1,c2,p2,left_genes,right_genes in EVENTS:
        refs=references(g,left_genes+right_genes)
        for sid,sample in [('7010dd389c50467d','T1'),('2bc1fc291308debb','T2')]:
            source=args.alignments/sid/'regions.bam'
            receipt=json.loads((source.parent/'receipt.json').read_text())
            selected=None
            for flank in (60,45,30,20):
                groups=collections.defaultdict(list)
                with pysam.AlignmentFile(source) as bam:
                    for r in bam:
                        row=extract(r,c1,p1,c2,p2,flank)
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
                version='1',parameters=dict(min_mapq=20,min_base_quality=10,flank=flank,primary_and_SA_only=True),
                source=receipt['url'],contig_id=name+'-'+sample+'-local-window',
                acquisition=receipt,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())),
                references=refs,reads=[dict(sample_id=sample,library_id=sample+'-ONT',fragment_id=r['fragment_id'],
                  read_id=r['read_id'],source=receipt['url'],source_query_start=r['source_query_start'],
                  cdna_start=0,sequence=r['sequence'],blocks=r['blocks']) for r in selected],
                original_records=[dict(sam=r['original_sam'],minimum_base_quality=r['minimum_base_quality']) for r in selected],
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
