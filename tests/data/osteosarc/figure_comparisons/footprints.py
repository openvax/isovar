"""Reproducible, bounded RNA footprint audit of nine osteosarc candidates."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import gzip
import json
import logging
from pathlib import Path
import subprocess
from collections import Counter, defaultdict

import pysam
from pyensembl import EnsemblRelease
from varcode import Variant
from isovar import ReadCollector
from isovar.read_identity import fragment_ids
from isovar.visualization import collect_visualization_data
from isovar import ProteinSequenceCreator
from tests.data.fusions.build_osteosarc import extract, linked_records, references, segment_key
from tests.data.osteosarc.expansion.references import apply_variant
from tests.data.osteosarc.expansion.runner import protein_check

INDELS = [
    dict(gene="GTF3C5", contig="chr9", start=133057893, ref="GGAGGAGGAGGAA", alt="G"),
    dict(gene="RNF213", contig="chr17", start=80327830, ref="ATAC", alt="A"),
    dict(gene="GLIS3", contig="chr9", start=3856149, ref="CTGATGTGG", alt="C"),
    dict(gene="KTN1", contig="chr14", start=55627965, ref="G", alt="GTT"),
]
DELETIONS = [
    dict(gene="DLG5", contig="chr10", start=77850921, end=77930452),
    dict(gene="AFF3", contig="chr2", start=99668095, end=99668223),
    dict(gene="KEAP1", contig="chr19", start=10492683, end=10492894),
]
FUSIONS = [
    ("GABBR1--SLC29A1", "chr6", 29612906, "chr6", 44218908, "-", "+"),
    ("OTUD7A--FMN1", "chr15", 31743670, "chr15", 33043479, "-", "+"),
]
LEGACY = [("chr7", 66205522), ("chr7", 66127704), ("chr6", 108611245),
          ("chr17", 63745980), ("chr2", 205310103), ("chr9", 22007648)]
ROOT = "https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/"


def source_products():
    for sample in ("T1", "T2"):
        stem = "IPISRC044_%s_sclrs_ONT" % sample
        base = ROOT + "ONT/IPISRC044_ONT_upload/IPISRC044_ONT/processed/" + stem + "/" + stem + "/"
        for product, suffix in (("ONT-dedup", stem + "_dedup/" + stem + "_dedup.bam"),
                                ("ONT-tagged", stem + ".tagged.bam")):
            yield dict(id=sample + "-" + product, sample=sample, product=product, url=base + suffix)
        yield dict(id=sample + "-short", sample=sample, product="short", url=ROOT +
                   "kamil/oncoanalyser/IPISRC044_%s_ucla/alignments/rna/IPISRC044_tumor_%s_ucla_rna.md.bam" % (sample, sample))


def regions():
    spans = [(v["contig"], v["start"] - 2000, v["start"] + 2000) for v in INDELS]
    spans += [(v["contig"], v["start"] - 5000, v["end"] + 5000) for v in DELETIONS]
    for _, c1, p1, c2, p2, _, _ in FUSIONS:
        spans += [(c1, p1 - 2000, p1 + 2000), (c2, p2 - 2000, p2 + 2000)]
    spans += [(c, p - 2000, p + 2000) for c, p in LEGACY]
    return ["%s:%d-%d" % row for row in spans]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def acquire(source, output):
    directory = output / source["id"]
    directory.mkdir(parents=True)
    bam = directory / "regions.bam"
    command = ["samtools", "view", "--no-PG", "-b", "-M", "-X", source["url"],
               source["url"] + ".bai", *regions(), "-o", str(bam)]
    subprocess.run(command, check=True, timeout=600, cwd=directory)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    receipt = dict(source, command=command, regions=regions(), assembly="GRCh38",
                   acquired_at=datetime.now(timezone.utc).isoformat(), bam_sha256=digest(bam),
                   index_sha256=digest(Path(str(bam) + ".bai")))
    (directory / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(source["id"], bam.stat().st_size, flush=True)
    return receipt


def transcript_models(genome, genes):
    models = []
    for gene in genes:
        for tid in genome.transcript_ids_of_gene_name(gene):
            t = genome.transcript_by_id(tid)
            coding_ranges = list(genome.db.query(
                select_column_names=['start', 'end'], filter_column='transcript_id',
                filter_value=tid, feature='CDS'))
            models.append(dict(id=tid, name=t.name, gene=gene, contig='chr'+t.contig,
                               strand=t.strand, biotype=t.biotype, complete=t.complete,
                               exons=[(a-1,b) for a,b in sorted(t.exon_intervals)],
                               start_codon_positions=list(t.start_codon_positions) if t.contains_start_codon else [],
                               coding_ranges=[(a-1,b) for a,b in coding_ranges]))
    return models


def deletion_footprint(records, deletion, models):
    """Keep RNA D operations, splice N operations, and aligned bases distinct."""
    start, end = deletion['start'], deletion['end']
    known = {(a[1],b[0]) for t in models for a,b in zip(t['exons'],t['exons'][1:])}
    junctions = defaultdict(set)
    exact = defaultdict(set)
    covered = defaultdict(set)
    for r in records:
        if r.flag & (4|256|512|1024|2048) or r.mapping_quality < 20:
            continue
        key = segment_key(r)[:2]
        reference = r.reference_start
        for op,n in r.cigartuples:
            if op in (0,7,8):
                for label,lo,hi in [('left',start-100,start),('inside',start,end),('right',end,end+100)]:
                    if reference < hi and reference+n > lo:
                        covered[label].add(key)
            elif op in (2,3):
                if reference < end+5000 and reference+n > start-5000:
                    junctions[(reference,reference+n,'D' if op==2 else 'N')].add(key)
                if abs(reference-start) <= 3 and abs(reference+n-end) <= 3:
                    exact['deletion' if op==2 else 'splice'].add(key)
            if op in (0,2,3,7,8):
                reference += n
    return dict(aligned_templates={k:len(covered[k]) for k in ('left','inside','right')},
                breakpoint_matching_templates={k:len(exact[k]) for k in ('deletion','splice')},
                junctions=[dict(start=a,end=b,operation=op,templates=len(keys),annotated=(a,b) in known)
                           for (a,b,op),keys in sorted(junctions.items(),key=lambda kv:(-len(kv[1]),kv[0]))
                           if op == 'N' or b-a >= 20])


def fusion_footprint(records, event):
    name,c1,p1,c2,p2,s1,s2 = event
    by_segment = defaultdict(list)
    for r in records:
        by_segment[segment_key(r)].append(r)
    paths=[]
    for group in by_segment.values():
        usable=[r for r in group if not r.flag & (4|256|512|1024) and r.mapping_quality>=20]
        sides=[[r for r in usable if r.reference_name==c and r.reference_start <= p+2000
                and r.reference_end >= p-2000] for c,p in ((c1,p1),(c2,p2))]
        if len(sides[0])!=1 or len(sides[1])!=1 or not linked_records(sides[0][0],sides[1][0],group):
            continue
        first=sides[0][0]
        row=dict(read_id=first.query_name, cell_umi=(first.get_tag('CB')+'/'+first.get_tag('UB')
                 if first.has_tag('CB') and first.has_tag('UB') else None),
                 original_records=[r.to_string() for r in group], windows=[])
        for flank in (60,45,30,20):
            for r in group:
                window=extract(r,c1,p1,c2,p2,flank,group,(s1,s2))
                if window:
                    row['windows'].append(dict(flank=flank,**window))
        paths.append(row)
    return dict(paths=paths, complete_paths=len(paths),
                cell_umi_labels=len({r['cell_umi'] for r in paths if r['cell_umi']}))


def split_deletion_paths(records,start,end):
    """Explicit supplementary adjacency at the listed genomic deletion ends."""
    groups=defaultdict(list)
    for r in records:
        if not r.flag & (4|256|512|1024) and r.mapping_quality>=20:
            groups[segment_key(r)].append(r)
    supported=[]
    for group in groups.values():
        left=[r for r in group if abs(r.reference_end-start)<=3]
        right=[r for r in group if abs(r.reference_start-end)<=3]
        if any(a.is_reverse==b.is_reverse and linked_records(a,b,group) for a in left for b in right):
            supported.append([r.to_string() for r in group])
    return supported


def audit(inputs, output):
    """Analyze original BAMs; retain unresolved outcomes rather than invent alleles."""
    output.mkdir(parents=True,exist_ok=False)
    genome=EnsemblRelease(87)
    from isovar import __version__
    from importlib.metadata import version
    data=dict(annotation='Ensembl 87 / GRCh38', indels=[], deletions=[], fusions=[], sources=[],
              software=dict(isovar=__version__,varcode=version('varcode'),pysam=pysam.__version__))
    sources=list(source_products())
    for source in sources:
        source['receipt']=json.loads((inputs/source['id']/'receipt.json').read_text())
    data['sources']=sources
    for row in INDELS:
        entry=dict(row,models=transcript_models(genome,[row['gene']]),products=[])
        variant=Variant(row['contig'].removeprefix('chr'),row['start'],row['ref'],row['alt'],ensembl=genome)
        expectations={}
        for tid in genome.transcript_ids_at_locus(variant.contig,variant.start,variant.end):
            t=genome.transcript_by_id(tid)
            if not t.complete or t.sequence is None:
                continue
            model=dict(contig=t.contig,exons=sorted(t.exon_intervals),strand=t.strand,cdna=t.sequence,
                       cds_start=min(t.start_codon_spliced_offsets),cds_end=max(t.stop_codon_spliced_offsets)+1,
                       genetic_code=1,transcript_version=tid,protein_version=t.protein_id)
            try:
                expectations[tid]=apply_variant(dict(chrom=row['contig'],pos=row['start'],ref=row['ref'],alt=row['alt']),model)
            except ValueError:
                continue  # Noncoding overlap remains in the transcript/effect tracks.
        for source in sources:
            if source['product']=='ONT-tagged':
                continue
            with pysam.AlignmentFile(inputs/source['id']/'regions.bam') as bam:
                evidence=ReadCollector(use_secondary_alignments=False).read_evidence_for_variant(variant,bam)
                default_evidence=ReadCollector().read_evidence_for_variant(variant,bam)
            counts={k:len(fragment_ids(getattr(evidence,k+'_reads'))) for k in ('ref','alt','other')}
            default_counts={k:len(fragment_ids(getattr(default_evidence,k+'_reads'))) for k in ('ref','alt','other')}
            visualization=collect_visualization_data(variant,evidence,compare_assembly=True)
            validation=[]
            for assembly in (True,False):
                creator=ProteinSequenceCreator(variant_sequence_assembly=assembly)
                proteins=creator.sorted_protein_sequences_for_variant(variant,evidence)
                checked=protein_check(proteins[0],expectations,creator.protein_sequence_length,
                                      creator.protein_context_peptide_length) if proteins else None
                if checked and checked['validation_status']!='ok':
                    raise ValueError('Independent protein validation failed: %s %s %s' % (row['gene'],source['id'],checked))
                validation.append(dict(assembly=assembly,result=checked))
            visualization.setdefault('provenance',{}).update(source=source, sample_label=source['id'],
                count_unit='RG/QNAME templates; not proven independent molecules',secondary_alignments=False,
                independent_validation=validation)
            entry['products'].append(dict(source=source['id'],counts=counts,default_counts=default_counts,visualization=visualization))
            print(row['gene'],source['id'],counts,[len(m['protein']['amino_acids']) if m['protein'] else 0
                                                 for m in visualization['modes']],flush=True)
        data['indels'].append(entry)
    for row in DELETIONS:
        models=transcript_models(genome,[row['gene']])
        entry=dict(row,models=models,products=[],supplementary_products=[])
        for t in models:
            t['deleted_exonic_bases']=sum(max(0,min(b,row['end'])-max(a,row['start'])) for a,b in t['exons'])
            t['deleted_coding_bases']=sum(max(0,min(b,row['end'])-max(a,row['start'])) for a,b in t['coding_ranges'])
        for source in sources:
            with pysam.AlignmentFile(inputs/source['id']/'regions.bam') as bam:
                records=list(bam.fetch(row['contig'],row['start']-5000,row['end']+5000))
            paths=split_deletion_paths(records,row['start'],row['end'])
            entry['supplementary_products'].append(dict(source=source['id'],complete_paths=len(paths),original_records=paths))
            if source['product']=='ONT-tagged':
                continue
            footprint=deletion_footprint(records,row,models)
            entry['products'].append(dict(source=source['id'],**footprint))
            print(row['gene'],source['id'],footprint['aligned_templates'],footprint['breakpoint_matching_templates'],flush=True)
        data['deletions'].append(entry)
    for event in FUSIONS:
        name,c1,p1,c2,p2,s1,s2=event
        genes=name.split('--')
        entry=dict(name=name,breakpoints=event[1:],models=transcript_models(genome,genes),
                   references=references(genome,genes),products=[])
        for source in sources:
            with pysam.AlignmentFile(inputs/source['id']/'regions.bam') as bam:
                records=list(bam.fetch(c1,p1-2000,p1+2000))+list(bam.fetch(c2,p2-2000,p2+2000))
            footprint=fusion_footprint(records,event)
            entry['products'].append(dict(source=source['id'],**footprint))
            print(name,source['id'],footprint['complete_paths'],Counter(w['flank'] for p in footprint['paths'] for w in p['windows']),flush=True)
        data['fusions'].append(entry)
    (output/'evidence.json').write_text(json.dumps(data,indent=2)+'\n')
    packed=output/'rna-footprints.json.gz'
    packed.write_bytes(gzip.compress((json.dumps(data,sort_keys=True)+'\n').encode(),mtime=0))
    (output/'rna-footprints-manifest.json').write_text(json.dumps(dict(file=packed.name,sha256=digest(packed),
        inputs=str(inputs),annotation=data['annotation'],candidates=9,
        scope='T1/T2, ONT and oncoanalyser short RNA; no claim about unqueried libraries',
        records='Complete regional BAMs retained separately; fusion path SAMs embedded; indel/count summaries require reacquisition'),indent=2)+'\n')
    return data


def main():
    logging.disable(logging.INFO)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path, help="New destination; never overwritten")
    parser.add_argument("--inputs", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.acquire:
        args.acquire.mkdir(parents=True, exist_ok=False)
        with ThreadPoolExecutor(max_workers=3) as pool:
            receipts = list(pool.map(lambda source: acquire(source, args.acquire), source_products()))
        (args.acquire / "manifest.json").write_text(json.dumps(receipts, indent=2) + "\n")
    else:
        audit(args.inputs,args.output)


if __name__ == "__main__":
    main()
