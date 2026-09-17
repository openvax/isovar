"""Nine additional indel/SV RNA footprints, including unresolved consequences."""

from collections import defaultdict
import gzip
from hashlib import sha256
import json
from pathlib import Path

from isovar.fusion import fusion_from_dict, reconstruct_fusion
from isovar.fusion_visualization import _canvas, save_fusion_figures
from isovar.visualization import BLUE, GRAY, ORANGE, _plot_imports, _side_note, save_variant_figures

CORPUS=Path(__file__).resolve().parents[1]/'tests/data/osteosarc/figure_comparisons/corpus'
LABELS={'T1-ONT-dedup':'T1 ONT', 'T2-ONT-dedup':'T2 ONT', 'T1-short':'T1 short RNA', 'T2-short':'T2 short RNA'}


def load():
    manifest=json.loads((CORPUS/'rna-footprints-manifest.json').read_text())
    raw=(CORPUS/manifest['file']).read_bytes()
    if sha256(raw).hexdigest()!=manifest['sha256']:
        raise ValueError('RNA footprint fixture digest mismatch')
    return json.loads(gzip.decompress(raw))


def indel_panel(entry):
    figure,ax=_canvas(entry['gene']+' | RNA allele support','Exact listed allele; processing choices retained explicitly',5.6)
    columns=[('ref','Reference'),('alt','Alternate'),('other','Other')]
    for x,(_,title) in enumerate(columns,1):
        ax.text(x,4.6,title,ha='center',fontsize=12,weight='bold',color=GRAY)
    for y,p in zip([4,3,2,1],entry['products']):
        ax.text(0,y,LABELS[p['source']],ha='right',va='center',fontsize=12)
        for x,(key,_) in enumerate(columns,1):
            ax.text(x,y,str(p['default_counts'][key]),ha='center',va='center',fontsize=21,color=BLUE if key=='alt' else GRAY)
        ax.axhline(y-.4,color='#eeeeee',lw=.7)
    ax.set(xlim=(.1,3.5),ylim=(.3,5))
    ax.axis('off')
    _side_note(ax,'Template IDs\nIsovar default filters\n\nONT: deduplicated\nNo cross-product pooling\n\nNot independent\nmolecule counts')
    size=len(entry['alt'])-len(entry['ref'])
    figure.text(.19,.11,'%s:%s  %s>%s  (%+d nt)' % (entry['contig'],format(entry['start'],','),entry['ref'],entry['alt'],size),fontsize=12)
    figure.text(.19,.055,'Protein panels use primary alignments only; alternate-placement sensitivity is recorded in evidence.json.',fontsize=10,color=GRAY)
    return figure


def deletion_panels(entry):
    start,end=entry['start'],entry['end']
    models=[m for m in entry['models'] if m['biotype']=='protein_coding' and m['complete']]
    models.sort(key=lambda m:(-m['deleted_coding_bases'],m['id']))
    models=models[:5]
    junctions=defaultdict(int)
    for product in entry['products']:
        for j in product['junctions']:
            if j['operation']=='N' and j['annotated'] and j['start']<=start and j['end']>=end:
                junctions[j['start'],j['end']]+=j['templates']
    common=max(junctions,key=junctions.get) if junctions else None
    lo,hi=(min(start-1000,common[0]-500),max(end+1000,common[1]+500)) if common else (start-5000,end+5000)
    unit=1000
    figure,ax=_canvas(entry['gene']+' | DNA interval in transcript context',
                     '%s:%s-%s  |  %s-bp listed deletion' % (entry['contig'],format(start,','),format(end,','),format(end-start,',')),6)
    _,_,rectangle=_plot_imports()
    ax.axvspan(0,(end-start)/unit,color=ORANGE,alpha=.14)
    for y,m in enumerate(models[::-1],1):
        exons=m['exons']
        a,b=max(lo,exons[0][0]),min(hi,exons[-1][1])
        if a<b:
            ax.plot([(a-start)/unit,(b-start)/unit],[y,y],color=GRAY,lw=.9)
        for a,b in exons:
            a,b=max(lo,a),min(hi,b)
            if a<b:
                ax.add_patch(rectangle(((a-start)/unit,y-.12),(b-a)/unit,.24,color=BLUE))
        ax.text(-.035,y,m['id']+'\n('+m['name']+')',ha='right',va='center',transform=ax.get_yaxis_transform(),fontsize=10)
    ax.set(xlim=((lo-start)/unit,(hi-start)/unit),ylim=(.4,len(models)+.7),xlabel='Genomic offset from deletion start (kb)')
    _side_note(ax,'Orange: DNA interval\nBlue: annotated exon\n\nEnsembl 87 / GRCh38\nModels run on - strand\n\nNot an RNA assembly')
    note=('Deletes the DLG5-001 coding-start region; no mutant RNA junction reconstructed.' if entry['gene']=='DLG5' else
          'Intronic in coding models; normal splicing can remove this interval with or without the DNA deletion.')
    figure.text(.19,.065,note,fontsize=11,color=GRAY)
    yield 'transcript-interval',figure

    figure,ax=_canvas(entry['gene']+' | Observed RNA footprint',
                     'Ordinary splice paths are not evidence of the DNA deletion',5.8)
    columns=['Inside interval','Normal skip','Exact-size D / N']
    for x,title in enumerate(columns,1):
        ax.text(x,4.65,title,ha='center',fontsize=11,weight='bold',color=GRAY)
    for y,p in zip([4,3,2,1],entry['products']):
        ax.text(0,y,LABELS[p['source']],ha='right',va='center',fontsize=12)
        skip=next((j['templates'] for j in p['junctions'] if common and (j['start'],j['end'])==common and j['operation']=='N'),0)
        values=[str(p['aligned_templates']['inside']),str(skip) if common else 'n/a',
                '%d / %d' % (p['breakpoint_matching_templates']['deletion'],p['breakpoint_matching_templates']['splice'])]
        for x,value in enumerate(values,1):
            ax.text(x,y,value,ha='center',va='center',fontsize=20,color=BLUE)
        ax.axhline(y-.4,color='#eeeeee',lw=.7)
    ax.set(xlim=(.1,3.5),ylim=(.3,5))
    ax.axis('off')
    _side_note(ax,'Template IDs\nPrimary / MAPQ >=20\n\nD: CIGAR deletion\nN: CIGAR skipped region\nEndpoints within 3 bp\n\nZero is not proof\nof biological absence')
    footer=('Normal skip: %s-%s; this observed annotated intron contains the DNA interval.' % tuple(format(x,',') for x in common)
            if common else 'No normal splice spanning the entire interval; RNA within it does not exclude a subclonal deletion.')
    figure.text(.19,.09,footer,fontsize=10,color=GRAY)
    figure.text(.19,.045,'No coverage-normalized expression or allele-specific causality is inferred. Other junctions remain in evidence.json.',fontsize=10,color=GRAY)
    yield 'rna-footprint',figure


def fusion_inputs(entry,sources):
    c1,p1,c2,p2,s1,s2=entry['breakpoints']
    for product in entry['products']:
        if not product['paths']:
            continue
        for flank in (60,45,30,20):
            groups=defaultdict(list)
            for path in product['paths']:
                for w in path['windows']:
                    if w['flank']==flank:
                        key=json.dumps({k:w[k] for k in ('sequence','junction_start','junction_end','blocks')},sort_keys=True)
                        groups[key].append(w)
            ranking=sorted(groups.values(),key=lambda rows:(-len({r['fragment_id'] for r in rows}),rows[0]['sequence']))
            if ranking and len({w['fragment_id'] for w in ranking[0]})>=2:
                break
        if not ranking or len({w['fragment_id'] for w in ranking[0]})<2:
            continue
        windows=ranking[0]
        w=windows[0]
        source=sources[product['source']]
        sample=source['sample']
        yield dict(fusion=dict(event_id=entry['name'],reference_name='GRCh38',sequence=w['sequence'],
            junction_start=w['junction_start'],junction_end=w['junction_end'],blocks=w['blocks'],
            donor=dict(contig=c1,position=p1-1 if s1=='-' else p1,strand=s1),
            acceptor=dict(contig=c2,position=p2-1 if s2=='+' else p2,strand=s2),
            provenance=dict(sample_id=sample,method='direct identical original RNA windows',version='2',
                parameters=dict(min_mapq=20,min_base_quality=10,flank=flank,observed_partner_CIGARs=True),
                source=source['url'],contig_id=entry['name']+'-'+sample)),references=entry['references'],
            reads=[dict(sample_id=sample,library_id=source['id'],fragment_id=r['fragment_id'],read_id=r['read_id'],
                        source=source['url'],source_query_start=r['source_query_start'],cdna_start=0,sequence=r['sequence'],blocks=r['blocks']) for r in windows],
            original_records=[dict(sam=r['original_sam'],partner_sam=r['partner_sam'],source_reverse_complement=r['source_reverse_complement']) for r in windows])


def fusion_support_panel(entry):
    figure,ax=_canvas(entry['name'].replace('--',' / ')+' | RNA paths',
                     'Both actual partner alignments required; SA alone does not count',5.6)
    products={p['source']:p for p in entry['products']}
    rows=[('T1 ONT tagged','T1-ONT-tagged'),('T1 short RNA','T1-short'),
          ('T2 ONT tagged','T2-ONT-tagged'),('T2 short RNA','T2-short')]
    for y,(label,sid) in zip([4,3,2,1],rows):
        p=products[sid]
        ax.barh(y,p['complete_paths'],height=.35,color=BLUE)
        ax.text(-.03,y,label,ha='right',va='center',transform=ax.get_yaxis_transform(),fontsize=12)
        ax.text(p['complete_paths']+.15,y,str(p['complete_paths']),va='center',fontsize=18,color=BLUE)
    ax.set(xlim=(0,max(4,max(p['complete_paths'] for p in products.values())+2)),ylim=(.3,4.7),xlabel='Observed high-MAPQ chimeric read paths')
    _side_note(ax,'Same-segment pieces\nNot mate-pair links\nNot alternate placements\n\nONT tagged predecessors\nNot added to dedup counts\n\nFrame remains unresolved')
    labels=', '.join('%s: %d CB/UMI labels' % (s[:2],products[s]['cell_umi_labels']) for s in ('T1-ONT-tagged','T2-ONT-tagged'))
    figure.text(.19,.09,labels+'. Labels do not establish molecular independence.',fontsize=10,color=GRAY)
    figure.text(.19,.045,'Junction panels select repeated exact Q10 windows; their support is a subset of these event paths.',fontsize=10,color=GRAY)
    return figure


def save_panels(directory,panels):
    from matplotlib.backends.backend_pdf import PdfPages
    directory.mkdir(parents=True,exist_ok=True)
    with PdfPages(directory/'all-figures.pdf') as pdf:
        for name,figure in panels:
            figure.savefig(directory/(name+'.png'),dpi=600,facecolor='white')
            figure.savefig(directory/(name+'.svg'),facecolor='white',metadata={'Date':None})
            pdf.savefig(figure,facecolor='white')
            figure.clear()


def generate(output):
    data=load()
    root=output/'rna-footprints'
    root.mkdir()
    (root/'evidence.json').write_text(json.dumps(data,indent=2)+'\n')
    for entry in data['indels']:
        directory=root/entry['gene']
        save_panels(directory,[('allele-support',indel_panel(entry))])
        for product in entry['products']:
            if any(m['protein'] for m in product['visualization']['modes']):
                save_variant_figures(product['visualization'],directory/product['source'])
    for entry in data['deletions']:
        save_panels(root/entry['gene'],deletion_panels(entry))
    sources={s['id']:s for s in data['sources']}
    for entry in data['fusions']:
        save_panels(root/entry['name'],[('path-support',fusion_support_panel(entry))])
        for supplied in fusion_inputs(entry,sources):
            fusion,refs,reads=fusion_from_dict(supplied)
            result=reconstruct_fusion(fusion,refs,reads)
            directory=save_fusion_figures(result,root/entry['name'],refs,{m['id']:m['name'] for m in entry['models']})
            (directory/'input.json').write_text(json.dumps(supplied,indent=2)+'\n')
    return root
