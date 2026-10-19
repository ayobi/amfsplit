#!/usr/bin/env python3
"""Reproducible fixed-backbone LSU placement and uncertainty reports.

Commands:
  build --v18 DIR --queries FASTA --out DIR --mafft PATH --raxml-ng PATH --epa-ng PATH
  summarize --jplace FILE --metadata TSV --queries FASTA --aligned-query FASTA --out DIR

The reference has no culture copies. Taxonomic labels in summarize are external
V18 annotations, not held-out culture predictions. Benchmark interpretation is
provided separately by compare_holdouts.py.
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from Bio import Phylo, SeqIO

EXPECTED_AMF_SHA='e6794cbbe166735a94e44b146aae2ecffb35ebfee61cf6ae8ed2e48ba9227198'
REFERENCE_COMMIT='f4a0014336f49be6aea48e8374dae052360a88bd'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def read_fasta(p):
    out={}
    for r in SeqIO.parse(p,'fasta'):
        if r.id in out:raise ValueError(f'Duplicate FASTA ID: {r.id}')
        out[r.id]=str(r.seq).upper()
    if not out:raise ValueError(f'Empty FASTA: {p}')
    return out

def write_fasta(p, seqs):
    with Path(p).open('w') as f:
        for x,s in seqs.items():f.write(f'>{x}\n{s}\n')

def read_tsv(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(p,rows):
    rows=list(rows)
    if not rows:Path(p).write_text('');return
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

def write_json(p,data):Path(p).write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')

def credible_edges(placements,mass=.95):
    """Smallest highest-weight set reaching mass, expanded for boundary ties."""
    if not 0<mass<=1:raise ValueError('Credible mass must be in (0,1]')
    ps=sorted(placements,key=lambda r:(-r['like_weight_ratio'],r['edge_num']))
    total=math.fsum(r['like_weight_ratio'] for r in ps)
    if not math.isclose(total,1.,abs_tol=1e-6):raise ValueError(f'Placement weights do not sum to 1: {total}')
    chosen=[];acc=0.;cut=None
    for r in ps:
        w=r['like_weight_ratio']
        if cut is not None and not math.isclose(w,cut,rel_tol=1e-12,abs_tol=0.):break
        chosen.append(r);acc+=w
        if cut is None and acc>=mass:cut=w
    return chosen


def jplace_tree(tree_string,metadata):
    """Parse edge-tagged Newick; orient descendant sets away from the outgroup.

    Rooting uses a specified outgroup tip solely to orient bipartitions. Edge
    identifiers remain those emitted by EPA; no rerooted IDs are substituted.
    """
    replaced=re.sub(r'\{(\d+)\}',r'[&edge=\1]',tree_string)
    tree=Phylo.read(io.StringIO(replaced),'newick')
    refs=set(metadata);tips={n.name for n in tree.get_terminals()}
    if tips!=refs:raise ValueError('Jplace tree and metadata taxon sets differ')
    outgroups=sorted(x for x,r in metadata.items() if r['is_amf']=='0')
    if not outgroups:raise ValueError('Missing outgroup metadata')
    sentinel=outgroups[0]
    nodes=list(tree.find_clades(order='postorder'));desc={}
    for n in nodes:
        desc[n]={n.name} if n.is_terminal() else set().union(*(desc[c] for c in n.clades))
    edges={}
    for n in nodes:
        if not n.comment:continue
        match=re.fullmatch(r'&edge=(\d+)',n.comment)
        if not match:raise ValueError(f'Unrecognized edge comment: {n.comment}')
        e=int(match[1]);side=desc[n]
        if sentinel in side:side=refs-side
        amf=side and all(metadata[t]['is_amf']=='1' for t in side)
        spp={metadata[t]['label'] for t in side} if amf else set()
        genera={metadata[t]['genus'] for t in side} if amf else set()
        edges[e]={'edge':e,'descendant_tips':sorted(side),'is_amf':bool(amf),'species':sorted(spp),'genera':sorted(genera),
                  'branch_length':n.branch_length or 0.,'unit':f'E{e}','genus_unit':f'E{e}'}
    # Collapse only nested pure-species components, not every occurrence of a
    # repeated species label on the tree. This preserves polyphyletic labels.
    for field,unit_field,prefix in [('species','unit','SC'),('genera','genus_unit','GC')]:
        pure=defaultdict(list)
        for e,d in edges.items():
            if len(d[field])==1:pure[d[field][0]].append((e,set(d['descendant_tips'])))
        for label,components in pure.items():
            for e,ds in components:
                supers=[(p,ps) for p,ps in components if ds<=ps]
                p,_=max(supers,key=lambda x:(len(x[1]),-x[0]))
                edges[e][unit_field]=f'{prefix}{p}'
    return edges


def summarize(jplace,metadata_path,query_path,aligned_path,out,mass=.95):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    meta={r['ref_id']:r for r in read_tsv(metadata_path)}
    data=json.loads(Path(jplace).read_text());edges=jplace_tree(data['tree'],meta)
    raw=read_fasta(query_path);aligned=read_fasta(aligned_path)
    if set(raw)!=set(aligned):raise ValueError('Aligned/raw query IDs differ')
    counts=Counter();rows=[];distributions={}
    for rec in data['placements']:
        ps=[dict(zip(data['fields'],p)) for p in rec['p']]
        for p in ps:
            p['edge_num']=int(p['edge_num'])
            if p['edge_num'] not in edges:raise ValueError('Unknown edge ID')
            if not math.isfinite(p['like_weight_ratio']) or not 0<=p['like_weight_ratio']<=1:raise ValueError('Invalid LWR')
        if len({p['edge_num'] for p in ps})!=len(ps):raise ValueError('Duplicate edge for a query')
        ids=rec.get('n') or [n for n,_ in rec['nm']]
        cs=credible_edges(ps,mass);best=max(ps,key=lambda p:p['like_weight_ratio'])
        species_mass=defaultdict(float);genus_mass=defaultdict(float);amf_mass=0.
        for p in ps:
            info=edges[p['edge_num']];w=p['like_weight_ratio']
            if info['is_amf']:amf_mass+=w
            if len(info['species'])==1:species_mass[info['unit']]+=w
            if len(info['genera'])==1:genus_mass[info['genus_unit']]+=w
        best_species_unit,maxsp=max(species_mass.items(),key=lambda x:x[1]) if species_mass else ('',0.)
        best_genus_unit,maxgn=max(genus_mass.items(),key=lambda x:x[1]) if genus_mass else ('',0.)
        best_species=next((d['species'][0] for d in edges.values() if d['unit']==best_species_unit and len(d['species'])==1),'')
        best_genus=next((d['genera'][0] for d in edges.values() if d['genus_unit']==best_genus_unit and len(d['genera'])==1),'')
        rank='species_clade' if maxsp>=mass else 'genus_clade' if maxgn>=mass else 'AMF_broad' if amf_mass>=mass else 'unresolved'
        spp=sorted({s for p in cs for s in edges[p['edge_num']]['species']});gns=sorted({g for p in cs for g in edges[p['edge_num']]['genera']})
        for x in ids:
            if x in distributions or x not in raw:raise ValueError(f'Duplicate/unknown placement query {x}')
            aligned_bases=sum(c not in '-.' for c in aligned[x]);raw_bases=sum(c not in '-.' for c in raw[x])
            retention=aligned_bases/raw_bases if raw_bases else 0.
            if aligned_bases>raw_bases:raise ValueError('Alignment adds query nucleotides')
            usable=retention>=.9 and aligned_bases>0
            rows.append({'accession':x,'aligned_bases':aligned_bases,'raw_bases':raw_bases,'retained_fraction':retention,'deleted_insertions':raw_bases-aligned_bases,'ambiguous_aligned_bases':sum(c not in 'ACGT-.' for c in aligned[x]),'bridge_alignment_eligible':int(usable),
                         'best_edge':best['edge_num'],'best_edge_lwr':best['like_weight_ratio'],'best_pendant_length':best['pendant_length'],'best_distal_length':best['distal_length'],
                         'lwr_sum':math.fsum(p['like_weight_ratio'] for p in ps),'lwr_entropy_nats':-math.fsum(p['like_weight_ratio']*math.log(p['like_weight_ratio']) for p in ps if p['like_weight_ratio']>0),
                         'credible_mass_target':mass,'credible_mass_retained':math.fsum(p['like_weight_ratio'] for p in cs),'credible_edge_count':len(cs),'credible_edges':','.join(str(p['edge_num']) for p in cs),
                         'credible_species_clade_units':','.join(sorted({edges[p['edge_num']]['unit'] for p in cs})),
                         'amf_edge_mass':amf_mass,'external_species_clade':best_species,'external_species_clade_unit':best_species_unit,'external_species_clade_mass':maxsp,'external_genus_clade':best_genus,'external_genus_clade_unit':best_genus_unit,'external_genus_clade_mass':maxgn,
                         'external_rank':rank,'external_species_in_credible_edges':';'.join(spp),'external_genera_in_credible_edges':';'.join(gns),
                         'interpretation':'conditional placement support; external labels are not validated culture predictions'})
            distributions[x]=[{k:p[k] for k in ('edge_num','like_weight_ratio','pendant_length','distal_length')} for p in ps]
            counts[rank]+=1
    if set(distributions)!=set(raw):raise ValueError('Not every query has a placement')
    rows.sort(key=lambda r:r['accession'])
    write_tsv(out/'lsu_placements.copies.tsv',rows)
    write_json(out/'edge_annotations.json',edges)
    import gzip
    with gzip.open(out/'placement_weights.json.gz','wt') as f:json.dump(distributions,f,sort_keys=True)
    summary={'queries':len(rows),'external_rank_counts':dict(counts),'alignment_eligible':sum(r['bridge_alignment_eligible'] for r in rows),'credible_mass':mass,'reference_tips':len(meta),'tree_edges':len(edges),'placements_sha256':sha(jplace),'metadata_sha256':sha(metadata_path),'raw_query_sha256':sha(query_path),'aligned_query_sha256':sha(aligned_path)}
    write_json(out/'placement_summary.json',summary)
    return summary


def run_logged(command,log,stdout=None):
    print('Running '+Path(command[0]).name,flush=True)
    with Path(log).open('w') as err:
        if stdout:
            with Path(stdout).open('w') as output:subprocess.run(command,stdout=output,stderr=err,check=True)
        else:subprocess.run(command,stdout=err,stderr=subprocess.STDOUT,check=True)


def build(a):
    out=a.out.resolve()
    if out.exists() and any(out.iterdir()):raise ValueError('Output must be new or empty; earlier results are preserved')
    out.mkdir(parents=True,exist_ok=True);ref=out/'reference';ref.mkdir();logs=out/'logs';logs.mkdir();epa=out/'epa';epa.mkdir()
    v18=a.v18.resolve()
    amf_path=v18/'V18_LSUDB_052025_AMFONLY.fasta';full_path=v18/'V18_LSUDB_052025.fasta';tree_path=v18/'Root_V18_LSUDB_052025.newick'
    if sha(amf_path)!=EXPECTED_AMF_SHA:raise ValueError('V18 version differs from frozen benchmark')
    full=read_fasta(full_path);amf=read_fasta(amf_path);tree=Phylo.read(tree_path,'newick')
    if {n.name for n in tree.get_terminals()}!=set(full):raise ValueError('Tree and reference sequences differ')
    sys.path.insert(0,str(Path(__file__).resolve().parent/'code/baseline'))
    from arms import parse_delavaux
    mapping={x:f'R{i+1:04d}' for i,x in enumerate(full)}
    metadata=[{'ref_id':mapping[x],'original_id':x,'is_amf':int(x in amf),**parse_delavaux(x)} for x in full]
    write_tsv(ref/'reference_metadata.tsv',metadata);write_fasta(ref/'reference.unaligned.fasta',{mapping[x]:s for x,s in full.items()})
    for n in tree.find_clades():
        n.name=mapping[n.name] if n.is_terminal() else None;n.confidence=None;n.comment=None
    Phylo.write(tree,ref/'backbone.topology.newick','newick',format_branch_length='%1.15f')
    queries=out/'queries.fasta';shutil.copyfile(a.queries,queries)
    commands=[]
    def run(cmd,log,stdout=None):commands.append(cmd);run_logged(cmd,logs/log,stdout)
    run([a.mafft,'--localpair','--maxiterate','1000','--thread',str(a.threads),'--threadit','0',str(ref/'reference.unaligned.fasta')],'reference_mafft.log',ref/'reference.aligned.fasta')
    run([a.raxml_ng,'--evaluate','--msa',str(ref/'reference.aligned.fasta'),'--tree',str(ref/'backbone.topology.newick'),'--model','GTR+G4','--threads',str(min(a.threads,4)),'--seed','20260929','--prefix',str(ref/'fit'),'--force','perf_threads'],'reference_model.log')
    run([a.mafft,'--auto','--addfragments',str(queries),'--keeplength','--mapout','--thread',str(min(a.threads,4)),str(ref/'reference.aligned.fasta')],'query_mafft.log',out/'combined.aligned.fasta')
    combined=read_fasta(out/'combined.aligned.fasta');rseq=read_fasta(ref/'reference.aligned.fasta');qseq=read_fasta(queries)
    if set(combined)!=set(rseq)|set(qseq) or any(combined[x]!=s for x,s in rseq.items()):raise ValueError('Reference alignment changed or IDs differ')
    write_fasta(out/'queries.aligned.fasta',{x:combined[x] for x in qseq})
    run([a.epa_ng,'--ref-msa',str(ref/'reference.aligned.fasta'),'--tree',str(ref/'fit.raxml.bestTree'),'--query',str(out/'queries.aligned.fasta'),'--model',str(ref/'fit.raxml.bestModel'),'--outdir',str(epa),'--threads',str(a.threads),'--no-heur','--no-pre-mask','--filter-min-lwr','0','--filter-max','2000','--precision','12'],'epa.log')
    summary=summarize(epa/'epa_result.jplace',ref/'reference_metadata.tsv',queries,out/'queries.aligned.fasta',out/'summary',a.mass)
    write_json(out/'run.json',{'commands':commands,'reference_commit':REFERENCE_COMMIT,'reference_sha256':{p.name:sha(p) for p in (amf_path,full_path,tree_path)},'script_sha256':sha(__file__),'summary':summary})


def main():
    ap=argparse.ArgumentParser(description=__doc__);sub=ap.add_subparsers(dest='command',required=True)
    b=sub.add_parser('build');b.add_argument('--v18',type=Path,required=True);b.add_argument('--queries',type=Path,required=True);b.add_argument('--out',type=Path,required=True)
    b.add_argument('--mafft',default='mafft');b.add_argument('--raxml-ng',default='raxml-ng');b.add_argument('--epa-ng',default='epa-ng');b.add_argument('--threads',type=int,default=8);b.add_argument('--mass',type=float,default=.95)
    s=sub.add_parser('summarize')
    for n in ['jplace','metadata','queries','aligned-query','out']:s.add_argument('--'+n,type=Path,required=True)
    s.add_argument('--mass',type=float,default=.95)
    a=ap.parse_args()
    if a.command=='build':build(a)
    else:print(json.dumps(summarize(a.jplace,a.metadata,a.queries,a.aligned_query,a.out,a.mass),indent=2))
if __name__=='__main__':main()
