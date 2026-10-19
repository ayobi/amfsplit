#!/usr/bin/env python3
"""Training-only comparison of EPA LSU units with the frozen v4 baseline.

Two prespecified development variants: exact credible edges and maximal pure
V18-species clades. External labels define geometry only; output names are
inferred from training culture holders. No threshold is fitted to test labels.
"""
import argparse
from collections import Counter, defaultdict
from dataclasses import replace
import gzip
import json
from pathlib import Path
import statistics
import sys
import time
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'code/baseline'))
from bridge import Data, ArmTrainer, Fit, reconcile
from decision_tree import evaluate
from lsu_placement import read_tsv,write_tsv,write_json,sha


def placement_fit(data,training,sets,arm='LSU'):
    training=frozenset(training);cs={};holders=defaultdict(set)
    for c in sorted(training):
        nonempty=[sets[x] for x in data.by_culture[c] if data.cop[x]['copy_class']=='main' and sets[x]]
        inter=frozenset.intersection(*nonempty) if nonempty else frozenset()
        cs[c]=inter or frozenset().union(*nonempty)
        for u in cs[c]:holders[u].add(c)
    return Fit(arm,training,{},0.,0.,0.,[],sets,cs,dict(holders),{},data.best['LSU'],[])


def metrics(rows):
    known=[r for r in rows if r['copy_class']=='main' and r['group']=='named, species elsewhere']
    held=[r for r in rows if r['copy_class']=='main' and r['group']=='named, species held out']
    short=[r for r in rows if r['div_kind']=='short_5.8S'];mains=[r for r in rows if r['copy_class']=='main']
    statuses=defaultdict(Counter);bysp=defaultdict(list);bycu=defaultdict(list)
    for r in rows:
        g=f"divergent ({r['div_kind']})" if r['copy_class']=='divergent' else r['group'];statuses[g][r['status']]+=1
    for r in held:
        bysp[r['organism']].append(r['status'].startswith('known'));bycu[r['culture']].append(r['status'].startswith('known'))
    return {'copies':len(rows),'known_n':len(known),'known_retained':sum(r['status'].startswith('known') for r in known),
            'known_compatible':sum(r['status'].startswith('known') and r['organism'] in r['placement'].split(';') for r in known),
            'known_wrong':sum(r['status'].startswith('known') and r['organism'] not in r['placement'].split(';') for r in known),
            'reconcile_exact':sum(r['rank']=='species' and r['correct']=='yes' for r in known),
            'reconcile_containing_group':sum(r['rank']=='species group' and r['correct']=='contains' for r in known),
            'held_n':len(held),'false_known':sum(r['status'].startswith('known') for r in held),
            'false_known_species_macro':statistics.mean(statistics.mean(x) for x in bysp.values()) if bysp else None,
            'false_known_culture_macro':statistics.mean(statistics.mean(x) for x in bycu.values()) if bycu else None,
            'short_n':len(short),'short_recognised':sum(r['status']=='divergent class' for r in short),
            'main_divergent':sum(r['status']=='divergent class' for r in mains),'artifacts':sum(r['status']=='artifact' for r in rows),
            'statuses':dict(statuses)}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--placements',type=Path,default=ROOT/'placement/summary/lsu_placements.copies.tsv');ap.add_argument('--inputs',type=Path,default=ROOT/'inputs/benchmark');ap.add_argument('--out',type=Path,default=ROOT/'comparison')
    a=ap.parse_args();out=a.out.resolve()
    if out.exists() and any(out.iterdir()):ap.error('Choose a new output directory')
    out.mkdir(parents=True,exist_ok=True);start=time.monotonic();data=Data(a.inputs)
    placements={r['accession']:r for r in read_tsv(a.placements)}
    if set(placements)!=set(data.cop):raise ValueError('Placement/culture reference IDs differ')
    psets={}
    for policy,field in [('epa_edges','credible_edges'),('epa_species_clades','credible_species_clade_units')]:
        psets[policy]={x:(frozenset(filter(None,r[field].split(','))) if r['bridge_alignment_eligible']=='1' else frozenset()) for x,r in placements.items()}
    trainers={n:ArmTrainer(data,n) for n in ['VT','SH','LSU']};summaries={};allchanges=[];foldstats=[]
    for mode in ['physical_culture','physical_species']:
        results={p:[] for p in ['baseline',*psets]};folds=data.folds(mode)
        for i,(drop,xs) in enumerate(folds,1):
            train=data.cultures-drop;fits={n:t.fit(train) for n,t in trainers.items()}
            for policy in results:
                model=dict(fits)
                if policy!='baseline':model['LSU']=placement_fit(data,train,psets[policy])
                if any(set(v)&drop for v in model['LSU'].holders.values()):raise AssertionError('Held culture entered LSU holders')
                rec={x:reconcile(data,x,model,drop) for x in xs}
                nov=evaluate(data,xs,model,drop,rec)
                for r in nov:
                    x=r['accession'];r={**rec[x],**r,'policy':policy,'mode':mode,'held_cultures':';'.join(sorted(drop))}
                    r.update({f'lsu_{k}':placements[x][k] for k in ['best_edge_lwr','credible_edge_count','best_pendant_length','retained_fraction','external_species_clade','external_species_clade_mass','external_genus_clade','external_genus_clade_mass']})
                    # LWR supported by the training-holder unit set; overlap is
                    # a compatibility measure, never a posterior on a name.
                    results[policy].append(r)
            foldstats.append({'mode':mode,'fold':i,'held_cultures':';'.join(sorted(drop)),'training_cultures':len(train),'queries':len(xs)})
            if i%20==0 or i==len(folds):print(f'{mode}: {i}/{len(folds)} folds, {time.monotonic()-start:.1f}s',flush=True)
        baseline={r['accession']:r for r in results['baseline']}
        expected={r['accession']:r for r in read_tsv(ROOT/'inputs/previous_results'/f'train_calibrated_{mode}.copies.tsv')}
        parity=[]
        for x,r in baseline.items():
            for k in ['status','placement','rank','call_species','call_genera','correct','reasons','s58']:
                if r[k]!=expected[x][k]:parity.append((x,k,r[k],expected[x][k]))
        if parity:raise AssertionError(f'Baseline mismatch: {parity[:5]}')
        for policy,rows in results.items():
            rows.sort(key=lambda r:r['accession']);assert len(rows)==903 and len({r['accession'] for r in rows})==903
            write_tsv(out/f'{policy}_{mode}.copies.tsv',rows);summaries[f'{policy}_{mode}']=metrics(rows)
            if policy=='baseline':continue
            for r in rows:
                b=baseline[r['accession']]
                if r['status']!=b['status'] or r['placement']!=b['placement']:
                    allchanges.append({'policy':policy,'mode':mode,'accession':r['accession'],'culture':r['culture'],'organism':r['organism'],'copy_class':r['copy_class'],'group':r['group'],'baseline_status':b['status'],'epa_status':r['status'],'baseline_placement':b['placement'],'epa_placement':r['placement'],'baseline_LSU_state':b['LSU_state'],'epa_LSU_state':r['LSU_state'],'reason':r['reasons']})
        write_json(out/'summary.json',summaries);write_tsv(out/'changes.tsv',allchanges)
    write_tsv(out/'folds.tsv',foldstats)
    write_json(out/'run.json',{'script_sha256':sha(__file__),'placements_sha256':sha(a.placements),'reference_fingerprint':data.manifest['fingerprint'],'seconds':time.monotonic()-start,'baseline_parity':True,'protocol':'PROTOCOL.md','protocol_sha256':sha(ROOT/'PROTOCOL.md'),'interpretation':'Prespecified development variants; external databases fixed; no independent validation claim.'})
    print(json.dumps(summaries,indent=2))
if __name__=='__main__':main()
