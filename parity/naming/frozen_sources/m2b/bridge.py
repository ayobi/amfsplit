"""Training-only arm merges and calibration from lossless frozen unit evidence.

No BLAST needed: external reference units and pairwise query/reference scores
remain fixed. This module does not remove taxa from the external databases.
"""
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import gzip
import hashlib
import json
import math
import statistics

from arms import read_fasta, read_tsv, GENUS_SYNONYMS, is_named
from arm_sets import ARMS, dist, union_find, vt_set

GRID = tuple(round(i * .1, 3) for i in range(31))
FIXED = {'VT': .6, 'SH': 0., 'LSU': 0.}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def split(s, sep=','):
    return frozenset(filter(None, s.split(sep)))


@dataclass
class Fit:
    arm: str
    training: frozenset
    unit_map: dict
    delta: float
    d_cons: float
    d_spread: float
    calibration: list
    sets: dict
    culture_sets: dict
    holders: dict
    anchors: dict
    best: dict
    proposals: list


class Data:
    def __init__(self, inputs):
        self.inputs = Path(inputs)
        ref = self.inputs / 'ref'
        self.manifest = json.loads((ref / 'MANIFEST.json').read_text())
        for name, expected in self.manifest['outputs'].items():
            if sha(ref / name) != expected:
                raise ValueError(f'Frozen reference hash mismatch: {name}')
        self.cop = {r['accession']: r for r in read_tsv(str(ref/'copies.tsv')) if r['status']=='retained'}
        self.cul = {r['culture']: r for r in read_tsv(str(ref/'cultures.tsv'))}
        self.cultures = frozenset(self.cul)
        self.genus = {c: GENUS_SYNONYMS.get(r['genus'],r['genus']) for c,r in self.cul.items()}
        self.species = {c:r['organism'] for c,r in self.cul.items() if is_named(r['organism']) and r.get('name_status','ok')=='ok'}
        self.by_species = defaultdict(set)
        self.by_label = defaultdict(set)
        self.by_culture = defaultdict(list)
        for c, sp in self.species.items(): self.by_species[sp].add(c)
        for c, r in self.cul.items(): self.by_label[r['culture_label']].add(c)
        for x, r in self.cop.items(): self.by_culture[r['culture']].append(x)
        self.fa = {m:read_fasta(str(ref/'fasta'/f'{m}.fasta')) for m in ('ssu_flank','full_its','s58','LSU_LROR-FLR2')}
        self.evidence = {}
        self.legacy = {}
        self.best = {}
        self.verification = []
        bundle = json.loads((self.inputs/'evidence/MANIFEST.json').read_text())
        for item in bundle['files']:
            p = self.inputs/'evidence'/item['path']
            if p.stat().st_size != item['bytes'] or sha(p) != item['sha256']: raise ValueError(f'Evidence digest mismatch: {p.name}')
        prior = json.loads((self.inputs/'baseline/crosswalk_v4.run.json').read_text())
        for arm, name in [('VT','vt_sets_v4'),('SH','sets_SH_v4'),('LSU','sets_LSU_v4')]:
            e = json.load(gzip.open(self.inputs/'evidence'/f'{arm}.evidence.json.gz','rt'))
            if e['schema']!='amf_arm_unit_evidence_v1' or e['arm']!=arm: raise ValueError('Unexpected evidence schema/arm')
            for key, path in [('copies_sha256',ref/'copies.tsv'),('cultures_sha256',ref/'cultures.tsv'),('query_sha256',ref/'fasta'/f"{e['region']}.fasta")]:
                if sha(path)!=e[key]: raise ValueError(f'{arm}: {key} differs')
            if e['reference_sha256']!=prior['arms'][arm]['sha256']: raise ValueError(f'{arm} external release mismatch')
            if e['culture_reference_fingerprint']!=self.manifest['fingerprint']: raise ValueError('Fingerprint mismatch')
            if set(e['queries'])!=set(self.cop): raise ValueError('Query ID mismatch')
            _, region, task, min_id, min_cov, spread = ARMS[arm]
            for k,v in [('region',region),('task',task),('min_id',min_id),('min_cov',min_cov),('max_evalue',1e-50),('default_spread_quantile',spread)]:
                if e[k]!=v: raise ValueError(f'Unexpected evidence setting: {arm} {k}')
            # Source hashes bind scores to the actual original filtering functions.
            for src in ('arms.py','arm_sets.py'):
                if sha(Path(__file__).parent/src)!=e['source_script_sha256'][src]: raise ValueError(f'Filtering code mismatch: {src}')
            rows = read_tsv(str(self.inputs/'baseline'/f'{name}.copies.tsv'))
            self.legacy[arm] = {r['accession']: r for r in rows}
            self.best[arm] = {}
            top = 'top_units' if 'top_units' in rows[0] else 'top_vts'
            ties = 0
            for x,r in self.legacy[arm].items():
                q = e['queries'][x]; ids = q['unit_pident']
                if any(not math.isfinite(p) or not 0<=p<=100 for p in ids.values()): raise ValueError('Invalid identity')
                recorded = float(r['anchor_pident']) if r['anchor_pident'] else None
                if q['anchor_pident'] != recorded: raise ValueError(f'Anchor changed: {arm} {x}')
                if r[top]:
                    u, old_id = r[top].split(',')[0].rsplit(':',1)
                    maximum = max(ids.values())
                    if ids.get(u) != maximum: raise ValueError(f'Recorded best is not full-precision maximum: {arm} {x}')
                    ties += sum(p==maximum for p in ids.values())>1
                    # Preserve the original raw-unit tie choice; this contains no
                    # culture labels or learned mergers. Use full precision now.
                    self.best[arm][x] = (u, ids[u])
                else:
                    if ids: raise ValueError('Missing legacy best despite evidence')
                    self.best[arm][x] = ('',0.)
            self.evidence[arm] = e
            self.verification.append({'arm':arm,'queries':len(e['queries']),'scores':sum(len(q['unit_pident']) for q in e['queries'].values()),'maximum_identity_ties':ties,'hashes_verified':True,'anchors_identical':True,'recorded_best_is_true_maximum':True})
        self.nearest_cache = {}
        self.mains = [x for x,r in self.cop.items() if r['copy_class']=='main' and x in self.fa['ssu_flank']]
        self.main58 = defaultdict(list)
        for x in self.mains:
            if self.cop[x].get('s58_len') and x in self.fa['s58']: self.main58[self.cop[x]['culture']].append(x)
        for x,r in self.cop.items():
            if x in self.fa['s58'] and int(r['s58_len']) != len(self.fa['s58'][x]): raise ValueError('5.8S length metadata mismatch')

    def folds(self, mode):
        grouped = {}
        for c in sorted(self.cultures):
            drop = set(self.by_species[self.species[c]]) if mode.endswith('species') and c in self.species else {c}
            if mode.startswith('physical_'):
                drop = set().union(*(self.by_label[self.cul[o]['culture_label']] for o in drop))
            drop = frozenset(drop)
            grouped.setdefault(drop,[]).extend(self.by_culture[c])
        return [(drop,sorted(xs)) for drop,xs in sorted(grouped.items(),key=lambda p:sorted(p[0]))]

    def truth_group(self, c, drop):
        sp = self.species.get(c)
        if not sp: return "unnamed ('sp.')"
        return 'named, species elsewhere' if self.by_species[sp]-drop else 'named, species held out'

    def nearest(self, x, drop):
        """SSU nearest training main copies, including all ties; cached per fold.

        No merged unit or taxonomic call enters the distance calculation.
        The cache includes the complete excluded-culture set in its key.
        """
        key = (x, frozenset(drop))
        if key in self.nearest_cache: return self.nearest_cache[key]
        from novelty_baseline_consistent import dist as bounded_dist
        near_d, nearest = None, []
        fa = self.fa['ssu_flank']
        if x in fa:
            # Raw sequence similarity only affects evaluation order, not eligibility.
            q=fa[x]
            candidates=[y for y in self.mains if self.cop[y]['culture'] not in drop]
            candidates.sort(key=lambda y:(self.best['VT'][y][0]!=self.best['VT'][x][0],self.cop[y]['culture'],y))
            for y in candidates:
                k=-1 if near_d is None else int(near_d*min(len(q),len(fa[y]))/100)+1
                d=bounded_dist(q,fa[y],k)
                if d is None: continue
                if near_d is None or d<near_d-1e-12: near_d,nearest=d,[y]
                elif math.isclose(d,near_d,rel_tol=0,abs_tol=1e-12): nearest.append(y)
        nearest.sort(key=lambda y:(self.cop[y]['culture'],y))
        result=(nearest[0] if nearest else None,near_d,nearest)
        self.nearest_cache[key]=result
        return result


class ArmTrainer:
    def __init__(self,data,arm):
        self.data,self.arm=data,arm
        self.e=data.evidence[arm]
        self.min_id=self.e['min_id']
        self.region=self.e['region']
        qseq=data.fa[self.region]
        self.main_of={c:[x for x in data.by_culture[c] if data.cop[x]['copy_class']=='main' and x in qseq] for c in data.cultures}
        self.raw={d:{x:vt_set(q['unit_pident'],q['anchor_pident'],d,self.min_id,{}) for x,q in self.e['queries'].items()} for d in GRID}
        self.spreads={c:[dist(qseq[p],qseq[q]) for p,q in combinations(xs,2)] for c,xs in self.main_of.items()}
        self.proposals={}
        self.rejected={}
        ssu=data.fa['ssu_flank']
        # Culture-local summaries only; fit() filters whole cultures before
        # combining proposals, quantiles, consistency, or species links.
        for c,xs in self.main_of.items():
            side=defaultdict(list)
            for x in xs:
                s=self.raw[0.][x]
                if s: side[s].append(x)
            if len(side)<2 or frozenset.intersection(*side): continue
            gap=max(min(dist(ssu[p],ssu[q]) if p in ssu and q in ssu else 100. for p in side[a] for q in side[b]) for a,b in combinations(side,2))
            row={'culture':c,'units':sorted(frozenset().union(*side)),'ssu_gap':gap}
            (self.proposals if gap<=1. else self.rejected)[c]=row
        self.fit_cache={}

    def fit(self, training, fixed_delta=None):
        training=frozenset(training)
        key=(training,fixed_delta)
        if key in self.fit_cache: return self.fit_cache[key]
        if not training <= self.data.cultures or not training: raise ValueError('Invalid training cultures')
        proposals=[self.proposals[c] for c in sorted(training) if c in self.proposals]
        unit_map={}
        for g in union_find([p['units'] for p in proposals]):
            if len(g)>1:
                name='+'.join(sorted(g));unit_map.update({u:name for u in g})
        spread=[d for c in sorted(training) for d in self.spreads[c]]
        qs=statistics.quantiles(spread,n=100) if len(spread)>1 else [0.]*99
        quantile=self.e['default_spread_quantile']
        d_spread=round(min(math.ceil(round(qs[quantile-1]/.1,6))*.1,3.),3) if quantile else 0.
        cal=[]
        for d in GRID:
            tested=consistent=0
            for c in training:
                ss=[frozenset(unit_map.get(u,u) for u in self.raw[d][x]) for x in self.main_of[c] if self.raw[d][x]]
                if len(ss)>=2:
                    tested+=1;consistent+=bool(frozenset.intersection(*ss))
            cal.append({'delta':d,'cultures_tested':tested,'consistent':consistent,'consistent_share':round(consistent/tested,4) if tested else ''})
        ok=[r['delta'] for r in cal if r['consistent_share']!='' and r['consistent_share']>=.99]
        d_cons=ok[0] if ok else GRID[-1]
        delta=max(d_cons,d_spread) if fixed_delta is None else fixed_delta
        sets={x:frozenset(unit_map.get(u,u) for u in us) for x,us in self.raw[delta].items()}
        culsets={}
        holders=defaultdict(set)
        for c in sorted(training):
            ss=[sets[x] for x in self.main_of[c] if sets[x]]
            inter=frozenset.intersection(*ss) if ss else frozenset()
            culsets[c]=inter or frozenset().union(*ss)
            for u in culsets[c]: holders[u].add(c)
        best=self.data.best[self.arm]
        anchors={x:next((u for u in sets[x] if best[x][0] in u.split('+')),'') for x in sets}
        result=Fit(self.arm,training,unit_map,delta,d_cons,d_spread,cal,sets,culsets,dict(holders),anchors,best,proposals)
        # Avoid retaining all training-set copies for hundreds of folds.
        return result


def reconcile(data, x, fits, drop):
    r=data.cop[x];c=r['culture']; per={}
    for n,f in fits.items():
        units=f.sets[x]
        others=set().union(*(f.holders.get(u,set()) for u in units))-drop
        sp={data.species[o] for o in others if o in data.species}
        gn={data.genus[o] for o in others if data.genus[o]}
        per[n]=('no call' if not units else 'matched' if gn else 'unmatched',sp,gn)
    matched=[n for n in per if per[n][0]=='matched']
    if not matched: rank,sp,gn='unplaced',set(),set()
    else:
        gn=set.intersection(*(per[n][2] for n in matched))
        sn=[n for n in matched if per[n][1]]
        sp=set.intersection(*(per[n][1] for n in sn)) if sn else set()
        sp={s for s in sp if any(data.genus[o] in gn for o in data.by_species[s]-drop)} if gn else set()
        rank='discordant' if not gn else 'species' if len(sp)==1 else 'species group' if sp else 'genus'
    row={'accession':x,'culture':c,'organism':r['organism'],'genus':data.genus[c],'copy_class':r['copy_class'],'div_kind':r.get('div_kind',''),'group':data.truth_group(c,drop)}
    for n in ('VT','SH','LSU'):
        st,s,g=per[n];row.update({f'{n}_state':st,f'{n}_species':';'.join(sorted(s)),f'{n}_genera':';'.join(sorted(g)),f'{n}_units':','.join(sorted(fits[n].sets[x]))})
    true=data.species.get(c,'')
    correct='yes' if rank=='species' and sp=={true} else 'contains' if rank=='species group' and true in sp else 'genus' if rank=='genus' and data.genus[c] in gn else 'genus only' if rank in ('species','species group') and data.genus[c] in gn else 'no' if rank in ('species','species group','genus') else ''
    row.update(rank=rank,call_species=';'.join(sorted(sp)),call_genera=';'.join(sorted(gn)),correct=correct,unmatched_arms=';'.join(n for n in ('VT','SH','LSU') if per[n][0]=='unmatched'))
    return row
