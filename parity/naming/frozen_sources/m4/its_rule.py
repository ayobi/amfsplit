"""Culture-ITS novelty rule (N6), fitted on training cultures only.

SH is silent for a third of copies because UNITE has no close sequence, so a held-out species that
shares its VT and LSU units with a sister species is called that sister species (G. chinense ->
G. rugosae). The culture reference holds the ITS of every culture, so it can answer the question SH
cannot: is this copy's ITS like the ITS of any training culture?

Threshold, per fold, from training cultures only: for every main-type copy of a named species with
>= 2 training cultures, the full-ITS distance to the nearest main-type copy of another training culture
of the same species. T is the --quantile percentile of those distances (default 99th).

Query statistic: the full-ITS distance to the nearest copy of any training culture. By default the
pool holds main-type and divergent-class copies (pool='all'), so a divergent copy the D1/D2 screen
missed still finds its paralog elsewhere instead of looking new; pool='main' uses main-type copies only.

Rule: a 'known species' or 'known species group' status becomes 'novel species, genus known' when the
query statistic exceeds T. Nothing else changes. Distances are % edit distance with the shorter
sequence aligned inside the longer (the edlib infix distance used everywhere in the project).
"""
from collections import defaultdict
import statistics

from arm_sets import dist


class ItsRule:
    def __init__(self, data, quantile=99):
        self.data, self.quantile = data, quantile
        its = data.fa['full_its']
        self.its_ok = sorted(x for x in data.cop if x in its)
        self.mains = {x for x in self.its_ok if data.cop[x]['copy_class'] == 'main'}
        # nearest copy of each culture to each query: pool 'all' (every copy) and 'main' (main-type only)
        self.dmin = {'all': defaultdict(dict), 'main': defaultdict(dict)}
        for i, x in enumerate(self.its_ok):
            for y in self.its_ok[i + 1:]:
                d = dist(its[x], its[y])
                for q, r in ((x, y), (y, x)):
                    c = data.cop[r]['culture']
                    if d < self.dmin['all'][q].get(c, float('inf')):
                        self.dmin['all'][q][c] = d
                    if r in self.mains and d < self.dmin['main'][q].get(c, float('inf')):
                        self.dmin['main'][q][c] = d
        self._thr = {}

    def threshold(self, training):
        """(T, n): quantile of nearest same-species other-culture distance, main-type copies of named
        species with >= 2 training cultures. Training cultures only."""
        training = frozenset(training)
        if training in self._thr:
            return self._thr[training]
        data, vals, used = self.data, [], set()
        for sp, cs in data.by_species.items():
            tc = cs & training
            if len(tc) < 2:
                continue
            for c in tc:
                for y in data.by_culture[c]:
                    if y not in self.mains:
                        continue
                    ds = [self.dmin['main'][y][o] for o in tc - {c} if o in self.dmin['main'][y]]
                    if ds:
                        vals.append(min(ds))
                        used.add(c)
        if not used <= training:
            raise AssertionError('held-out culture entered the ITS threshold')
        t = statistics.quantiles(vals, n=100)[self.quantile - 1] if len(vals) >= 2 else None
        self._thr[training] = (t, len(vals))
        return self._thr[training]

    def nearest(self, x, training, pool='all'):
        own = self.data.cop[x]['culture']
        ds = [d for c, d in self.dmin[pool].get(x, {}).items() if c in training and c != own]
        return min(ds) if ds else None

    def apply(self, rows, rec, training, pool='all'):
        """Post-process decision-tree rows in place; returns the number demoted."""
        t, n = self.threshold(training)
        fired = 0
        for r in rows:
            x = r['accession']
            d = self.nearest(x, training, pool)
            r['its_nearest_training_pct'] = round(d, 3) if d is not None else ''
            r['its_threshold_pct'] = round(t, 3) if t is not None else ''
            r['its_threshold_n'] = n
            r['its_pool'] = pool
            hit = d is not None and t is not None and d > t and r['status'].startswith('known')
            r['its_rule_fired'] = int(hit)
            if hit:
                fired += 1
                r['status'] = 'novel species, genus known'
                r['placement'] = rec[x]['call_genera']
                r['novelty_claim'] = 'candidate_only'
                r['reasons'] = '; '.join(filter(None, [r['reasons'], (
                    f"N6 full ITS {d:.1f}% from the nearest training culture, above {t:.1f}% "
                    f"(the {self.quantile}th percentile of within-species distance between training cultures)")]))
        return fired
