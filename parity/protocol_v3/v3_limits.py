#!/usr/bin/env python3
"""Within-species and within-genus limits on culture reference v4, for protocol v3 (draft).

Two within-species pools are frozen in release 0.3.1 and both are reproduced here:
  - N6 pool (stage3/its_rule.py): every main-type copy of a named species with >= 2 cultures;
    its distance to the nearest main-type copy of another culture of the same species.
  - complex pool (paper/analysis/species_groups.py): unique (sequence, culture) pairs per species.
The within-genus limit of R3a uses unique (sequence, species) pairs per genus: the distance to the
nearest pair of another named species of the same genus, over genera with >= 2 named species.
Percentiles: statistics.quantiles exclusive method (species_groups.percentile for 97.5).
Distance: edlib infix (HW) edit distance, % of the shorter sequence, first argument as query on ties.

Usage: python3 v3_limits.py --repo /path/to/amf-novelty-benchmark [--out v3_limits.tsv]
Exits 1 if the complex-pool limits differ from the paper's (5.98803, 3.02698, 0.80652; NS31-AM1 1.7
to one decimal). The NS31-AM1 limit is the radius measurement E uses for MaarjAM records.
"""
import argparse, csv, os, statistics, sys
from collections import defaultdict
import edlib

PAPER = {"ITS": 5.98803, "LSU": 3.02698, "SSU": 0.80652, "NS31-AM1": 1.7}  # NS31-AM1: Methods S8, one decimal
MARKERS = (("ITS", "full_its"), ("LSU", "LSU_LROR-FLR2"), ("SSU", "ssu_flank"), ("NS31-AM1", "SSU_VT_NS31-AM1"))


def dist(a, b):
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    return 100 * edlib.align(q, t, mode="HW", task="distance")["editDistance"] / len(q)


def pct(vals, p):  # species_groups.percentile: exclusive method, any real p
    d = sorted(vals); n = len(d); h = (n + 1) * p / 100.0
    j = min(max(int(h), 1), n - 1)
    return d[j - 1] + (d[j] - d[j - 1]) * (h - j)


def fasta(path):
    seqs, key = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            key = line[1:].split()[0]; seqs[key] = []
        elif key:
            seqs[key].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--repo", required=True)
    ap.add_argument("--out", default="v3_limits.tsv")
    a = ap.parse_args()
    sys.path.insert(0, os.path.join(a.repo, "stage3"))
    from arms import is_named
    ref = os.path.join(a.repo, "reference", "culture_ref_v4")
    cul = {r["culture"]: r for r in csv.DictReader(open(os.path.join(ref, "cultures.tsv")), delimiter="\t")}
    sp_of = {c: m["organism"] for c, m in cul.items() if is_named(m["organism"]) and m.get("name_status", "ok") == "ok"}
    cop = [r for r in csv.DictReader(open(os.path.join(ref, "copies.tsv")), delimiter="\t")
           if r["status"] == "retained" and r["copy_class"] == "main" and r["culture"] in sp_of]
    genus_of = {sp_of[r["culture"]]: r["genus"] for r in cop}
    rows, ok = [], True
    for lab, fn in MARKERS:
        seq = fasta(os.path.join(ref, "fasta", fn + ".fasta"))
        by_sp = defaultdict(list)
        for r in cop:
            if r["accession"] in seq:
                by_sp[sp_of[r["culture"]]].append(r)
        # complex pool (species_groups.py)
        uniq = {sp: sorted({(seq[r["accession"]], r["culture"]) for r in rs}) for sp, rs in by_sp.items()}
        cx = []
        for sp in sorted(uniq):
            items = uniq[sp]
            if len({c for _, c in items}) < 2:
                continue
            for s1, c1 in items:
                cx.append(min(dist(s1, s2) for s2, c2 in items if c2 != c1))
        # N6 pool (its_rule.py), ITS only: every main copy, pairs in sorted accession order
        n6 = None
        if lab == "ITS":
            acc = sorted(r["accession"] for rs in by_sp.values() for r in rs)
            culture = {r["accession"]: r["culture"] for rs in by_sp.values() for r in rs}
            species = {r["accession"]: sp_of[r["culture"]] for rs in by_sp.values() for r in rs}
            near = defaultdict(dict)
            by_g = defaultdict(list)
            for x in acc:
                by_g[genus_of[species[x]]].append(x)
            for xs in by_g.values():  # same-species pairs are always same-genus
                for i, x in enumerate(xs):
                    for y in xs[i + 1:]:
                        if species[x] != species[y]:
                            continue
                        d = dist(seq[x], seq[y])
                        for q, r in ((x, y), (y, x)):
                            c = culture[r]
                            if c != culture[q] and d < near[q].get(c, float("inf")):
                                near[q][c] = d
            vals = []
            sp_cult = defaultdict(set)
            for x in acc:
                sp_cult[species[x]].add(culture[x])
            for x in acc:
                if len(sp_cult[species[x]]) >= 2 and near[x]:
                    vals.append(min(near[x].values()))
            n6 = statistics.quantiles(vals, n=100)[98]
        # within-genus pool (R3a)
        gp = defaultdict(set)
        for sp, items in uniq.items():
            for s, _ in items:
                gp[genus_of[sp]].add((s, sp))
        gv = []
        for g in sorted(gp):
            items = sorted(gp[g])
            if len({sp for _, sp in items}) < 2:
                continue
            for s1, sp1 in items:
                gv.append(min(dist(s1, s2) for s2, sp2 in items if sp2 != sp1))
        t = pct(cx, 99)
        good = abs(t - PAPER[lab]) < 5e-6 if lab != "NS31-AM1" else round(t, 1) == PAPER[lab]
        ok &= good
        rows.append({"marker": lab, "within_species_complex_pool_p99": round(t, 5),
                     "paper": PAPER[lab], "matches_paper": good,
                     "within_species_n6_pool_p99": "" if n6 is None else round(n6, 5),
                     "within_genus_p95": round(pct(gv, 95), 3), "within_genus_p97_5": round(pct(gv, 97.5), 3),
                     "within_genus_p99": round(pct(gv, 99), 3), "n_complex_pool": len(cx), "n_genus_pool": len(gv)})
        print(rows[-1])
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
