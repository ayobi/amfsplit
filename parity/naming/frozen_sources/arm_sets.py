#!/usr/bin/env python3
"""Phase 2: compatible unit sets for any arm (VT, SH or LSU), with a calibrated margin δ.

Supersedes vt_sets.py: `--arm VT` gives the same result. Arms and defaults:
  VT   SSU flank  vs MaarjAM               megablast  >= 97%    >= 95%   spread floor 95th percentile
  SH   full ITS   vs UNITE (Glomeromycota) blastn     >= 98.5%  >= 90%   no spread floor (see below)
  LSU  LROR-FLR2  vs Delavaux V18          blastn     >= 97%    >= 90%   no spread floor (see below)
e-value <= 1e-50 throughout; coverage is of the shorter sequence with reference Ns not counted; hits
to culture-reference accessions are dropped.

Each copy gets a set instead of a best hit. The anchor is the copy's best hit (highest bit score
among hits passing identity, coverage and e-value). The set holds every unit whose best identity is
at least (anchor identity - δ) and passes the identity threshold, so the best-hit unit is always in
it. A copy with no passing hit gets an empty set (unassigned).

Merged units (--merge auto, the default): units that split the main-type copies of one culture at
δ = 0 are joined (e.g. VTX00113+VTX00115) when the closest copies across the split are within
--merge-max-dist in the SSU flank (default 1.0%, about the 99th percentile of within-culture SSU
spread in the Stefani set). The SSU check is used on every arm: it asks whether the two sides are one
genome, which a mixed culture (such as an unsplit BR105_2435A) fails. On the SH arm the merged units
are the SH groups of a genome.

Calibration, on main-type copies. A culture with >= 2 main-type copies that have a set is consistent
when those sets share a unit. For each δ the report shows consistency and the cost: copies left
with one unit, mean set size, and named species still resolved (no other named species shares any
of their units). The chosen δ is the larger of
  (a) the smallest δ at which >= --target of cultures are consistent (default 0.99), and
  (b) the --spread-quantile percentile of within-culture spread in the arm's own region, rounded up:
      two copies d% apart can differ by up to d points in identity to any reference.
The SH and LSU arms have no spread floor by default. Within-genome ITS and LSU spread is wide and
uneven across genomes (Stefani v4, 95th percentile: ITS 8.6%, LROR-FLR2 3.9%, driven by genomes such
as Gigaspora rosea whose rDNA types differ by 6% in LSU), so a global floor would inflate every set to
cover a few genomes. The merged units carry each genome's observed variation instead. SSU spread is
small and even (95th percentile 0.54%), so the VT arm keeps its floor. --delta fixes δ.

Reports:
  1. within-culture spread in the arm's region
  2. merged units
  3. calibration table and the chosen δ
  4. at the chosen δ: multi-unit culture sets, cultures still inconsistent, units shared by species
  5. divergent-class copies: does their set contain their culture's unit?
Writes {out}.copies.tsv, {out}.cultures.tsv and {out}.calibration.tsv.

Usage (from amfsplit_stage3):
  python3 arm_sets.py --arm VT --ref ../ref/culture_ref/v4 --reference ../ref/maarjam/maarjam.fasta \\
      --out ../bench/stefani/sets_VT_v4 --threads 8
  python3 arm_sets.py --arm SH --ref ../ref/culture_ref/v4 \\
      --reference ../ref/unite/UNITE_public_19.02.2025.fasta.gz \\
      --crosswalk ../bench/stefani/crosswalk_v4.copies.tsv --out ../bench/stefani/sets_SH_v4 --threads 8
"""
import argparse
import collections
import itertools
import math
import os
import shutil
import statistics
import sys

import edlib

from arms import base_acc, is_named, load_reference, read_fasta, read_tsv, run_blast, write_tsv


def dist(a, b):
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    return 100 * edlib.align(q, t, mode="HW", task="distance")["editDistance"] / len(q)


ARMS = {  # arm: loader, query region, BLAST task, min identity, min coverage, default spread quantile
    "VT": ("maarjam", "ssu_flank", "megablast", 97.0, 0.95, 95),
    "SH": ("unite", "full_its", "blastn", 98.5, 0.90, 0),
    "LSU": ("delavaux", "LSU_LROR-FLR2", "blastn", 97.0, 0.90, 0),
}


def vt_hits(hs, refs, min_cov, max_ev, min_id, drop):
    """({unit: best identity}, anchor identity): identities over hits passing coverage and e-value; the
    anchor is the identity of the highest-bit-score hit that also passes min_id (None if none does)."""
    ids, best = {}, None
    for h in hs:
        info = refs[h["ref"]]
        if base_acc(info["acc"]) in drop:
            continue
        cov = h["length"] / min(h["qlen"], max(1, h["slen"] - info.get("n_N", 0)))
        if cov < min_cov or h["evalue"] > max_ev:
            continue
        v = info["label"]
        if h["pident"] > ids.get(v, -1):
            ids[v] = h["pident"]
        if h["pident"] >= min_id and (best is None or (h["bitscore"], h["pident"]) > best):
            best = (h["bitscore"], h["pident"])
    return ids, (best[1] if best else None)


def vt_set(ids, anchor, delta, min_id, unit):
    if anchor is None:
        return frozenset()
    thr = max(min_id, anchor - delta) - 1e-9
    return frozenset(unit.get(v, v) for v, p in ids.items() if p >= thr)


def union_find(groups):
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for g in groups:
        g = sorted(g)
        for v in g:
            find(v)
        for v in g[1:]:
            parent[find(v)] = find(g[0])
    out = collections.defaultdict(set)
    for v in parent:
        out[find(v)].add(v)
    return list(out.values())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm", required=True, choices=list(ARMS))
    ap.add_argument("--ref", required=True, help="frozen culture reference directory")
    ap.add_argument("--reference", "--maarjam", dest="reference", required=True,
                    help="MaarjAM (VT), UNITE (SH) or Delavaux (LSU) FASTA")
    ap.add_argument("--unite-keep", default="Glomeromycota", help="SH arm: keep UNITE records matching this regex")
    ap.add_argument("--crosswalk", help="crosswalk copies table, to compare SH/LSU sets with the best-hit calls")
    ap.add_argument("--min-id", type=float, help="default per arm (VT 97, SH 98.5, LSU 97)")
    ap.add_argument("--min-cov", type=float, help="default per arm (VT 0.95, SH 0.90, LSU 0.90)")
    ap.add_argument("--max-evalue", type=float, default=1e-50)
    ap.add_argument("--merge", choices=["auto", "none"], default="auto",
                    help="join VTs that split one culture's main-type copies at δ = 0 (default auto)")
    ap.add_argument("--merge-max-dist", type=float, default=1.0,
                    help="max SSU %% distance between the closest copies across a split for it to merge (default 1.0)")
    ap.add_argument("--target", type=float, default=0.99, help="share of cultures that must be consistent")
    ap.add_argument("--spread-quantile", type=int,
                    help="percentile of within-culture spread in the arm's region that δ must cover "
                         "(0 = ignore; default VT 95, SH 0, LSU 0)")
    ap.add_argument("--max-delta", type=float, default=3.0)
    ap.add_argument("--step", type=float, default=0.1)
    ap.add_argument("--delta", type=float, help="use this δ instead of calibrating")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    kind, region, task, min_id, min_cov, spread_q = ARMS[a.arm]
    a.min_id = min_id if a.min_id is None else a.min_id
    a.min_cov = min_cov if a.min_cov is None else a.min_cov
    a.spread_quantile = spread_q if a.spread_quantile is None else a.spread_quantile

    rows = read_tsv(os.path.join(a.ref, "copies.tsv"))
    cop = {r["accession"]: r for r in rows if r["status"] == "retained"}
    drop = {base_acc(r["accession"]) for r in rows}
    cul_meta = {r["culture"]: r for r in read_tsv(os.path.join(a.ref, "cultures.tsv"))}
    ssu = read_fasta(os.path.join(a.ref, "fasta", "ssu_flank.fasta"))
    qpath = os.path.join(a.ref, "fasta", f"{region}.fasta")
    qseq = read_fasta(qpath)
    refs, fasta, _ = load_reference(a.reference, kind, a.unite_keep if kind == "unite" else None)
    hits = run_blast(qpath, fasta, a.out + ".blastdb", a.arm.lower(), task, a.threads)
    shutil.rmtree(a.out + ".blastdb", ignore_errors=True)
    hv = {x: vt_hits(hits.get(x, []), refs, a.min_cov, a.max_evalue, a.min_id, drop) for x in cop}

    by_cul = collections.defaultdict(list)
    for x, r in cop.items():
        by_cul[r["culture"]].append(x)
    main_of = {cu: [x for x in xs if cop[x]["copy_class"] == "main" and x in qseq] for cu, xs in by_cul.items()}
    names = {cu: cul_meta[cu]["organism"] for cu in by_cul}
    named = {cu for cu in by_cul if is_named(names[cu]) and cul_meta[cu].get("name_status", "ok") == "ok"}
    taxon = collections.defaultdict(collections.Counter)
    for v in refs.values():
        taxon[v["label"]][v["taxon"]] += 1
    show_u = (lambda u: "+".join(f"{p} ({taxon[p].most_common(1)[0][0]})" for p in u.split("+"))) \
        if a.arm == "SH" else (lambda u: u)
    print(f"arm {a.arm}: {len(cop)} copies, {len(by_cul)} cultures; {os.path.basename(a.reference)} "
          f"{len(refs)} sequences, {len({v['label'] for v in refs.values()})} units; query {region}")

    # ---- 1. within-culture spread in the arm's region
    spread = []
    for cu, xs in main_of.items():
        spread += [dist(qseq[p], qseq[q]) for p, q in itertools.combinations(xs, 2)]
    qs = statistics.quantiles(spread, n=100) if len(spread) > 1 else [0.0] * 99
    print(f"\n1. within-culture {region} spread, main-type copy pairs (n = {len(spread)}): median "
          f"{statistics.median(spread):.2f}%, 90th {qs[89]:.2f}%, 95th {qs[94]:.2f}%, 99th {qs[98]:.2f}%, "
          f"max {max(spread):.2f}%")

    # ---- 2. merged units from within-culture splits at δ = 0
    unit, why, kept_split = {}, collections.defaultdict(list), []
    merge_max = a.merge_max_dist
    if a.merge == "auto":
        splits = []
        for cu, xs in main_of.items():
            side = collections.defaultdict(list)
            for x in xs:
                s0 = vt_set(*hv[x], 0.0, a.min_id, {})
                if s0:
                    side[s0].append(x)
            if len(side) < 2 or frozenset.intersection(*side):
                continue
            gap = max(min((dist(ssu[p], ssu[q]) if p in ssu and q in ssu else 100.0)
                          for p in side[s1] for q in side[s2])
                      for s1, s2 in itertools.combinations(side, 2))
            if gap <= merge_max:
                splits.append((cu, frozenset().union(*side), gap))
            else:
                kept_split.append((cu, sorted("/".join(sorted(k)) for k in side), gap))
        for g in union_find([s for _, s, _ in splits]):
            if len(g) > 1:
                name = "+".join(sorted(g))
                for v in g:
                    unit[v] = name
        for cu, s, gap in splits:
            why[unit[next(iter(s))]].append(f"{cu} ({gap:.2f}%)")
    print(f"\n2. merged units: units that split one culture's main-type copies at δ = 0, closest copies across "
          f"the split <= {merge_max:.2f}% apart in the SSU flank" + (": none" if not why else ""))
    for u, cus in sorted(why.items(), key=lambda kv: (-len(kv[1]), kv[0]))[:15]:
        print(f"      {show_u(u)}: {', '.join(sorted(cus)[:4])}" + (f" +{len(cus) - 4} more" if len(cus) > 4 else ""))
    if len(why) > 15:
        print(f"      ... {len(why) - 15} more merged units (listed per culture in {a.out}.cultures.tsv)")
    for cu, sides, gap in kept_split[:10]:
        print(f"      not merged, {gap:.2f}% apart in SSU (mixed culture?): {cu}: {' vs '.join(sides)}")
    if len(kept_split) > 10:
        print(f"      ... {len(kept_split) - 10} more splits not merged")

    # ---- 3. calibration
    grid = [round(i * a.step, 3) for i in range(int(round(a.max_delta / a.step)) + 1)]
    if a.delta is not None and a.delta not in grid:
        grid = sorted(grid + [a.delta])
    all_named_species = {names[cu] for cu in named}
    cal, state = [], {}
    for d in grid:
        sets = {x: vt_set(*hv[x], d, a.min_id, unit) for x in cop}
        testable = consistent = 0
        culset = {}
        for cu, xs in main_of.items():
            ss = [sets[x] for x in xs if sets[x]]
            if not ss:
                culset[cu] = frozenset()
                continue
            inter = frozenset.intersection(*ss)
            culset[cu] = inter if inter else frozenset().union(*ss)
            if len(ss) >= 2:
                testable += 1
                consistent += bool(inter)
        filled = [s for s in sets.values() if s]
        groups = union_find([s for s in culset.values() if s])
        group_of = {u: i for i, g in enumerate(groups) for u in g}
        g_species = collections.defaultdict(set)
        for cu in named:
            for u in culset[cu]:
                g_species[group_of[u]].add(names[cu])
        sp_groups = collections.defaultdict(set)
        for cu in named:
            for u in culset[cu]:
                sp_groups[names[cu]].add(group_of[u])
        resolved = sum(1 for sp, gs in sp_groups.items() if all(len(g_species[g]) == 1 for g in gs))
        row = {"delta": d, "cultures_tested": testable, "consistent": consistent,
               "consistent_share": round(consistent / testable, 4) if testable else "",
               "copies_with_set": len(filled),
               "single_unit_share": round(sum(len(s) == 1 for s in filled) / len(filled), 4) if filled else "",
               "mean_set_size": round(statistics.mean(len(s) for s in filled), 3) if filled else "",
               "unit_groups": len(groups), "groups_with_several_species": sum(len(s) > 1 for s in g_species.values()),
               "species_resolved": resolved, "named_species": len(all_named_species)}
        cal.append(row)
        state[d] = (sets, culset, groups, g_species, group_of)
    ok = [r["delta"] for r in cal if r["consistent_share"] != "" and r["consistent_share"] >= a.target]
    d_cons = ok[0] if ok else grid[-1]
    d_spread = (math.ceil(round(qs[a.spread_quantile - 1] / a.step, 6)) * a.step) if a.spread_quantile else 0.0
    d_spread = round(min(d_spread, grid[-1]), 3)
    chosen = a.delta if a.delta is not None else max(d_cons, d_spread)
    print(f"\n3. calibration: a culture is consistent when its main-type copies' sets share a unit")
    print(f"  {'δ':>5s}{'consistent':>14s}{'one unit':>10s}{'mean set':>10s}{'species resolved':>18s}")
    show = {0.0, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, d_cons, d_spread, chosen}
    for r in cal:
        if r["delta"] in show:
            tags = [t for t, v in (("consistency", d_cons), ("spread", d_spread), ("chosen", chosen)) if r["delta"] == v]
            share = 100 * r["consistent_share"] if r["consistent_share"] != "" else float("nan")
            print(f"  {r['delta']:>5.1f}{r['consistent']:>6d}/{r['cultures_tested']:<3d}{share:>4.0f}%"
                  f"{100 * r['single_unit_share']:>9.0f}%{r['mean_set_size']:>10.2f}"
                  f"{r['species_resolved']:>11d} of {r['named_species']:<3d}" + (f"  <- {', '.join(tags)}" if tags else ""))
    if a.delta is not None:
        print(f"  chosen δ = {chosen:g} (fixed by --delta)")
    else:
        print(f"  consistency: smallest δ with >= {100 * a.target:.0f}% of cultures consistent = {d_cons:g}; "
              + (f"spread: {a.spread_quantile}th percentile of within-culture {region} spread, rounded up = {d_spread:g}"
                 if a.spread_quantile else "no spread floor"))
        print(f"  chosen δ = {chosen:g} (the larger)")
    sets0 = state[grid[0]][0]
    if a.arm == "VT":
        call = {x: (cop[x]["vt_call"], cop[x]["vt"]) for x in cop}
        what = "Stage 1 best-hit VT"
    elif a.crosswalk:
        cw = {r["accession"]: r for r in read_tsv(a.crosswalk)}
        call = {x: (cw[x][f"{a.arm}_call"], cw[x][f"{a.arm}_label"]) for x in cop if x in cw}
        what = f"crosswalk best-hit {a.arm}"
    else:
        call, what = {}, ""
    if call:
        n_as = sum(c == "assigned" for c, _ in call.values())
        inset = sum(c == "assigned" and unit.get(u, u) in sets0[x] for x, (c, u) in call.items())
        print(f"  at δ = 0 the {what} is in the set for {inset} of {n_as} assigned copies")

    # ---- 4. at the chosen δ
    sets, culset, groups, g_species, group_of = state[chosen]
    print(f"\n4. at δ = {chosen:g}")
    multi = collections.defaultdict(list)
    for cu, s in culset.items():
        if len(s) > 1:
            multi[s].append(cu)
    print(f"  culture sets with > 1 unit: {sum(len(v) for v in multi.values())} cultures in {len(multi)} distinct sets")
    for s, cus in sorted(multi.items(), key=lambda kv: -len(kv[1]))[:10]:
        orgs = collections.Counter(names[c] for c in cus)
        print(f"      {'/'.join(sorted(s))[:110]}: {', '.join(f'{o} ({n})' for o, n in orgs.most_common())}")
    bad = [cu for cu, xs in main_of.items() if len([x for x in xs if sets[x]]) >= 2
           and not frozenset.intersection(*[sets[x] for x in xs if sets[x]])]
    print(f"  cultures still inconsistent: {len(bad)}")
    for cu in sorted(bad)[:10]:
        cnt = collections.Counter("/".join(sorted(sets[x])) or "-" for x in main_of[cu])
        print(f"      {cu}: " + ", ".join(f"{k} x{n}" for k, n in cnt.most_common()))
    shared = [(i, sp) for i, sp in g_species.items() if len(sp) > 1]
    print(f"  unit groups shared by > 1 named species: {len(shared)} of {len(groups)}")
    for i, sp in sorted(shared, key=lambda t: -len(t[1]))[:10]:
        g = sorted(groups[i])
        print(f"      {'/'.join(g[:6])}{f' +{len(g) - 6} more' if len(g) > 6 else ''}: {', '.join(sorted(sp))}")

    # ---- 5. divergent-class copies
    print(f"\n5. non-main copies at δ = {chosen:g}: does the copy's set contain a unit of its culture's set?")
    for kind in ("short_5.8S", "same_length", "longer", ""):
        xs = [x for x in cop if cop[x]["copy_class"] == "divergent" and cop[x].get("div_kind", "") == kind]
        if not xs:
            continue
        yes = sum(bool(sets[x] & culset[cop[x]["culture"]]) for x in xs)
        empty = sum(not sets[x] for x in xs)
        nocul = sum(bool(sets[x]) and not culset[cop[x]["culture"]] for x in xs)
        other = [x for x in xs if sets[x] and culset[cop[x]["culture"]] and not sets[x] & culset[cop[x]["culture"]]]
        print(f"  divergent ({kind or 'no 5.8S length'}): {len(xs)} copies; shares a culture unit {yes}, "
              f"other units only {len(other)}, no unit {empty}, culture has no unit {nocul}")
        for x in other[:6]:
            print(f"      {x} {cop[x]['culture']}: {'/'.join(sorted(sets[x]))} (culture {'/'.join(sorted(culset[cop[x]['culture']]))})")

    # ---- tables
    crow = []
    for x in sorted(cop):
        ids, anchor = hv[x]
        top = sorted(ids.items(), key=lambda kv: -kv[1])[:5]
        crow.append({"accession": x, "culture": cop[x]["culture"], "organism": cop[x]["organism"],
                     "copy_class": cop[x]["copy_class"], "div_kind": cop[x].get("div_kind", ""),
                     "best_hit_call": call.get(x, ("", ""))[0], "best_hit_unit": call.get(x, ("", ""))[1],
                     "anchor_pident": anchor if anchor is not None else "",
                     "unit_set": ",".join(sorted(sets[x])), "set_size": len(sets[x]),
                     "culture_set": ",".join(sorted(culset[cop[x]["culture"]])),
                     "top_units": ",".join(f"{v}:{p:.2f}" for v, p in top)})
    write_tsv(f"{a.out}.copies.tsv", crow, list(crow[0]))
    urow = []
    for cu in sorted(by_cul):
        xs = [x for x in main_of[cu] if sets[x]]
        urow.append({"culture": cu, "organism": names[cu], "named": cu in named, "main_with_set": len(xs),
                     "culture_set": ",".join(sorted(culset[cu])),
                     "consistent": (bool(frozenset.intersection(*[sets[x] for x in xs])) if len(xs) >= 2 else ""),
                     "copy_sets": ";".join(f"{k}:{n}" for k, n in collections.Counter(
                         "/".join(sorted(sets[x])) for x in xs).most_common())})
    write_tsv(f"{a.out}.cultures.tsv", urow, list(urow[0]))
    write_tsv(f"{a.out}.calibration.tsv", cal, list(cal[0]))
    print(f"\ntables: {a.out}.copies.tsv, .cultures.tsv, .calibration.tsv (arm {a.arm}, δ = {chosen:g}, merge {a.merge})")


if __name__ == "__main__":
    sys.exit(main())
