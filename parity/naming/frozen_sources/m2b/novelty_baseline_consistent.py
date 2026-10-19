#!/usr/bin/env python3
"""Phase 4: reviewed software patch to the uploaded novelty decision tree.

This is not a newly calibrated biological model. Legacy status values are
preserved for compatibility; novelty_claim and artifact_claim make their
provisional interpretation explicit. The real v4 rerun remains pending.

Every copy gets one status and a reason code, testing the cheapest explanation first:

  1. artifact     selection signature: against the nearest main-type copy of another culture (by SSU
                  flank), the SSU differs by >= --artifact-min-ssu % and by at least --artifact-ratio of
                  the full-ITS difference. All SSU-tied nearest copies must agree. This is a heuristic
                  inspired by PX215117.1, not a diagnostic test of function or an artifact mechanism.
                  Reason A1.
  2. divergent    short-5.8S class, read from the copy itself. The comparison is the main-type 5.8S of the
                  other cultures that share the copy's best VT (or, failing that, of the nearest culture if
                  it is within 1% in SSU).
                  D1: the 5.8S is >= --short-58s bp shorter AND >= --subst-58s % different in sequence.
                  A shorter 5.8S below the sequence-distance cutoff is only flagged (D4); its boundaries
                  need inspection before a biological interpretation.
                  D2: the copy's best-hit VT or SH is a paralog unit, reached in other cultures only by
                  short-5.8S copies, from at least two cultures, and by no main-type copy.
                  Placed by the VT arm alone as a compatible placement, not proof of genomic origin.
                  A same-length 5.8S that differs by >= 3% is only flagged (D3): such copies turned out to
                  be second lineages in Septoglomus and co-dominant types in Gigaspora, so they are not an
                  explanation and go on through the tree.
  3. boundary     no VT call, but the recorded best VT identity is just below the threshold (B1).
                  This identifies a threshold-adjacent candidate; coverage still needs inspection.
  4. placement / novelty, from the Phase 3 call:
       known species / known species group   all matched arms agree and neither VT nor LSU is unmatched (K1)
       novel species, genus known            a species call demoted because VT or LSU is unmatched (N1, N2),
                                             or no species survives while a genus does (N3)
       discordant                            matched arms share no genus (X1)
       novel lineage                         candidate only: no arm matched (N4; N5 when a recorded best
                                             VT identity is below the boundary)

Scoring against what the culture reference knows (the reconcile run's holdout decides the truth groups):
  main-type copies, species held by another culture    expected known species or group (the M3 control)
  main-type copies, species held out                   expected novel species / novel lineage, not a
                                                       known species (false species)
  divergent-class copies (reference labels)            expected divergent (short 5.8S) or, for same-length
                                                       copies, anything but a false call of novelty
Writes {out}.copies.tsv with status, reason codes and placement for every copy.

Usage (from amfsplit_stage3):
  python3 novelty.py --ref ../ref/culture_ref/v4 --reconcile ../bench/stefani/reconcile_v4 \\
      --vt ../bench/stefani/vt_sets_v4 --sh ../bench/stefani/sets_SH_v4 --lsu ../bench/stefani/sets_LSU_v4 \\
      --out ../bench/stefani/novelty_v4
  (and the same with --reconcile ../bench/stefani/reconcile_v4_species --holdout species)
"""
import argparse
import collections
import datetime
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import statistics
import sys

import edlib
import arms

from arms import is_named, read_fasta, read_tsv, write_tsv

REVIEW_VERSION = "phase4-baseline-consistency-experiment-20260929.1"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dist(a, b, k=-1):
    if not a or not b:
        raise ValueError("Cannot compare empty marker sequences")
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    d = edlib.align(q, t, mode="HW", task="distance", k=k)["editDistance"]
    return None if d < 0 else 100 * d / len(q)


def split(s, sep=","):
    return [u for u in s.split(sep) if u] if s else []


def load_sets(prefix):
    """Per copy: set of units, best raw unit with its identity, and the anchor unit (the unit of the set
    that holds the best raw unit, merged or not). Per culture: culture set."""
    rows = read_tsv(f"{prefix}.copies.tsv")
    col = "unit_set" if "unit_set" in rows[0] else "vt_set"
    top = "top_units" if "top_units" in rows[0] else "top_vts"
    units = {r["accession"]: set(split(r[col])) for r in rows}
    best, anchor = {}, {}
    for r in rows:
        t = split(r[top])
        b = (t[0].rsplit(":", 1)[0], float(t[0].rsplit(":", 1)[1])) if t else ("", 0.0)
        best[r["accession"]] = b
        anchor[r["accession"]] = next((u for u in units[r["accession"]] if b[0] in u.split("+")), "")
    cul = {r["culture"]: set(split(r["culture_set"])) for r in read_tsv(f"{prefix}.cultures.tsv")}
    return units, best, anchor, cul


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--reconcile", required=True, help="prefix of a reconcile.py run")
    ap.add_argument("--vt", required=True)
    ap.add_argument("--sh")
    ap.add_argument("--lsu")
    ap.add_argument("--holdout", choices=["culture", "species"], default="culture",
                    help="must match the reconcile run")
    ap.add_argument("--artifact-min-ssu", type=float, default=1.0)
    ap.add_argument("--artifact-ratio", type=float, default=0.5)
    ap.add_argument("--short-58s", type=int, default=3)
    ap.add_argument("--subst-58s", type=float, default=3.0)
    ap.add_argument("--vt-min-id", type=float, default=97.0)
    ap.add_argument("--boundary", type=float, default=1.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if not (0 < a.artifact_min_ssu <= 100 and 0 < a.artifact_ratio
            and math.isfinite(a.artifact_ratio) and a.short_58s > 0
            and 0 < a.subst_58s <= 100 and 0 < a.vt_min_id <= 100
            and 0 <= a.boundary <= a.vt_min_id):
        ap.error("Invalid threshold parameters")
    if Path(f"{a.out}.copies.tsv").exists() or Path(f"{a.out}.run.json").exists():
        ap.error("Output already exists; use a new --out prefix to preserve earlier results")

    ref_rows = [r for r in read_tsv(os.path.join(a.ref, "copies.tsv")) if r["status"] == "retained"]
    cop = {r["accession"]: r for r in ref_rows}
    cul_meta = {r["culture"]: r for r in read_tsv(os.path.join(a.ref, "cultures.tsv"))}
    rec = {r["accession"]: r for r in read_tsv(f"{a.reconcile}.copies.tsv")}
    groups = collections.Counter(r["group"] for r in rec.values())
    if (a.holdout == "species") != (groups.get("named, species elsewhere", 0) == 0):
        raise SystemExit(f"--holdout {a.holdout} does not match the reconcile run ({dict(groups)})")
    fa = {m: read_fasta(os.path.join(a.ref, "fasta", f"{m}.fasta")) for m in ("ssu_flank", "full_its", "s58")}
    vt_units, vt_best, vt_anchor, vt_cul = load_sets(a.vt)
    arm_units, arm_anchor = {"VT": vt_units}, {"VT": vt_anchor}
    for name, p in (("SH", a.sh), ("LSU", a.lsu)):
        if p:
            arm_units[name], _, arm_anchor[name], _ = load_sets(p)
    vt_holders = collections.defaultdict(set)
    for cu, us in vt_cul.items():
        for u in us:
            vt_holders[u].add(cu)

    species = {c: m["organism"] for c, m in cul_meta.items()
               if is_named(m["organism"]) and m.get("name_status", "ok") == "ok"}
    by_species = collections.defaultdict(set)
    for c, sp in species.items():
        by_species[sp].add(c)

    def held(c):
        return by_species[species[c]] if a.holdout == "species" and c in species else {c}

    mains = [x for x, r in cop.items() if r["copy_class"] == "main" and x in fa["ssu_flank"]]
    main_58 = collections.defaultdict(list)
    for x in mains:
        if cop[x].get("s58_len") and x in fa["s58"]:
            main_58[cop[x]["culture"]].append(x)

    # anchor units reached by main-type copies vs short-5.8S copies, per culture, for the paralog-unit
    # test (VT and SH only: the paralog units found in Phase 1 are VTs and SHs)
    para_arms = [n for n in ("VT", "SH") if n in arm_anchor]
    reach = {n: (collections.defaultdict(set), collections.defaultdict(set)) for n in para_arms}
    for n in para_arms:
        for x, u in arm_anchor[n].items():
            if x not in cop or not u:
                continue
            if cop[x]["copy_class"] == "main":
                reach[n][0][u].add(cop[x]["culture"])
            elif cop[x].get("div_kind") == "short_5.8S":
                reach[n][1][u].add(cop[x]["culture"])

    out = []
    for x, r in cop.items():
        c = r["culture"]
        drop = held(c)
        reasons, flags = [], []
        # nearest main-type copy of another culture, by SSU flank (bounded edlib, same-VT candidates first)
        near, near_d, nearest = None, None, []
        if x in fa["ssu_flank"]:
            q = fa["ssu_flank"][x]
            mine = vt_units.get(x, set())
            cands = [y for y in mains if cop[y]["culture"] not in drop]
            cands.sort(key=lambda y: (not (vt_units.get(y, set()) & mine), cop[y]["culture"], y))
            for y in cands:
                k = -1 if near_d is None else int(near_d * min(len(q), len(fa["ssu_flank"][y])) / 100) + 1
                d = dist(q, fa["ssu_flank"][y], k)
                if d is None:
                    continue
                if near_d is None or d < near_d - 1e-12:
                    near_d, nearest = d, [y]
                elif math.isclose(d, near_d, rel_tol=0, abs_tol=1e-12):
                    nearest.append(y)
            nearest.sort(key=lambda y: (cop[y]["culture"], y))
            near = nearest[0] if nearest else None
        its_d = dist(fa["full_its"][x], fa["full_its"][near]) if near and x in fa["full_its"] and near in fa["full_its"] else None

        # 1. artifact
        # SSU often has ties. Require complete, concordant evidence across ALL
        # tied nearest copies before the A1 screen fires. The representative
        # pair remains in the legacy report columns; it is not the sole test.
        artifact = False
        artifact_its = []
        if near_d is not None and near_d >= a.artifact_min_ssu:
            artifact_its = [dist(fa["full_its"][x], fa["full_its"][y])
                            for y in nearest if x in fa["full_its"] and y in fa["full_its"]]
            votes = [near_d >= a.artifact_ratio * d for d in artifact_its]
            artifact = len(votes) == len(nearest) and bool(votes) and all(votes)
            if votes and any(votes) and (not all(votes) or len(votes) != len(nearest)):
                flags.append("A2 artifact screen unresolved across SSU-tied nearest copies")
        if artifact:
            reasons.append(f"A1 screening flag: SSU {near_d:.1f}% vs ITS "
                           f"{min(artifact_its):.1f}-{max(artifact_its):.1f}% across "
                           f"{len(nearest)} SSU-tied nearest copies; not proof of an artifact")

        # 2. divergent class from the copy's own sequence
        short = paralog = False
        ncul = cop[near]["culture"] if near else ""
        comp = sorted(cu for cu in vt_holders.get(vt_anchor.get(x, ""), set()) - drop if main_58.get(cu))
        comp_src = "cultures sharing its VT"
        if not comp and near_d is not None and near_d <= 1.0:
            comp = sorted({cop[y]["culture"] for y in nearest if main_58.get(cop[y]["culture"])})
            comp_src = f"SSU-tied nearest cultures ({near_d:.2f}% SSU)"
        s58_info = ""
        if r.get("s58_len") and comp and x in fa["s58"]:
            ref_len = statistics.median(statistics.median(int(cop[y]["s58_len"]) for y in main_58[cu]) for cu in comp)
            delta = int(r["s58_len"]) - ref_len
            # Development ablation: when testing a short query against ref_len,
            # the sequence comparator must represent that length baseline too.
            # Do not let a shorter main-labelled sequence supply a different
            # baseline for the sequence-distance half of the SAME D1 test.
            compare_ids = [y for cu in comp for y in main_58[cu]]
            if delta <= -a.short_58s:
                compare_ids = [y for y in compare_ids if len(fa["s58"][y]) >= ref_len]
            d58 = min(dist(fa["s58"][x], fa["s58"][y]) for y in compare_ids)
            s58_info = f"5.8S {r['s58_len']} bp vs {ref_len:g} in {comp_src}, {d58:.1f}% different"
            if delta <= -a.short_58s and d58 >= a.subst_58s:
                short = True
                reasons.append(f"D1 {s58_info}")
            elif delta <= -a.short_58s:
                flags.append(f"D4 {s58_info} (shorter, below divergence threshold: inspect boundaries)")
            elif abs(delta) < a.short_58s and d58 >= a.subst_58s:
                flags.append(f"D3 {s58_info}")
        for n in para_arms:
            u = arm_anchor[n].get(x, "")
            if not u:
                continue
            m_c = reach[n][0].get(u, set()) - drop
            d_c = reach[n][1].get(u, set()) - drop
            if len(d_c) >= 2 and not m_c:
                paralog = True
                reasons.append(f"D2 paralog unit {u} ({n}; short-5.8S copies of {len(d_c)} other cultures)")
        divergent = short or paralog

        # 3. boundary
        rc = rec[x]
        vt_state = rc["VT_state"]
        best_vt, best_id = vt_best.get(x, ("", 0.0))
        boundary = bool(best_vt) and vt_state == "no call" and a.vt_min_id - a.boundary <= best_id < a.vt_min_id
        off_vt = bool(best_vt) and vt_state == "no call" and best_id < a.vt_min_id - a.boundary
        if vt_state == "no call" and best_vt and best_id >= a.vt_min_id:
            flags.append("B2 VT identity passes threshold but no call: inspect coverage, filtering and provenance")
        if vt_state == "no call" and not best_vt:
            flags.append("Q1 no recorded best VT hit; distance from all VTs is not established")

        # 4. placement and novelty
        rank, um = rc["rank"], set(split(rc["unmatched_arms"], ";"))
        if artifact:
            status = "artifact"
        elif divergent:
            status = "divergent class"
            place = rc["VT_species"] or rc["VT_genera"]
            reasons.append(f"placed by VT: {place or 'no VT match'}")
        elif boundary:
            status = "boundary"
            reasons.append(f"B1 best VT {best_vt} {best_id:.2f}%")
        elif rank in ("species", "species group"):
            if um & {"VT", "LSU"}:
                status = "novel species, genus known"
                reasons += [f"{'N1' if n == 'VT' else 'N2'} {n} unmatched" for n in ("VT", "LSU") if n in um]
            else:
                status = "known " + rank
                reasons.append("K1")
        elif rank == "genus":
            status = "novel species, genus known"
            reasons.append("N3 no species shared by the arms")
        elif rank == "discordant":
            status = "discordant"
            reasons.append("X1 arms share no genus")
        else:
            status = "novel lineage"
            reasons.append("N5 recorded best VT identity below boundary" if off_vt else "N4 no arm matched")
        placement = (rc["call_species"] if status.startswith("known") else
                     (rc["VT_species"] or rc["VT_genera"]) if status == "divergent class" else
                     rc["call_genera"] if status == "novel species, genus known" else "")
        out.append({"accession": x, "culture": c, "organism": r["organism"], "copy_class": r["copy_class"],
                    "div_kind": r.get("div_kind", ""), "group": rc["group"], "status": status,
                    "novelty_claim": "candidate_only" if status.startswith("novel") else "not_asserted",
                    "artifact_claim": "screen_flag_only" if status == "artifact" else "not_asserted",
                    "placement": placement, "reasons": "; ".join(reasons), "flags": "; ".join(flags),
                    "s58": s58_info, "nearest_culture": ncul, "nearest_ssu_pct": round(near_d, 2) if near_d is not None else "",
                    "nearest_its_pct": round(its_d, 2) if its_d is not None else "",
                    "nearest_ssu_copy_count": len(nearest),
                    "nearest_ssu_cultures": ";".join(sorted({cop[y]["culture"] for y in nearest})),
                    "artifact_comparator_its_pct": ";".join(f"{d:.8g}" for d in sorted(artifact_its)),
                    "best_vt": best_vt, "best_vt_pident": best_id, "reconcile_rank": rank,
                    "unmatched_arms": rc["unmatched_arms"]})

    # ---- report
    def truth(r):
        if r["copy_class"] == "divergent":
            return f"divergent ({r['div_kind'] or '?'})"
        if r["copy_class"] != "main":
            return r["copy_class"]
        return {"named, species elsewhere": "main, species known", "named, species held out": "main, species held out",
                "unnamed ('sp.')": "main, 'sp.' culture"}[r["group"]]
    statuses = ["known species", "known species group", "novel species, genus known", "novel lineage",
                "boundary", "divergent class", "discordant", "artifact"]
    heads = ["known sp.", "known grp", "novel sp.", "novel lin.", "boundary", "divergent", "discord.", "artifact"]
    print(f"{len(out)} copies; holdout {a.holdout}; arms {', '.join(arm_units)}")
    print(f"\n1. status by what the reference knows (rows) ")
    print(f"  {'':30s}{'copies':>7s}" + "".join(f"{h:>11s}" for h in heads))
    rows_order = ["main, species known", "main, species held out", "main, 'sp.' culture",
                  "divergent (short_5.8S)", "divergent (same_length)", "divergent (longer)", "divergent (?)", "outlier"]
    for t in rows_order:
        rs = [r for r in out if truth(r) == t]
        if rs:
            cnt = collections.Counter(r["status"] for r in rs)
            print(f"  {t:30s}{len(rs):>7d}" + "".join(f"{cnt[s]:>11d}" for s in statuses))

    known = [r for r in out if truth(r) == "main, species known"]
    heldo = [r for r in out if truth(r) == "main, species held out"]
    short = [r for r in out if truth(r) == "divergent (short_5.8S)"]
    print("\n2. checks")
    if known:
        bad = [r for r in known if not r["status"].startswith("known")]
        print(f"  M3 negative control: main-type copies of known species not called known: {len(bad)} of {len(known)}")
        for cu, n in collections.Counter(r["culture"] for r in bad).most_common(10):
            r0 = next(r for r in bad if r["culture"] == cu)
            print(f"      {cu} ({n}): {r0['status']} | {r0['reasons']}"[:220])
    if heldo:
        fs = [r for r in heldo if r["status"].startswith("known")]
        print(f"  false species: held-out species called a known species or group: {len(fs)} of {len(heldo)}"
              f" ({100 * len(fs) / len(heldo):.0f}%)")
        for cu, n in collections.Counter(r["culture"] for r in fs).most_common(8):
            r0 = next(r for r in fs if r["culture"] == cu)
            print(f"      {cu} ({n}): {r0['status']} {r0['placement']}"[:200])
    if short:
        det = sum(r["status"] == "divergent class" for r in short)
        how = collections.Counter(code for r in short if r["status"] == "divergent class"
                                  for code in ("D1", "D2") if code in r["reasons"])
        print(f"  short-5.8S copies recognised from their own sequence: {det} of {len(short)} "
              f"(5.8S {how['D1']}, paralog unit {how['D2']}); called novel: "
              f"{sum(r['status'].startswith('novel') for r in short)}"
              + (" (species holdout tests reference absence; copy type and novelty remain separate questions)"
                 if a.holdout == "species" else ""))
        for r in [r for r in short if r["status"] != "divergent class"][:10]:
            print(f"      not recognised: {r['accession']} {r['culture']}: {r['s58'] or 'no 5.8S comparison'}; "
                  f"{r['status']}"[:220])
    fp = [r for r in out if r["copy_class"] == "main" and r["status"] == "divergent class"]
    print(f"  main-type copies called divergent: {len(fp)}")
    for r in fp[:6]:
        print(f"      {r['accession']} {r['culture']}: {r['reasons']}"[:200])
    art = [r for r in out if r["status"] == "artifact"]
    print(f"  artifacts: {len(art)}" + (" (the culture reference should have none left)" if not art else ""))
    for r in art[:6]:
        print(f"      {r['accession']} {r['culture']}: {r['reasons']}"[:200])
    d3 = [r for r in out if "D3" in r["flags"]]
    print(f"  same-length 5.8S flags (D3): {len(d3)} copies "
          f"({sum(r['copy_class'] == 'divergent' for r in d3)} of them divergent-class in the reference)")
    d4 = [r for r in out if "D4" in r["flags"]]
    print(f"  shorter 5.8S below the divergence threshold (D4, inspect boundaries): {len(d4)} copies "
          f"({sum(r['copy_class'] == 'main' for r in d4)} main-type)")
    for r in d4[:6]:
        print(f"      {r['accession']} {r['culture']} ({r['copy_class']}): {r['s58']}"[:200])
    ratio = [r["nearest_ssu_pct"] / r["nearest_its_pct"] for r in out
             if r["copy_class"] == "main" and r["nearest_its_pct"] not in ("", 0, 0.0) and r["nearest_ssu_pct"] != ""]
    if len(ratio) >= 2:
        q = statistics.quantiles(ratio, n=100)
        print(f"  SSU/ITS divergence ratio to the nearest culture, main-type copies: median {statistics.median(ratio):.3f}, "
              f"95th {q[94]:.3f}, max {max(ratio):.3f} (artifact needs >= {a.artifact_ratio} with SSU >= "
              f"{a.artifact_min_ssu}%)")

    if not out:
        raise ValueError("No retained copies to report")
    write_tsv(f"{a.out}.copies.tsv", out, list(out[0]))
    inputs = [Path(a.ref)/n for n in ("copies.tsv", "cultures.tsv")]
    inputs += [Path(a.ref)/"fasta"/f"{m}.fasta" for m in ("ssu_flank", "full_its", "s58")]
    inputs.append(Path(f"{a.reconcile}.copies.tsv"))
    for prefix in (a.vt, a.sh, a.lsu):
        if prefix:
            inputs += [Path(f"{prefix}.{suffix}.tsv") for suffix in ("copies", "cultures")]
    manifest = {
        "review_version": REVIEW_VERSION,
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256(__file__),
        "python_version": sys.version,
        "edlib_version": importlib.metadata.version("edlib"),
        "arms_module_sha256": sha256(arms.__file__) if getattr(arms, "__file__", None) else None,
        "parameters": vars(a),
        "inputs": [{"path": str(p.resolve()), "sha256": sha256(p)} for p in inputs],
        "output_sha256": sha256(f"{a.out}.copies.tsv"),
        "interpretation": "Legacy novelty statuses denote candidates; artifact is a screening flag. Independent-library and upstream artifact checks are not implemented.",
        "validation": "Development ablation on already examined v4 data. Short-query sequence comparators are restricted to length >= the median length baseline. Thresholds unchanged. Global arm merges are not refitted inside this script.",
    }
    Path(f"{a.out}.run.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print("  Interpretation: novelty statuses are candidates; artifact is a screening flag, not a confirmed diagnosis.")
    print(f"\ntable: {a.out}.copies.tsv")


if __name__ == "__main__":
    sys.exit(main())
